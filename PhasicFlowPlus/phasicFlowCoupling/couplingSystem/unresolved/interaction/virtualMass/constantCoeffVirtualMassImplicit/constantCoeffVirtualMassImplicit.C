/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): Antonio Di Mauro
  email: antoniodimauro03@gmail.com
------------------------------------------------------------------------------
Licence:
  This file is part of PhasicFlow, a CFD-DEM stack built on phasicFlow and
  PhasicFlowPlus (www.cemf.ir). It is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License v3 or
  any later version.

  It is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.
-----------------------------------------------------------------------------*/

// Maxey & Riley (1983) Phys. Fluids 26, 883

#include <algorithm>

#include "fvc.H"

#include "constantCoeffVirtualMassImplicit.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "momentumSemiImplicitSphereUnresolvedCouplingSystem.hpp"
#include "distributionBase.hpp"
#include "fluidAcceleration.hpp"

pFlow::coupling::constantCoeffVirtualMassImplicit::constantCoeffVirtualMassImplicit
(
    const unresolvedCouplingSystem& uCS,
    const porosity& prsty
)
:
    virtualMassImplicit(uCS, prsty),
    Cvm_(this->dict().getOrDefault<Foam::scalar>("Cvm", 0.5)),
    isCompressible_
    (
      this->mesh().template lookupObject<Foam::volScalarField>("p").dimensions() 
      == Foam::dimPressure
    )

    
{
    tmpVirtualMassForce_ = Foam::tmp<Foam::volVectorField>::New
    (
        Foam::IOobject
        (
            "virtualMassForce",
            Foam::timeName(this->mesh().time()),
            this->mesh(),
            Foam::IOobject::READ_IF_PRESENT,
            (this->printVirtualMass()?Foam::IOobject::AUTO_WRITE:Foam::IOobject::NO_WRITE)
        ),
        this->mesh(),
        Foam::dimensionedVector
        (
            "virtualMassForce",
            Foam::dimensionSet(1,-2,-2,0,0),
            Foam::vector(0,0,0)
        )
    );
}

void pFlow::coupling::constantCoeffVirtualMassImplicit::calculateVirtualMassForceImplicit
(
    const Foam::volVectorField& U,
    const Foam::volVectorField& DDtUField,
    sphereQuadrature* quadPtr,
    const Plus::realProcCMField& diameter,
    const distributionBase& cellDistribution,
    Plus::realx3ProcCMField& particleForce,
    Plus::realProcCMField& addedMass
)
{
    auto& vMF = tmpVirtualMassForce_.ref();

    forAll(vMF, celli)
    {
        vMF[celli] = Foam::vector(0,0,0);
    }

    const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
    const auto& parCellInd = this->parCellIndex();
    const auto& Vcells = this->mesh().V();
    
    const auto& tDDtU = DDtUField;
    
    const auto& DDtU = this->sampleDDtU(tDDtU);
    
    const size_t nPar = diameter.size();
    Plus::realx3ProcCMField& vmPar = this->vmScratch();
    std::fill(vmPar.begin(), vmPar.end(), realx3(0.0));

    Foam::scalar amFaxenRatio = 0.0;

    // Faxen volume average
    if(quadPtr) quadPtr->bind(tDDtU);

    #pragma omp parallel for schedule(dynamic) reduction(max:amFaxenRatio)
    for(size_t i=0; i<nPar; ++i)
    {
        const Foam::label cellI = parCellInd[i];

        if(cellI < 0) continue;

        const Foam::scalar dp = diameter[i];
        const Foam::scalar Vp = Foam::constant::mathematical::pi/6.0 * Foam::pow(dp,3.0);
        
        const Foam::scalar ma = Cvm_ * rho[cellI] * Vp;
               
        addedMass[i] += static_cast<real>(ma);
        
        Foam::vector accEff = DDtU[i];

        if(quadPtr)
        {
            const Foam::vector aV = quadPtr->volumeAverage(i, tDDtU);
            amFaxenRatio = Foam::max
            (
                amFaxenRatio,
                Foam::mag(aV - accEff)/Foam::max(Foam::mag(aV), SMALL)
            );
            accEff = aV;
        }
        
        const Foam::vector vmForce = ma * accEff;
    
        particleForce[i] += realx3(vmForce.x(), vmForce.y(), vmForce.z());

        vmPar[i] = realx3(vmForce.x(), vmForce.y(), vmForce.z());
    }

    if(quadPtr)
    {
        quadPtr->unbind();
        this->reportAddedMassFaxen(amFaxenRatio);
    }

    semiImplicitCoupling(this->Porosity().uCS()).globalParticleSum(vmPar);

    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<nPar; ++i)
    {
        const Foam::vector vmForce(vmPar[i].x(), vmPar[i].y(), vmPar[i].z());

        if(Foam::magSqr(vmForce) == 0) continue;
        if(!cellDistribution.hasSupport(i) && parCellInd[i] < 0) continue;

        cellDistribution.distributeValue_OMP(i, parCellInd[i], vMF, vmForce);
    }

    forAll(Vcells, celli)
    {
        vMF[celli] /= Vcells[celli];
    }

    vMF.correctBoundaryConditions();
}
