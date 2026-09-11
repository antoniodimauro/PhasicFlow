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


#include "lowReynoldsImplicit.hpp"
#include "wallDist.H"
#include "fvc.H"
#include "lift.hpp"

pFlow::coupling::lowReynoldsImplicit::lowReynoldsImplicit
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    surfaceRotationTorqueImplicit(uCS, prsty)
{

}


void pFlow::coupling::lowReynoldsImplicit::calculateSurfaceTorqueImplicit
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  parVel,
    const Plus::realx3ProcCMField&  parRotVel,
    const Plus::realProcCMField&    diameter,
    Plus::realx3ProcCMField&        particleTorque,
    Plus::realProcCMField&          rotDragCoeff
)
{

    const size_t nPar = diameter.size();
    const auto& parCellInd =  this->parCellIndex();
    const auto& nu = this->mesh().template lookupObject<Foam::volScalarField>("nu");
    const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
    
    auto curlUPtr = lift::fluidVorticity(U);
    const auto& curlU = this->sampleVorticity(curlUPtr.ref());

    Foam::tmp<Foam::volScalarField> tYw;
    Foam::tmp<Foam::volVectorField> tNw;

    if(this->wallCorrection())
    {
        tYw = Foam::wallDist::New(this->mesh()).y();
        tNw = Foam::fvc::grad(tYw());
    }

    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0; i<nPar; ++i)
    {
        const Foam::label cellI = parCellInd[i];
        
        if(cellI == -1 )continue;

        const Foam::vector up{
            parVel[i].x(), 
            parVel[i].y(), 
            parVel[i].z()};

        const Foam::vector wp{
            parRotVel[i].x(),
            parRotVel[i].y(),
            parRotVel[i].z()
        };

        const Foam::scalar dp = diameter[i];
        
        const Foam::scalar f_spin = 1.0;
        const Foam::scalar f_shear = 1.0;
            
        const Foam::scalar piRhoNuD3 =
            Foam::constant::mathematical::pi * rho[cellI] * nu[cellI] * dp * dp * dp;

        const Foam::scalar K_omega = piRhoNuD3 * f_spin;
        const Foam::vector torque_shear =
            0.5 * piRhoNuD3 * f_shear * curlU[i];

        Foam::vector torque_total = torque_shear;

        if(this->wallCorrection())
        {
            Foam::vector nw = tNw()[cellI];
            const Foam::scalar nMag = Foam::mag(nw);

            if(nMag > Foam::SMALL)
            {
                nw /= nMag;

                const Foam::scalar a = 0.5*dp;
                const Foam::scalar hOverA = Foam::max
                (
                    tYw()[cellI]/Foam::max(a, Foam::SMALL),
                    this->minWallRadii()
                );

                torque_total += wallCouple
                (
                    hOverA, nw, wp, curlU[i], up,
                    K_omega, 0.5*piRhoNuD3*f_shear, dp
                );
            }
        }

        particleTorque[i] += realx3
        (
            torque_total.x(),
            torque_total.y(),
            torque_total.z()
        );

        rotDragCoeff[i] += static_cast<real>(K_omega);

        this->history().publish
        (
            i,
            (0.5*f_shear/f_spin)*curlU[i],
            0.25*dp*dp/nu[cellI],
            piRhoNuD3/3.0,
            K_omega
        );

       
    }

}