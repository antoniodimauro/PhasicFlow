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


#include "surfTorqueShi2019Implicit.hpp"
#include "wallDist.H"
#include "fvc.H"
#include "lift.hpp"
#include "rotationalDrag.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "pairInteraction.hpp"

pFlow::coupling::surfTorqueShi2019Implicit::surfTorqueShi2019Implicit
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    surfaceRotationTorqueImplicit(uCS, prsty),
    residualRe_(this->dict().get<Foam::scalar>("residualRe"))
{
    const Foam::word rot = this->dict().template getOrDefault<Foam::word>
    (
        "rotationLift", "none"
    );

    if(rot == "Bluemink2010")
    {
        rotationSpin_  = true;
        strainSampler_ = makeUnique<fluidFieldSampler>(uCS, "torqueStrainRate");

        Foam::Info
            << "    Spin in rotation: the rotation part (|omega| - sigma) of the"
               " vorticity drives the spin\n"
               "      to (1 + 0.0045 Re_p) times the fluid rotation (Bluemink"
               " et al. 2010 eq. 2.13);\n"
               "      the shear part keeps the Shi & Rzehak (2019) ratio"
               " f_shear/f_spin.\n";
    }
}


void pFlow::coupling::surfTorqueShi2019Implicit::calculateSurfaceTorqueImplicit
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
    const auto& Uav   = this->sampleFluidVel(U);

    // strain rate sqrt(2 S:S), in the x component, for the shear/rotation split
    const Plus::procCMField<Foam::vector>* strainPtr = nullptr;
    Foam::tmp<Foam::volVectorField> tStrain;
    if(rotationSpin_)
    {
        tStrain = Foam::sqrt(2.0)*Foam::mag(Foam::symm(Foam::fvc::grad(U)))
                * Foam::dimensionedVector("ex", Foam::dimless, Foam::vector(1, 0, 0));
        strainPtr = &strainSampler_().sample(tStrain());
    }

    const auto& parCentre = this->Porosity().uCS().centerMass();
    const Foam::volVectorField& cellCentre = this->mesh().C();

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

        const Foam::vector uf = Uav[i];
        const Foam::vector up{
            parVel[i].x(), 
            parVel[i].y(), 
            parVel[i].z()};

        const Foam::vector wp{
            parRotVel[i].x(),
            parRotVel[i].y(),
            parRotVel[i].z()
        };

        // pairVorticity: the neighbours' vorticity the kernel smears out
        Foam::vector wf = curlU[i];
        if(this->pair() && this->pair()->pairVorticity())
        {
            wf += this->pair()->dVort(i);
        }

        const Foam::scalar dp = diameter[i];

        const Foam::vector uRelVec = up - uf;
        const Foam::scalar uRel = Foam::mag(uRelVec);

        const Foam::scalar Rep = Foam::max(uRel*dp/nu[cellI], residualRe_);
        const Foam::scalar Rew = Foam::max(Foam::mag(wp)*dp*dp/nu[cellI], Foam::SMALL);
        const Foam::scalar Res = Foam::max(Foam::mag(wf)*dp*dp/nu[cellI], Foam::SMALL);
       
        
        const Foam::scalar f_spin = rotationalDragFactor(Rew);
        const Foam::scalar f_shear = f_spin * 
            (
                (1 + 0.4 * (Foam::exp(-0.0135*Res) - 1)) * 
                (1-0.07026*Foam::pow(Rep, 0.455))
            );
            
        const Foam::scalar piRhoNuD3 =
            Foam::constant::mathematical::pi * rho[cellI] * nu[cellI] * dp * dp * dp;

        const Foam::scalar K_omega = piRhoNuD3 * f_spin;

        // spin at which the fluid torque vanishes
        Foam::vector spinTarget = (0.5*f_shear/f_spin)*wf;
        if(rotationSpin_)
        {
            const Foam::scalar wMag  = Foam::mag(wf);
            const Foam::scalar sigma = (*strainPtr)[i].x();
            if(wMag > Foam::SMALL && sigma < wMag)
            {
                const Foam::vector wShear = (sigma/wMag)*wf;
                const Foam::vector wRot   = wf - wShear;
                // Bluemink et al. (2010) eq. (2.13), fitted for Re_p <= 200
                const Foam::scalar spinRatio =
                    1.0 + 0.0045*Foam::min(Rep, Foam::scalar(200));
                spinTarget = (0.5*f_shear/f_spin)*wShear + 0.5*spinRatio*wRot;
            }
        }

        const Foam::vector torque_shear = K_omega*spinTarget;

        Foam::vector torque_total = torque_shear;

        if(this->wallCorrection())
        {
            Foam::vector nw = tNw()[cellI];
            const Foam::scalar nMag = Foam::mag(nw);

            if(nMag > Foam::SMALL)
            {
                nw /= nMag;

                const Foam::scalar a = 0.5*dp;
                // wall distance at the particle centre, as in the drag
                const Foam::vector xp
                (
                    parCentre[i].x(), parCentre[i].y(), parCentre[i].z()
                );
                const Foam::scalar h = Foam::max
                (
                    tYw()[cellI] + (nw & (xp - cellCentre[cellI])),
                    Foam::scalar(0)
                );
                const Foam::scalar hOverA = Foam::max
                (
                    h/Foam::max(a, Foam::SMALL),
                    this->minWallRadii()
                );

                torque_total += wallCouple
                (
                    hOverA, nw, wp, wf, up,
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
            spinTarget,
            0.25*dp*dp/nu[cellI],
            piRhoNuD3/3.0,
            K_omega
        );

       
    }

}