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

// Shi & Rzehak (2019) Chem. Eng. Sci. 208, 115145

#include <algorithm>

#include "Shi2019Implicit.hpp"
#include "shearLiftJ.hpp"
#include "wallDist.H"
#include "unresolvedCouplingSystem.hpp"
#include "momentumSemiImplicitSphereUnresolvedCouplingSystem.hpp"
#include "distributionBase.hpp"

pFlow::coupling::Shi2019Implicit::Shi2019Implicit
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    liftImplicit(uCS, prsty),
    residualRe_
    (
        this->dict().template get<Foam::scalar>("residualRe")
    ),
    rotationLift_(false)
{
    this->setLiftActive(true);

    const Foam::word rot = this->dict().template getOrDefault<Foam::word>
    (
        "rotationLift", "none"
    );

    if(rot == "Bluemink2010")
    {
        rotationLift_ = true;
        strainSampler_ = makeUnique<fluidFieldSampler>(uCS, "strainRate");
        rotationLiftEkman_ = this->dict().template getOrDefault<Foam::Switch>
        (
            "rotationLiftEkman", false
        );
    }
    else if(rot != "none")
    {
        FatalErrorInFunction
            << "rotationLift " << rot << " unknown; valid: none, Bluemink2010"
            << Foam::exit(Foam::FatalError);
    }

    Foam::Info
        << "    Rotation lift: " << Green_Text(rot) << "\n";
    if(rotationLiftEkman_)
    {
        Foam::Info
            << "      rotationLiftEkman: the fit at Re_p + 2 Ta, Ta = Omega a^2/nu,"
               " Omega = (|omega| - sigma)/2\n"
               "      (the Ekman layer thins the boundary layer as inertia"
               " does).\n";
    }

    {
        const Foam::word sl = this->dict().template getOrDefault<Foam::word>
        (
            "spinLiftHighRotation", "none"
        );
        if(sl == "OesterleDinh1998")
        {
            spinLiftOesterleDinh_ = true;
            Foam::Info
                << "    Spin lift at high Rr: " << Green_Text("OesterleDinh1998")
                << ", C = 0.45 + (Rr - 0.45) exp(-0.05684 Rr^0.4 Re_p^0.7)\n"
                   "      (Oesterle & Dinh 1998, fitted on Rr 2-12; Rubinow & Keller"
                   " as Re_p -> 0);\n      Shi & Rzehak (2019) below Rr = 2, smoothstep"
                   " blend of the two over Rr 2-5.\n";
        }
        else if(sl != "none")
        {
            FatalErrorInFunction
                << "spinLiftHighRotation " << sl << " unknown; valid: none, OesterleDinh1998"
                << Foam::exit(Foam::FatalError);
        }
    }

    tmpLiftForce_ = Foam::tmp<Foam::volVectorField>::New
    (
        Foam::IOobject
        (
            "liftForce",
            Foam::timeName(this->mesh().time()),
            this->mesh(),
            Foam::IOobject::READ_IF_PRESENT,
            (this->printLift()?Foam::IOobject::AUTO_WRITE:Foam::IOobject::NO_WRITE)
        ),
        this->mesh(),
        Foam::dimensionedVector
        (
            "liftForce",
            Foam::dimensionSet(1,-2,-2,0,0),
            Foam::vector(0,0,0)
        )
    );
}


void pFlow::coupling::Shi2019Implicit::calculateLiftForceImplicit
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  parVel,
    const Plus::realx3ProcCMField&  parRotVel,
    const Plus::realProcCMField&    diameter,
    const distributionBase&         cellDistribution,
    Plus::realx3ProcCMField&        particleForce,
    Plus::realx3ProcCMField&        particleTorque,
    Plus::realx3ProcCMField&        liftCrossVec
) 
{
    auto& liftForce = tmpLiftForce_.ref();

    // initialize lift to zero
	forAll(liftForce, celli)
	{
		liftForce[celli] = Foam::vector(0,0,0);
	}

    const size_t nPar = diameter.size();
    const auto& parCellInd =  this->parCellIndex();
    Plus::realx3ProcCMField& lfPar = this->liftScratch();
    std::fill(lfPar.begin(), lfPar.end(), realx3(0.0));
    const auto& nu = this->mesh().template lookupObject<Foam::volScalarField>("nu");
    const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
    
    auto curlUPtr = this->fluidVorticity(U);
    const auto& curlU = this->sampleVorticity(curlUPtr.ref());
    const auto& Uav   = this->sampleFluidVel(U);

    // strain rate sqrt(2 S:S), in the x component, for the shear/rotation split
    const Plus::procCMField<Foam::vector>* strainPtr = nullptr;
    Foam::tmp<Foam::volVectorField> tStrain;
    if(rotationLift_)
    {
        tStrain = Foam::sqrt(2.0)*Foam::mag(Foam::symm(Foam::fvc::grad(U)))
                * Foam::dimensionedVector("ex", Foam::dimless, Foam::vector(1, 0, 0));
        strainPtr = &strainSampler_().sample(tStrain());
    }

    Foam::tmp<Foam::volScalarField> tYw;
    Foam::tmp<Foam::volVectorField> tNw;

    if(this->wallCorrection())
    {
        tYw = Foam::wallDist::New(this->mesh()).y();
        tNw = Foam::fvc::grad(tYw());
    }

    Foam::scalar wallLiftObs = 0.0;

    const auto& parCentre = this->Porosity().uCS().centerMass();
    const Foam::volVectorField& cellCentre = this->mesh().C();

    #pragma omp parallel for schedule(dynamic) reduction(max:wallLiftObs)
    for(size_t i=0; i<nPar; ++i)
    {
        const Foam::label cellI = parCellInd[i];
        
        if(cellI == -1 )continue;

        const Foam::vector uf = this->fluidVelForLift(i, Uav[i]);
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

        const Foam::vector uRelVec = up - uf;
        const Foam::scalar uRel = Foam::mag(uRelVec);

        const Foam::scalar Rep = Foam::max(uRel*dp/nu[cellI], residualRe_);
        const Foam::scalar Rew = Foam::max(Foam::mag(wp)*dp*dp/nu[cellI], Foam::SMALL);
        // Bluemink et al. (2010): a linear flow as a simple shear of rate
        // sigma plus a solid-body rotation of vorticity |omega| - sigma
        const Foam::scalar wMag = Foam::mag(curlU[i]);
        Foam::scalar wShear = wMag;
        Foam::scalar wRot = 0.0;
        if(rotationLift_)
        {
            const Foam::scalar sigma = (*strainPtr)[i].x();
            if(sigma < wMag)
            {
                wShear = sigma;
                wRot = wMag - sigma;
            }
        }

        const Foam::scalar Res = Foam::max(wShear*dp*dp/nu[cellI], Foam::SMALL);
        const Foam::scalar Sr = Res/Rep;
        const Foam::scalar Rr = Rew/Rep;
        
        // Shi & Rzehak (2019) eq. (13), proposed for 0.1 <= Rr <= 10
        Foam::scalar Cl_spin = Rr *
        (
            1 - 
            0.62 * Foam::tanh(0.3*Foam::sqrt(Rep)) -
            0.24 * Foam::tanh(0.01*Rep) / Foam::tanh(0.8*Foam::sqrt(Rr)) * Foam::atan(0.47*(Rr-1))
        );

        // with spinLiftHighRotation OesterleDinh1998: Oesterle & Dinh (1998),
        // fitted on 2 <= Rr <= 12, blended in over 2 <= Rr <= 5, where both
        // correlations hold, with the smoothstep weight w = t^2 (3 - 2t),
        // t = (Rr - 2)/3: the lift and its slope stay continuous
        if(spinLiftOesterleDinh_ && Rr > 2.0)
        {
            const Foam::scalar Cl_OD =
                0.45 + (Rr - 0.45)*Foam::exp(-0.05684*Foam::pow(Rr, 0.4)*Foam::pow(Rep, 0.7));
            const Foam::scalar t = Foam::min((Rr - 2.0)/3.0, 1.0);
            const Foam::scalar w = t*t*(3.0 - 2.0*t);
            Cl_spin = (1.0 - w)*Cl_spin + w*Cl_OD;
        }
        
        // Shi & Rzehak (2019) eqs. (23), (27)
        Foam::scalar Cl_shear;
        if(Rep <= 50.0)
        {
            const Foam::scalar eps = Foam::sqrt(Res)/Rep;

            const Foam::scalar J = shearLiftJ(eps);

            const Foam::scalar K = 18.0/Foam::sqr(Foam::constant::mathematical::pi);
            Cl_shear = K*
                Foam::sqrt(Sr/Rep)*J - 11.0/8.0*Sr*Foam::exp(-0.5*Rep);
        }
        else
        {
            Cl_shear = -0.064 * Foam::exp(0.525*Sr) * 
            ( 0.49 + 0.51 * Foam::tanh
                (
                    5*Foam::log10(Rep * Foam::pow(Sr,0.08)/120.0)
                )
            );
        }
        
        // Bluemink et al. (2010) eq. (2.8), non-spinning sphere, Auton form
        // F = rho V C_L (u - v) x curl(u): C_L(pi/8 form) = 4/3 Sr C_L,B
        if(wRot > 0)
        {
            const Foam::scalar SrRot = wRot*dp*dp/nu[cellI]/Rep;
            // Re_a + Ta with Ta = (wRot/2) (dp/2)^2/nu, in diameter units
            const Foam::scalar ReB = rotationLiftEkman_
                ? Rep + 0.25*wRot*dp*dp/nu[cellI]
                : Rep;
            const Foam::scalar ClB = Foam::max
            (
                0.51*Foam::log10(Foam::min(ReB, 200.0)) - 0.22, 0.0
            );
            Cl_shear += 4.0/3.0*SrRot*ClB;
        }

        const Foam::vector crossShear = curlU[i] ^ uRelVec;
        const Foam::vector crossSpin  = wp ^ uRelVec;

        const Foam::scalar magShear = Foam::max(Foam::mag(crossShear), Foam::SMALL);
        const Foam::scalar magSpin  = Foam::max(Foam::mag(crossSpin),  Foam::SMALL);

        const Foam::scalar piD2Rho = Foam::constant::mathematical::pi
                                   / 8.0 * dp * dp * rho[cellI];
        const Foam::scalar uRel2 = uRel * uRel;
            
        const Foam::scalar alphaShear = piD2Rho * Cl_shear * uRel2 / magShear;
        const Foam::scalar alphaSpin  = piD2Rho * Cl_spin  * uRel2 / magSpin;
        
        const Foam::scalar sinFloor = 0.1;

        const Foam::scalar magShearCond = Foam::max
            ( 
            magShear, sinFloor*Foam::mag(curlU[i])*uRel
            );
            
        const Foam::scalar magSpinCond = Foam::max
        (
            magSpin, sinFloor*Foam::mag(wp)*uRel
        );
        
        const Foam::scalar alphaShearImp =
            piD2Rho * Cl_shear * uRel2 / Foam::max(magShearCond, SMALL);

        const Foam::scalar alphaSpinImp =
            piD2Rho * Cl_spin  * uRel2 / Foam::max(magSpinCond, SMALL);
        
        const Foam::vector vecShear = -alphaShearImp * curlU[i];
        const Foam::vector vecSpin  = -alphaSpinImp  * wp;

        const Foam::vector liftCrossVecLocal = vecShear + vecSpin;

        const Foam::vector liftResidual =
            (alphaShear - alphaShearImp)*crossShear
          + (alphaSpin  - alphaSpinImp )*crossSpin;

        particleForce[i] += realx3
        (
            liftResidual.x(), liftResidual.y(), liftResidual.z()
        );

        liftCrossVec[i] += realx3
        (
            liftCrossVecLocal.x(),
            liftCrossVecLocal.y(),
            liftCrossVecLocal.z()
        );

        Foam::vector lF = piD2Rho * uRel2 *
            (Cl_shear * crossShear / magShear + Cl_spin * crossSpin / magSpin);

        if(this->wallCorrection())
        {
            Foam::vector nw = tNw()[cellI];

            const Foam::scalar nMag = Foam::mag(nw);

            if(nMag > Foam::SMALL)
            {
                nw /= nMag;

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
                    h/Foam::max(0.5*dp, Foam::SMALL),
                    this->minWallRadii()
                );

                const Foam::vector wallLF = wallLiftForce
                (
                    nw, uf - up, dp, rho[cellI], nu[cellI], hOverA
                );

                particleForce[i] += realx3
                (
                    wallLF.x(), wallLF.y(), wallLF.z()
                );

                lF += wallLF;

                const Foam::scalar dragScale = Foam::max
                (
                    3.0*Foam::constant::mathematical::pi
                   *rho[cellI]*nu[cellI]*dp*uRel,
                    Foam::SMALL
                );

                wallLiftObs = Foam::max
                (
                    wallLiftObs, Foam::mag(wallLF)/dragScale
                );
            }
        }

        lfPar[i] = realx3(lF.x(), lF.y(), lF.z());
    }

    this->reportWallLift(wallLiftObs);

    semiImplicitCoupling(this->Porosity().uCS()).globalParticleSum(lfPar);

    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<nPar; ++i)
    {
        const Foam::vector lF(lfPar[i].x(), lfPar[i].y(), lfPar[i].z());

        if(Foam::magSqr(lF) == 0) continue;
        if(!this->liftReaction()) continue;
        if(!cellDistribution.hasSupport(i) && parCellInd[i] < 0) continue;

        cellDistribution.distributeValue_OMP(i, parCellInd[i], liftForce, lF);
    }

    const auto& Vcells = this->mesh().V();

    forAll(Vcells, celli)
    {
        liftForce[celli] /= Vcells[celli];
    }

    liftForce.correctBoundaryConditions();

}