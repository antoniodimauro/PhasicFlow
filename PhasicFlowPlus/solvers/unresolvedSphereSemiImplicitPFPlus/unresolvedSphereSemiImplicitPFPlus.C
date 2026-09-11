/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Original file: PhasicFlowPlus/solvers/unresolvedSpherePFPlus/unresolvedSpherePFPlus.C,
  by Hamidreza Norouzi and Bahram Haddadi (PhasicFlowPlus as distributed, commit bbf2110).
  Copyright (C): Antonio Di Mauro
  email: antoniodimauro03@gmail.com
  Derived from the original file in 2026 by Antonio Di Mauro
  (T-junction CFD-DEM study).
------------------------------------------------------------------------------
License
    This file is partially part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This file is part of PhasicFlow, a CFD-DEM stack built on phasicFlow and
    PhasicFlowPlus (www.cemf.ir), and is distributed under the same licence,
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    unresolvedSphereSemiImplicitPFPlus

Description
    Unresolved sphere-fluid CFD-DEM solver with a semi-implicit
    fluid-particle exchange: one exchange per time step from the start-of-step
    state, the drag implicit in the fluid (fvm::Sp) and in the particle
    equation, the added and Basset masses on the particle inertia. The step is
    repeated from the rewound state until the lagged terms (pair forces, lift
    feedback) settle: the step is accepted when two successive passes agree,
    with the velocity of the last pass. A staggered explicit exchange is
    unstable at these density ratios (Causin, Gerbeau & Nobile, Comput.
    Methods Appl. Mech. Engrg. 194, 4506, 2005). The exchange runs before the
    fluid solve of the step, so the fluid time derivative in Du/Dt is the
    lagged (U^n - U^(n-1))/dt^(n-1) of the last accepted step.

\*---------------------------------------------------------------------------*/

// OpenFOAM
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "uniformDimensionedFields.H"

// phasicFlow
#include "momentumSemiImplicitSphereUnresolvedCouplingSystem.hpp"
#include "alphaTurbulentTransportModel.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"

    pFlow::Plus::processor::initMPI(argc, argv);

    #include "createFields.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    volVectorField U_save
    (
        IOobject
        (
            "U_couplingSave",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        U
    );
    volScalarField p_save
    (
        IOobject
        (
            "p_couplingSave",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        p
    );

    // C9: the face fluxes are rewound with U and p, so every pass of a step, the rollback and the write repeat start
    // PIMPLE from the same state (QUESTIONS.md G9)
    surfaceScalarField phi_save
    (
        IOobject("phi_couplingSave", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        phi
    );
    surfaceScalarField alphaRhoPhi_save
    (
        IOobject("alphaRhoPhi_couplingSave", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        alphaRhoPhi
    );
    volVectorField U_prev
    (
        IOobject
        (
            "U_couplingPrev",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        U
    );

    scalar couplingDtCap = GREAT;

    volVectorField r_curr
    (
        IOobject
        (
            "U_couplingResidualCurr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        U
    );

    // Fluid time derivative for the exchange (21 Oct 2026). The exchange runs before the fluid solve of the step, where
    // fvc::ddt(U) = (U^n - U^n)/dt = 0; it uses (U^n - U^(n-1))/dt^(n-1) of the last accepted step instead, the same in
    // every pass, the rollback and the write pass of a step. U^(n-1) and dt^(n-1) change only when a step is accepted,
    // are written with the fields of every save (runTime.write) and read at a restart; without them (t = 0, or a save of
    // an earlier code) the step uses 0. fluidAcceleration.C (materialDerivative) reads "ddtULagged" when it exists.
    IOobject U_nm1Header
    (
        "U_nm1", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE
    );
    const bool haveU_nm1 =
        returnReduce(U_nm1Header.typeHeaderOk<volVectorField>(true), andOp<bool>());
    autoPtr<volVectorField> U_nm1Ptr
    (
        haveU_nm1
      ? new volVectorField(U_nm1Header, mesh)
      : new volVectorField                  // U's patch types, not a copy of U: a copy would carry U's old-time level
        (                                   // and write it as U_nm1_0 at every save
            IOobject("U_nm1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            mesh,
            dimensionedVector(U.dimensions(), Zero),
            U.boundaryField().types()
        )
    );
    volVectorField& U_nm1 = U_nm1Ptr.ref();
    if (!haveU_nm1)
    {
        U_nm1 == U;
    }
    uniformDimensionedScalarField dt_nm1
    (
        IOobject
        (
            "dt_nm1", runTime.timeName(), "uniform", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        dimensionedScalar(dimTime, -1)
    );
    bool lagValid = haveU_nm1 && dt_nm1.value() > 0;
    volVectorField ddtULagged
    (
        IOobject("ddtULagged", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector(dimVelocity/dimTime, Zero)
    );
    if (lagValid)
    {
        Info<< "Lagged ddt(U): U_nm1 and dt_nm1 read, dt_nm1 = " << dt_nm1.value() << endl;
    }
    else
    {
        Info<< "Lagged ddt(U): no U_nm1/dt_nm1 at the start time, the first step uses 0" << endl;
    }
    bool stepRetry = false;   // true while a rolled-back step is repeated: ddtULagged is kept

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;


    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        if (couplingDtCap < runTime.deltaT().value())
        {
            Info<< "Coupling cap: shrinking dt from "
                << runTime.deltaT().value() << " to " << couplingDtCap << endl;
            runTime.setDeltaT(couplingDtCap);
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;        

        coupling.cfdTimers().start();        

        const scalar t      = runTime.time().value();
        const scalar dt     = runTime.deltaT().value();
        const scalar tStart = t - dt;

        U_save = U;
        p_save = p;
        phi_save = phi;
        alphaRhoPhi_save = alphaRhoPhi;
        coupling.saveDEMState();
        U_prev = U;

        if (!stepRetry)
        {
            if (lagValid)
            {
                ddtULagged == (U_save - U_nm1)/dimensionedScalar("dt_nm1", dimTime, dt_nm1.value());
            }
            else
            {
                ddtULagged == dimensionedVector(dimVelocity/dimTime, Zero);
            }
            Info<< "  Lagged ddt(U): max |dU/dt| = " << gMax(mag(ddtULagged.primitiveField())())
                << " m/s^2" << (lagValid ? "" : " (no previous accepted step: 0)") << endl;
        }
        stepRetry = false;

        scalar Uref = gMax(mag(U_save.primitiveField())());

        forAll(U_save.boundaryField(), patchI)
        {
            const fvPatchVectorField& pf = U_save.boundaryField()[patchI];

            if (pf.size() > 0)
            {
                Uref = max(Uref, max(mag(pf)()));
            }
        }

        reduce(Uref, maxOp<scalar>());

        Uref = max(Uref, SMALL);

        bool  converged = false;
        label iter      = 0;

        for (iter = 0; iter < couplingMaxIters; ++iter)
        {
            if (iter > 0)
            {
                U = U_save;  U.correctBoundaryConditions();
                p = p_save;  p.correctBoundaryConditions();
                // The face fluxes are kept from the previous pass: they are the latest iterate of the step's nonlinear
                // solution (their old time, used by ddtCorr, is phi^n either way; the exchange does not read them), and
                // a step that PIMPLE does not converge within its outer correctors goes on converging over the passes.
                // Rewinding them too made every pass repeat the first and accepted such a step unconverged (pairs
                // Re1-single, 15 Oct 2026). The rollback restores them: they come from a failed attempt.

                coupling.restoreDEMState(tStart);
            }
                
        while (pimple.loop())
        {
                if (pimple.firstIter())
            {
                coupling.cfdTimers().pause();
                
                coupling.getDataFromDEM(t, dt);
                coupling.calculatePorosity();
                coupling.calculateMomentumCoupling();
                coupling.sendDataToDEM(t, dt);
                
                    coupling.applyFluidVelocityToNewParticles();
                    coupling.iterate(t, false, runTime.timeName());
                
                coupling.cfdTimers().start();
            }
                
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.controlledUpdate();
                
                if (mesh.changing())
                {
                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }


            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

            r_curr.primitiveFieldRef() =
                U.primitiveField() - U_prev.primitiveField();

            // The scale is the larger of the largest velocities at the start of the step and after this pass. From rest
            // (U_ref ~ 0) the start-of-step value alone made a pass difference at round-off of the new velocity look
            // like 1e-4, and converged steps ran to couplingMaxIters and were rolled back (ten Cate E2, E3, 15 Oct 2026)
            const scalar Ucur = gMax(mag(U.primitiveField())());
            const scalar deltaU =
                gMax(mag(r_curr.primitiveField())())/max(Uref, Ucur);

            // The accepted velocity is the PIMPLE solution of the last pass, consistent with its p and phi. An Aitken
            // blend U_prev + alpha r_curr was applied here until 21 Oct 2026: every pass restarts from u^n, so the blend
            // never reached the next pass and only shifted the accepted field, by less than the coupling tolerance.
            Info<< "  Coupling iter " << iter
                << " : ||dU_f||_inf/|U|max = " << deltaU << endl;

            if (iter >= 1 && deltaU < couplingTolerance)
            {
                converged = true;
                ++iter;
                break;
            }
            U_prev = U;
        }

        bool particleLoss = false;
        {
            const label currentCount = static_cast<label>(coupling.numParticlesMaster());
            const label savedCount   = static_cast<label>(coupling.numSavedParticlesMaster());
            const label minCount     = returnReduce(currentCount, minOp<label>());
            const label minSaved     = returnReduce(savedCount,   minOp<label>());
            if (currentCount < savedCount)
            {
                particleLoss = true;
                Info<< "  Particle count dropped: " << savedCount
                    << " -> " << currentCount << " on master." << endl;
            }
            (void)minCount; (void)minSaved;
        }

        if (!converged || particleLoss)
        {
            Info<< "  Coupling "
                << (converged ? "converged but particles lost" :
                    "NOT converged after " + Foam::name(iter) + " iters")
                << "; rolling back step and halving dt." << endl;

            U = U_save;  U.correctBoundaryConditions();
            p = p_save;  p.correctBoundaryConditions();
            phi = phi_save;  alphaRhoPhi = alphaRhoPhi_save;
            coupling.restoreDEMState(tStart);

            const scalar newDt =
                max(couplingDtCutbackFactor * dt, scalar(1e-12));
            runTime.setTime(tStart, runTime.timeIndex() - 1);
            runTime.setDeltaT(newDt);
            couplingDtCap = newDt;

            Info<< "  Step rolled back to t = " << tStart
                << " with new dt = " << newDt << endl;

            stepRetry = true;   // the repeat keeps this step's ddtULagged; U_nm1, dt_nm1 untouched
            coupling.cfdTimers().end();
            continue;
        }
        else
        {
            // accepted: U^n and this step become the lagged level of the next step (forced assignment: also the
            // fixed-value patches)
            U_nm1 == U_save;
            dt_nm1.value() = dt;
            lagValid = true;

            if (couplingDtCap < GREAT)
            {
                couplingDtCap = min
                (
                    GREAT,
                    couplingDtRecoverFactor * couplingDtCap
                );
                if (couplingDtCap >= dt * 1.5)
                {
                    couplingDtCap = GREAT;
                }
            }
        }

        if (runTime.writeTime())
        {
            // The DEM output is written by a repeat of the accepted exchange:
            // same start-of-step fluid fields and rewound DEM state as every
            // coupling iteration, then the converged fluid fields back.
            coupling.cfdTimers().pause();

            volVectorField U_conv
            (
                IOobject
                (
                    "U_couplingConv",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                U
            );
            volScalarField p_conv
            (
                IOobject
                (
                    "p_couplingConv",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                p
            );

            // The exchange reads U and p, not the face fluxes, which keep their converged values of t^{n+1} (setting
            // them back to phi^n here left every step after a write with the fluxes of the step before)
            U = U_save;  U.correctBoundaryConditions();
            p = p_save;  p.correctBoundaryConditions();
            coupling.restoreDEMState(tStart);

            coupling.getDataFromDEM(t, dt);
            coupling.calculatePorosity();
            coupling.calculateMomentumCoupling();
            coupling.sendDataToDEM(t, dt);
            coupling.applyFluidVelocityToNewParticles();
            coupling.iterate(t, true, runTime.timeName());

            U = U_conv;  U.correctBoundaryConditions();
            p = p_conv;  p.correctBoundaryConditions();

            coupling.cfdTimers().start();
        }

        coupling.commitHistoryForces();

        coupling.commitNewParticleVelocity();


        coupling.cfdTimers().end();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    pFlow::Plus::processor::finalizeMPI();

    return 0;
}


// ************************************************************************* //
