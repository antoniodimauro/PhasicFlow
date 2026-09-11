/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
  Copyright (C): Antonio Di Mauro
  email: antoniodimauro03@gmail.com
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
    feedback) settle; the residual test is Aitken-relaxed (Kuettler & Wall,
    Comput. Mech. 43, 61, 2008). A staggered explicit exchange is unstable
    at these density ratios (Causin, Gerbeau & Nobile, Comput. Methods Appl.
    Mech. Engrg. 194, 4506, 2005).

\*---------------------------------------------------------------------------*/

// Kuettler & Wall (2008) Comput. Mech. 43, 61

// OpenFOAM
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

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
    volVectorField r_prev
    (
        IOobject
        (
            "U_couplingResidualPrev",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        U
    );

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
        coupling.saveDEMState();
        U_prev = U;

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
        scalar alpha    = couplingRelaxInit;
        r_prev.primitiveFieldRef() = vector::zero;

        for (iter = 0; iter < couplingMaxIters; ++iter)
        {
            if (iter > 0)
            {
                U = U_save;  U.correctBoundaryConditions();
                p = p_save;  p.correctBoundaryConditions();

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

            const scalar deltaU =
                gMax(mag(r_curr.primitiveField())())/Uref;

            if (iter >= 1)
            {
                volVectorField dr("dr", r_curr - r_prev);
                const scalar dr_dot_dr =
                    gSum((dr.primitiveField() & dr.primitiveField())());
                if (dr_dot_dr > SMALL)
                {
                    const scalar r_prev_dot_dr =
                        gSum((r_prev.primitiveField() & dr.primitiveField())());
                    alpha = -alpha * r_prev_dot_dr / dr_dot_dr;
                    alpha = max(couplingRelaxFloor, min(scalar(1.0), alpha));
                }
            }

            U.primitiveFieldRef() =
                U_prev.primitiveField() + alpha * r_curr.primitiveField();
            U.correctBoundaryConditions();

            Info<< "  Coupling iter " << iter
                << " : ||dU_f||_inf/|U|max = " << deltaU
                << ", alpha = " << alpha << endl;

            if (iter >= 1 && deltaU < couplingTolerance)
            {
                converged = true;
                ++iter;
                break;
            }
            U_prev = U;
            r_prev = r_curr;
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
            coupling.restoreDEMState(tStart);

            const scalar newDt =
                max(couplingDtCutbackFactor * dt, scalar(1e-12));
            runTime.setTime(tStart, runTime.timeIndex() - 1);
            runTime.setDeltaT(newDt);
            couplingDtCap = newDt;

            Info<< "  Step rolled back to t = " << tStart
                << " with new dt = " << newDt << endl;

            coupling.cfdTimers().end();
            continue;
        }
        else
        {
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
