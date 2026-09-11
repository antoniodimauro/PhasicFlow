/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Original file: PhasicFlowPlus/solvers/unresolvedSpherePFPlus/createFields.H,
  by Hamidreza Norouzi and Bahram Haddadi (PhasicFlowPlus as distributed, commit bbf2110).
  Copyright (C): Antonio Di Mauro
  email: antoniodimauro03@gmail.com
  Derived from the original file in 2026 by Antonio Di Mauro
  (T-junction CFD-DEM study).
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

#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating sphere coupling system for unresolved solver \n" << endl;
auto couplingPtr = pFlow::coupling::unresolvedCouplingSystem::create("sphere", "momentumSemiImplicit", mesh, argc, argv);
auto& coupling =
    dynamic_cast<pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem&>(couplingPtr());

Info<< "Reading/calculating face volumetica flux field alphaRhoPhi\n" << endl;
surfaceScalarField alphaRhoPhi
(
    IOobject
    (
        "alphaRhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U * coupling.alpha()) & mesh.Sf()
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    // The conservative fluid-velocity flux, as set after each pressure solve (pEqn.H), from alphaRhoPhi, which is read
    // at a restart: the Rhie-Chow ddtCorr of the first step after a restart is then that of a continuous run
    alphaRhoPhi/fvc::interpolate(coupling.alpha())
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

const scalar couplingTolerance =
    pimple.dict().getOrDefault<scalar>("couplingTolerance", 1e-3);
label couplingMaxIters =
    pimple.dict().getOrDefault<label>("couplingMaxIters", 8);

if (couplingMaxIters < 2)
{
    WarningInFunction
        << "couplingMaxIters = " << couplingMaxIters
        << " cannot converge: the coupling test compares two successive"
        << " iterations, so fewer than 2 always rolls the step back and"
        << " halves dt. Raising it to 2." << endl;
    couplingMaxIters = 2;
}
const scalar couplingDtCutbackFactor =
    pimple.dict().getOrDefault<scalar>("couplingDtCutbackFactor", 0.5);
const scalar couplingDtRecoverFactor =
    pimple.dict().getOrDefault<scalar>("couplingDtRecoverFactor", 1.2);
Info<< "Coupling loop: maxIters = " << couplingMaxIters
    << ", tolerance(||dU_f||_inf/|U|max) = " << couplingTolerance
    << ", dtCutback = " << couplingDtCutbackFactor
    << ", dtRecover = " << couplingDtRecoverFactor << nl << endl;

singlePhaseTransportModel laminarTransport(U, alphaRhoPhi);


autoPtr<incompressible::alphaTurbulenceModel> turbulence
(
    incompressible::alphaTurbulenceModel::New
    (
        coupling.alpha(),
        U,
        alphaRhoPhi,
        phi,
        laminarTransport
    )
);

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("rho", dimDensity, laminarTransport)
);

#include "createFvOptions.H"

autoPtr<surfaceVectorField> Uf;

if (mesh.dynamic())
{
    Info<< "Constructing face velocity Uf\n" << endl;
    
    Uf = Foam::autoPtr<Foam::surfaceVectorField>::New
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U * coupling.alpha())
    );
}
