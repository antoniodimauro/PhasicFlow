/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Original file: PhasicFlowPlus/solvers/unresolvedSpherePFPlus/pEqn.H,
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

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA * coupling.alpha())
    // Rhie-Chow ddtCorr as in DPMFoam
  + rAUf * fvc::interpolate(coupling.alpha())*fvc::ddtCorr(U, phi, Uf)
);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

// The consistent (SIMPLEC) correction is not available: with the pressure term -alpha*grad(p) and the alpha-weighted
// flux phiHbyA it needs the factors alpha^2 on the flux and alpha on HbyA, and it has not been tested.
if (pimple.consistent())
{
    FatalErrorInFunction
        << "PIMPLE consistent (SIMPLEC) is not supported by this solver; set consistent no."
        << exit(FatalError);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

volScalarField rAUAlpha("rAUAlpha",rAtU() * coupling.alpha());

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAUAlpha);

volScalarField rAUAlpha2("rAUAlpha2",rAtU() * coupling.alpha() * coupling.alpha());

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUAlpha2, p) == fvc::div(phiHbyA) + fvc::ddt(coupling.alpha())  
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        alphaRhoPhi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - coupling.alpha() * rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// phi is the conservative fluid-velocity flux after the pressure solve, alphaRhoPhi over the face volume fraction (phic
// in DPMFoam). The Rhie-Chow term fvc::ddtCorr(U, phi, Uf) above is built on phi.oldTime() - (interpolate(U.oldTime())
// & Sf); with phi = interpolate(U) & Sf that difference was only (1 - omega) times the last pass residual of the
// Aitken relaxation, below the coupling tolerance, and the transient correction was inactive.
phi = alphaRhoPhi/fvc::interpolate(coupling.alpha());

// Correct Uf if the mesh is moving
if (mesh.dynamic())
{
    Uf() = fvc::interpolate(U * coupling.alpha());
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf() += n*(alphaRhoPhi/mesh.magSf() - (n & Uf()));
}

// Make the fluxes relative to the mesh motion
if (alphaRhoPhi.mesh().moving())
{
    alphaRhoPhi -= fvc::meshPhi(U * coupling.alpha());
}
