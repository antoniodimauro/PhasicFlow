/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Original file: PhasicFlowPlus/solvers/unresolvedSpherePFPlus/UEqn.H,
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

// Solve the Momentum equation
// The viscous term is alpha*div(tau) (model A: Anderson & Jackson 1967, eq. 38; El Geitani, Golshan & Blais 2022,
// eq. 2). The fluid volume fraction multiplies the stress divergence as it multiplies grad(p), so the viscous force
// of the undisturbed flow on the particles, V_p div(tau), is taken from the fluid without a source term. This is the
// laminar operator of divDevReff with alpha outside the divergence: each matrix row is scaled by its cell alpha. The
// explicit part uses the divSchemes entry of the former term, the same Gauss linear scheme.
const volScalarField nuEff(turbulence->nuEff());

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(coupling.alpha(), U)
  + fvm::div(alphaRhoPhi, U)
  - coupling.alpha()()
   *(
        fvm::laplacian(nuEff, U)
      + fvc::div
        (
            nuEff*dev2(T(fvc::grad(U))),
            "div(((alpha*nuEff)*dev2(T(grad(U)))))"
        )
    )
  + fvm::Sp(coupling.Sp()/rho, U)
 ==
  - coupling.Su()/rho  
  + fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn ==  - coupling.alpha() * fvc::grad(p));

    fvOptions.correct(U);
}
