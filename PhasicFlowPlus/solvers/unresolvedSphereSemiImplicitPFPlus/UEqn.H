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

// Solve the Momentum equation
volScalarField alphaNuEff(coupling.alpha() * turbulence->nuEff());

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(coupling.alpha(), U)
  + fvm::div(alphaRhoPhi, U)
  + turbulence->divDevReff(U)
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
