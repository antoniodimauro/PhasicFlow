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

// Calculate absolute flux
alphaRhoPhi = mesh.Sf() & Uf();

correctUphiBCs(U, alphaRhoPhi);

CorrectPhi
(
    U,
    alphaRhoPhi,
    p,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"

// Make the flux relative to the mesh motion
fvc::makeRelative(alphaRhoPhi, U);
