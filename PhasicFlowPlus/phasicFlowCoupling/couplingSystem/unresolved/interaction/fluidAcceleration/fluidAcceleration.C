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

#include "fvc.H"
#include "fvMesh.H"
#include "volFields.H"

#include "fluidAcceleration.hpp"

Foam::tmp<Foam::volVectorField> pFlow::coupling::materialDerivative
(
    const Foam::volVectorField& U
)
{
    // unresolvedSphereSemiImplicitPFPlus exchanges before the fluid solve of the step, where fvc::ddt(U) is
    // (U^n - U^n)/dt = 0; it registers the lagged time derivative of the last accepted step as "ddtULagged"
    const Foam::fvMesh& mesh = U.mesh();
    if (mesh.foundObject<Foam::volVectorField>("ddtULagged"))
    {
        return mesh.lookupObject<Foam::volVectorField>("ddtULagged") + (U & Foam::fvc::grad(U));
    }
    return Foam::fvc::ddt(U) + (U & Foam::fvc::grad(U));
}
