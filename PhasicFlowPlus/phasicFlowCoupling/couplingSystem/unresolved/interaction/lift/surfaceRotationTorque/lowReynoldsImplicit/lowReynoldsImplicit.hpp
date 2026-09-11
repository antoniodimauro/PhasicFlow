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

// Happel & Brenner (1973) Low Reynolds Number Hydrodynamics

#ifndef __lowReynoldsImplicit_hpp__
#define __lowReynoldsImplicit_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "porosity.hpp"

#include "surfaceRotationTorqueImplicit.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class lowReynoldsImplicit
:
    public surfaceRotationTorqueImplicit
{

public:

    TypeInfo("lowReynoldsImplicit");

    lowReynoldsImplicit(
        const unresolvedCouplingSystem& uCS, 
	    const porosity&                 prsty
    );

    virtual ~lowReynoldsImplicit() = default;

    add_vCtor
    (
        surfaceRotationTorque,
        lowReynoldsImplicit,
        couplingSystem
    );

    void calculateSurfaceTorqueImplicit(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque,
        Plus::realProcCMField&          rotDragCoeff) override;

};

}

#endif // __lowReynoldsImplicit_hpp__