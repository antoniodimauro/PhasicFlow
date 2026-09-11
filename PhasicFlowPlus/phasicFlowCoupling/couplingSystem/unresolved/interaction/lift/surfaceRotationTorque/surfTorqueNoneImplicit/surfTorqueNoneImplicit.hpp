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

#ifndef __surfTorqueNoneImplicit_hpp__
#define __surfTorqueNoneImplicit_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "porosity.hpp"

#include "surfaceRotationTorqueImplicit.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class surfTorqueNoneImplicit
:
    public surfaceRotationTorqueImplicit
{


public:

    TypeInfo("noneImplicit");

    surfTorqueNoneImplicit(
        const unresolvedCouplingSystem& uCS, 
	    const porosity&                 prsty
    );

    virtual ~surfTorqueNoneImplicit() = default;

    add_vCtor
    (
        surfaceRotationTorque,
        surfTorqueNoneImplicit,
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

#endif // __surfTorqueNoneImplicit_hpp__