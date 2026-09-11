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

#ifndef __surfTorqueShi2019Implicit_hpp__
#define __surfTorqueShi2019Implicit_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "porosity.hpp"

#include "surfaceRotationTorqueImplicit.hpp"
#include "fluidFieldSampler.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class surfTorqueShi2019Implicit
:
    public surfaceRotationTorqueImplicit
{

    Foam::scalar    residualRe_;

    // rotationLift Bluemink2010: the rotation part of the vorticity drives the
    // spin to the torque-free rate of Bluemink et al. (2010) eq. (2.13)
    bool            rotationSpin_ = false;

    uniquePtr<fluidFieldSampler> strainSampler_;

public:

    TypeInfo("Shi2019Implicit");

    surfTorqueShi2019Implicit(
        const unresolvedCouplingSystem& uCS, 
	    const porosity&                 prsty
    );

    virtual ~surfTorqueShi2019Implicit() = default;

    add_vCtor
    (
        surfaceRotationTorque,
        surfTorqueShi2019Implicit,
        couplingSystem
    );

    // Shi & Rzehak (2019) Chem. Eng. Sci. 208, 115145
    void calculateSurfaceTorqueImplicit(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque,
        Plus::realProcCMField&          rotDragCoeff) override;

};

}

#endif // __surfTorqueShi2019Implicit_hpp__