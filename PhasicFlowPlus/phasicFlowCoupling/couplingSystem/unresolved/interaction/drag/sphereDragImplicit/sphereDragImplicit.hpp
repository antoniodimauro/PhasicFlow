/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Original file: PhasicFlowPlus/phasicFlowCoupling/couplingSystem/unresolved/interaction/drag/sphereDrag/sphereDrag.hpp,
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

// Schiller & Naumann (1933) Z. Ver. Dtsch. Ing. 77, 318

#ifndef __sphereDragImplicit_hpp__ 
#define __sphereDragImplicit_hpp__

#include <algorithm>

#include "dragImplicit.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"
#include "distributionBase.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "momentumSemiImplicitSphereUnresolvedCouplingSystem.hpp"
#include "particleGaussian.hpp"

#include "wallDist.H"
#include "fvcCurl.H"
#include "fvcGrad.H"

namespace pFlow::coupling
{

template<typename DragClosureType>
class sphereDragImplicit
:
	public dragImplicit
{
public:

	using SphereDragImplicitType = sphereDragImplicit<DragClosureType>;

private:

    DragClosureType 			dragClosure_;

public:

    // type info
    TypeInfoTemplate11("sphereDragImplicit",DragClosureType);

    sphereDragImplicit(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);

    virtual ~sphereDragImplicit() override = default ;

    add_vCtor
    (
        drag,
        SphereDragImplicitType,
        couplingSystem	
    );

    void calculateDragForceImplicit(
        const fluidAveraging& 	        fluidVelocity,
        const solidAveraging& 	        parVelocity,
        const Plus::realx3ProcCMField&  parRotVelocity,
        const Plus::realProcCMField& 	diameter,
        const distributionBase&         cellDistribution,    
        Plus::realx3ProcCMField& 		particleForce,
        Plus::realProcCMField&          dragCoeff,
        Plus::realx3ProcCMField&        dragFluidVel)override;

}; 

} // pFlow::coupling


#include "sphereDragImplicit.C"

#endif // __sphereDragImplicit_hpp__
