/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Original file: phasicFlow/src/Interaction/sphereInteraction/sphereInteractionsLinearModels.cpp,
  by Hamidreza Norouzi (phasicFlow as distributed, commit 454ba4d2).
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

// Costa et al. (2015) Phys. Rev. E 92, 053012

#include "sphereInteraction.hpp"
#include "geometryMotions.hpp"
#include "contactForceModels.hpp"
#include "unsortedContactList.hpp"
#include "sortedContactList.hpp"
#include "createBoundarySphereInteraction.hpp"

#define createInteraction(ForceModel,GeomModel) 	\
													\
	template class pFlow::sphereInteraction< 		\
		ForceModel,									\
		GeomModel,									\
		pFlow::unsortedContactList>;				\
													\
	template class pFlow::sphereInteraction< 		\
		ForceModel,									\
		GeomModel,									\
		pFlow::sortedContactList>;					\
	createBoundarySphereInteraction(ForceModel, GeomModel)

createInteraction(pFlow::cfModels::limitedCostaNormalRolling, pFlow::stationaryGeometry);
createInteraction(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::stationaryGeometry);

createInteraction(pFlow::cfModels::limitedCostaNormalRolling, pFlow::rotationAxisMotionGeometry);
createInteraction(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::rotationAxisMotionGeometry);

createInteraction(pFlow::cfModels::limitedCostaNormalRolling, pFlow::vibratingMotionGeometry);
createInteraction(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::vibratingMotionGeometry);

createInteraction(pFlow::cfModels::limitedCostaNormalRolling, pFlow::conveyorBeltMotionGeometry);
createInteraction(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::conveyorBeltMotionGeometry);

createInteraction(pFlow::cfModels::limitedCostaNormalRolling, pFlow::multiRotationAxisMotionGeometry);
createInteraction(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::multiRotationAxisMotionGeometry);
