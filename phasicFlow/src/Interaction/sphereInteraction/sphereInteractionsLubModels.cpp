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

#include "sphereInteractionLub.hpp"
#include "geometryMotions.hpp"
#include "contactForceModels.hpp"
#include "unsortedContactList.hpp"
#include "sortedContactList.hpp"

#define createInteractionLub(ForceModel,GeomModel) 	\
													\
	template class pFlow::sphereInteractionLub< 	\
		ForceModel,									\
		GeomModel,									\
		pFlow::unsortedContactList>;				\
													\
	template class pFlow::sphereInteractionLub< 	\
		ForceModel,									\
		GeomModel,									\
		pFlow::sortedContactList>;

createInteractionLub(pFlow::cfModels::limitedCostaNormalRolling, pFlow::stationaryGeometry);
createInteractionLub(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::stationaryGeometry);
createInteractionLub(pFlow::cfModels::limitedCostaNormalRolling, pFlow::rotationAxisMotionGeometry);
createInteractionLub(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::rotationAxisMotionGeometry);
createInteractionLub(pFlow::cfModels::limitedCostaNormalRolling, pFlow::vibratingMotionGeometry);
createInteractionLub(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::vibratingMotionGeometry);
createInteractionLub(pFlow::cfModels::limitedCostaNormalRolling, pFlow::conveyorBeltMotionGeometry);
createInteractionLub(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::conveyorBeltMotionGeometry);
createInteractionLub(pFlow::cfModels::limitedCostaNormalRolling, pFlow::multiRotationAxisMotionGeometry);
createInteractionLub(pFlow::cfModels::nonLimitedCostaNormalRolling,pFlow::multiRotationAxisMotionGeometry);
