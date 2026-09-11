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

#ifndef __DEMSystemSemiImplicit_hpp__
#define __DEMSystemSemiImplicit_hpp__

#include <vector>

#include "types.hpp"
#include "span.hpp"

namespace pFlow
{

class DEMSystemSemiImplicit
{
public:

	virtual ~DEMSystemSemiImplicit() = default;

	virtual span<real>   parDragCoeff() = 0;
	virtual span<realx3> parDragFluidVel() = 0;
	virtual span<realx3> parLiftCrossVec() = 0;
	virtual span<real>   parRotDragCoeff() = 0;
	virtual span<real>   parAddedMass() = 0;
	virtual span<real>   parAddedInertia() = 0;
	virtual span<real>   parDisplacedMass() = 0;

	virtual bool sendDragCoeffToDEM() = 0;
	virtual bool sendDragFluidVelToDEM() = 0;
	virtual bool sendLiftCrossVecToDEM() = 0;
	virtual bool sendRotDragCoeffToDEM() = 0;
	virtual bool sendAddedMassToDEM() = 0;
	virtual bool sendAddedInertiaToDEM() = 0;
	virtual bool sendDisplacedMassToDEM() = 0;

	virtual void resetIntegrationHistory() = 0;
	virtual void saveIntegrationHistory() = 0;
	virtual void restoreIntegrationHistory() = 0;

	virtual void restoreParticleKinematics(
		const std::vector<realx3>& pos,
		const std::vector<realx3>& vel,
		const std::vector<realx3>& rVel) = 0;

	virtual void setControlTime(real t) = 0;

	virtual void setNewParticleVelocityFromFluid() = 0;
	virtual span<realx3> parFluidSpin() = 0;
	virtual bool sendFluidSpinToDEM() = 0;
	virtual void setSpinFromFluid(bool on) = 0;
	virtual void commitNewParticleVelocityFromFluid() = 0;

	virtual std::vector<uint32> insertionNumInserted() = 0;
	virtual void setInsertionNumInserted(const std::vector<uint32>& vals) = 0;

	virtual uint32 particleMaxId() = 0;
	virtual void setParticleMaxId(uint32 id) = 0;
};

}

#endif
