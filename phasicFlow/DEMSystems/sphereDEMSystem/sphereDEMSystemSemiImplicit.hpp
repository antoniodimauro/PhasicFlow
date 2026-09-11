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

#ifndef __sphereDEMSystemSemiImplicit_hpp__
#define __sphereDEMSystemSemiImplicit_hpp__

#include <vector>

#include "sphereDEMSystem.hpp"
#include "DEMSystemSemiImplicit.hpp"
#include "sphereFluidParticlesImplicit.hpp"

namespace pFlow
{

class sphereDEMSystemSemiImplicit
:
	public sphereDEMSystem,
	public DEMSystemSemiImplicit
{
private:

	sphereFluidParticlesImplicit*	implicit_ = nullptr;

	int64 							velFromFluidLastId_ = -1;

	int64 							velFromFluidPending_ = -1;

	bool 							velFromFluidFirst_ = true;

	bool 							spinFromFluid_ = false;

	uniquePtr<interaction> createLubInteraction();

public:

	TypeInfo("sphereSemiImplicitDEMSystem");

	sphereDEMSystemSemiImplicit(
		word  demSystemName,
		const std::vector<box>& domains,
		int argc,
		char* argv[],
		bool requireRVel=false);

	~sphereDEMSystemSemiImplicit() override = default;

	add_vCtor(
		DEMSystem,
		sphereDEMSystemSemiImplicit,
		word);

	span<uint32> particleId() override;

	span<real>   parDragCoeff() override;
	span<realx3> parDragFluidVel() override;
	span<realx3> parLiftCrossVec() override;
	span<real>   parRotDragCoeff() override;
	span<real>   parAddedMass() override;
	span<real>   parAddedInertia() override;
	span<real>   parDisplacedMass() override;

	bool sendDragCoeffToDEM() override;
	bool sendDragFluidVelToDEM() override;
	bool sendLiftCrossVecToDEM() override;
	bool sendRotDragCoeffToDEM() override;
	bool sendAddedMassToDEM() override;
	bool sendAddedInertiaToDEM() override;
	bool sendDisplacedMassToDEM() override;

	void resetIntegrationHistory() override;
	void saveIntegrationHistory() override;
	void restoreIntegrationHistory() override;

	void restoreParticleKinematics(
		const std::vector<realx3>& pos,
		const std::vector<realx3>& vel,
		const std::vector<realx3>& rVel) override;

	void setControlTime(real t) override
	{
		Control().time().setTime(t);
	}

	void setNewParticleVelocityFromFluid() override;
	span<realx3> parFluidSpin() override;
	bool sendFluidSpinToDEM() override;
	void setSpinFromFluid(bool on) override
	{
		spinFromFluid_ = on;
	}
	void commitNewParticleVelocityFromFluid() override;

	std::vector<uint32> insertionNumInserted() override;
	void setInsertionNumInserted(const std::vector<uint32>& vals) override;

	uint32 particleMaxId() override;
	void setParticleMaxId(uint32 id) override;
};

}

#endif
