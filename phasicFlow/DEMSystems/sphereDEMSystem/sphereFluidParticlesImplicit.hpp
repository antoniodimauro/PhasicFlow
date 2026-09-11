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

#ifndef __sphereFluidParticlesImplicit_hpp__
#define __sphereFluidParticlesImplicit_hpp__

#include <vector>

#include "sphereFluidParticles.hpp"

namespace pFlow
{

class sphereFluidParticlesImplicit
:
	public sphereFluidParticles
{
private:

	realPointField_D        dragCoeff_;

	hostViewType1D<real>    dragCoeffHost_;

	realx3PointField_D      dragFluidVel_;

	hostViewType1D<realx3>  dragFluidVelHost_;

	realx3PointField_D      liftCrossVec_;

	hostViewType1D<realx3>  liftCrossVecHost_;

	realPointField_D        rotDragCoeff_;

	hostViewType1D<real>    rotDragCoeffHost_;

	realPointField_D        addedMass_;

	hostViewType1D<real>    addedMassHost_;

	realPointField_D        addedInertia_;

	hostViewType1D<real>    addedInertiaHost_;

	realPointField_D        displacedMass_;

	hostViewType1D<real>    displacedMassHost_;

	realx3PointField_D      fluidSpin_;

	hostViewType1D<realx3>  fluidSpinHost_;

	void checkHostMemoryImplicit();

public:

	TypeInfo("sphereFluidParticlesImplicit");

	sphereFluidParticlesImplicit(systemControl &control, const sphereShape& shpShape);

	bool beforeIteration() override;

	bool iterate() override;

	auto& dragCoeff()          { return dragCoeff_; }
	auto& dragCoeffHost()      { return dragCoeffHost_; }
	auto& dragFluidVel()       { return dragFluidVel_; }
	auto& dragFluidVelHost()   { return dragFluidVelHost_; }
	auto& liftCrossVec()       { return liftCrossVec_; }
	auto& liftCrossVecHost()   { return liftCrossVecHost_; }
	auto& rotDragCoeff()       { return rotDragCoeff_; }
	auto& rotDragCoeffHost()   { return rotDragCoeffHost_; }
	auto& addedMass()          { return addedMass_; }
	auto& addedMassHost()      { return addedMassHost_; }
	auto& addedInertia()       { return addedInertia_; }
	auto& addedInertiaHost()   { return addedInertiaHost_; }
	auto& displacedMass()      { return displacedMass_; }
	auto& displacedMassHost()  { return displacedMassHost_; }
	auto& fluidSpin()          { return fluidSpin_; }
	auto& fluidSpinHost()      { return fluidSpinHost_; }

	deviceViewType1D<real> addedMassDeviceView()const
	{
		return addedMass_.deviceViewAll();
	}

	void dragCoeffHostUpdatedSync();
	void dragFluidVelHostUpdatedSync();
	void liftCrossVecHostUpdatedSync();
	void rotDragCoeffHostUpdatedSync();
	void addedMassHostUpdatedSync();
	void addedInertiaHostUpdatedSync();
	void displacedMassHostUpdatedSync();
	void fluidSpinHostUpdatedSync();

	void setNewParticleVelocityFromFluid(int64 lastMaxId);

	void setNewParticleSpinFromFluid(int64 lastMaxId);

	int64 maxActiveId();

	void setMaxId(uint32 id);

	void resetIntegrationHistory();

	void saveIntegrationHistory();

	void restoreIntegrationHistory();

	void restoreKinematics(
		const std::vector<realx3>& pos,
		const std::vector<realx3>& vel,
		const std::vector<realx3>& rVel);
};

}

#endif
