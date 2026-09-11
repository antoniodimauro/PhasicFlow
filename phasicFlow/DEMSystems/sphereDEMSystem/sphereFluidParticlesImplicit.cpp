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

#include <vector>

#include "sphereFluidParticlesImplicit.hpp"
#include "sphereFluidParticlesImplicitKernels.hpp"
#include "AdamsMoulton4PEC.hpp"
#include "regularParticleIdHandler.hpp"

void pFlow::sphereFluidParticlesImplicit::checkHostMemoryImplicit()
{
	if(dragCoeff_.size()!=dragCoeffHost_.size())
	{
		resizeNoInit(dragCoeffHost_, dragCoeff_.size());
		resizeNoInit(dragFluidVelHost_, dragFluidVel_.size());
		resizeNoInit(liftCrossVecHost_, liftCrossVec_.size());
		resizeNoInit(rotDragCoeffHost_, rotDragCoeff_.size());
		resizeNoInit(addedMassHost_, addedMass_.size());
		resizeNoInit(addedInertiaHost_, addedInertia_.size());
		resizeNoInit(displacedMassHost_, displacedMass_.size());
		resizeNoInit(fluidSpinHost_, fluidSpin_.size());
	}
}

pFlow::sphereFluidParticlesImplicit::sphereFluidParticlesImplicit(
    systemControl &control,
    const sphereShape& shpShape)
    : sphereFluidParticles(control, shpShape),
      dragCoeff_(
          objectFile(
              "dragCoeff",
              "",
              objectFile::READ_IF_PRESENT,
              objectFile::WRITE_NEVER),
          dynPointStruct(),
          static_cast<real>(0)),
      dragFluidVel_(
          objectFile(
              "dragFluidVel",
              "",
              objectFile::READ_IF_PRESENT,
              objectFile::WRITE_NEVER),
          dynPointStruct(),
          zero3),
      liftCrossVec_(
          objectFile(
              "liftCrossVec",
              "",
              objectFile::READ_IF_PRESENT,
              objectFile::WRITE_NEVER),
          dynPointStruct(),
          zero3),
      rotDragCoeff_(
          objectFile(
              "rotDragCoeff",
              "",
              objectFile::READ_IF_PRESENT,
              objectFile::WRITE_NEVER),
          dynPointStruct(),
          static_cast<real>(0)),
      addedMass_(
          objectFile(
              "addedMass",
              "",
              objectFile::READ_IF_PRESENT,
              objectFile::WRITE_NEVER),
          dynPointStruct(),
          static_cast<real>(0)),
      addedInertia_(
          objectFile(
              "addedInertia",
              "",
              objectFile::READ_IF_PRESENT,
              objectFile::WRITE_NEVER),
          dynPointStruct(),
          static_cast<real>(0)),
      displacedMass_(
          objectFile(
              "displacedMass",
              "",
              objectFile::READ_IF_PRESENT,
              objectFile::WRITE_NEVER),
          dynPointStruct(),
          static_cast<real>(0)),
      fluidSpin_(
          objectFile(
              "fluidSpin",
              "",
              objectFile::READ_IF_PRESENT,
              objectFile::WRITE_NEVER),
          dynPointStruct(),
          zero3)
{
	checkHostMemoryImplicit();
}

bool pFlow::sphereFluidParticlesImplicit::beforeIteration()
{
	sphereFluidParticles::beforeIteration();
	checkHostMemoryImplicit();

	return true;
}

bool pFlow::sphereFluidParticlesImplicit::iterate() 
{
	const auto ti = this->TimeInfo();

	auto* amPos  = dynamic_cast<AdamsMoulton4PEC*>(&dynPointStruct().integrationPos());
	auto* amVel  = dynamic_cast<AdamsMoulton4PEC*>(&dynPointStruct().integrationVel());
	auto* amRVel = dynamic_cast<AdamsMoulton4PEC*>(&rVelIntegration());

	if(amPos && amVel && amRVel)
	{
		accelerationTimer().start();
		intCorrectTimer().start();

		auto hP = amPos->getDeviceHistory();
		auto hV = amVel->getDeviceHistory();
		auto hW = amRVel->getDeviceHistory();

		pFlow::sphereFluidParticlesImplicitKernels::correctFusedAM4(
			ti.dt(),
			dynPointStruct().dampingFactor(ti),
			control().g(),
			mass().deviceViewAll(),
			addedMass_.deviceViewAll(),
			displacedMass_.deviceViewAll(),
			I().deviceViewAll(),
			addedInertia_.deviceViewAll(),
			contactForce().deviceViewAll(),
			fluidForce_.deviceViewAll(),
			contactTorque().deviceViewAll(),
			fluidTorque_.deviceViewAll(),
			acceleration().deviceViewAll(),
			rAcceleration().deviceViewAll(),
			pStruct().pointPosition().deviceView(),
			dynPointStruct().velocity().deviceViewAll(),
			rVelocity().deviceViewAll(),
			hP.dy0, hP.dy1, hP.dy2, hP.dy3,
			hV.dy0, hV.dy1, hV.dy2, hV.dy3,
			hW.dy0, hW.dy1, hW.dy2, hW.dy3,
			dragCoeff_.deviceViewAll(),
			dragFluidVel_.deviceViewAll(),
			liftCrossVec_.deviceViewAll(),
			rotDragCoeff_.deviceViewAll(),
			pStruct().activePointsMaskDevice()
			);

		amPos->correctPStructBoundaryOnly(ti.dt(), dynPointStruct(), dynPointStruct().velocity());
		amVel->correctBoundaryOnly(ti.dt(), dynPointStruct().velocity(), acceleration());
		amRVel->correctBoundaryOnly(ti.dt(), rVelocity(), rAcceleration());

		intCorrectTimer().end();
		accelerationTimer().end();
		return true;
	}

	accelerationTimer().start();
		pFlow::sphereFluidParticlesImplicitKernels::acceleration(
			control().g(),
			mass().deviceViewAll(),
			addedMass_.deviceViewAll(),
			displacedMass_.deviceViewAll(),
			contactForce().deviceViewAll(),
			fluidForce_.deviceViewAll(),
			I().deviceViewAll(),
			addedInertia_.deviceViewAll(),
			contactTorque().deviceViewAll(),
			fluidTorque_.deviceViewAll(),
			pStruct().activePointsMaskDevice(),
			acceleration().deviceViewAll(),
			rAcceleration().deviceViewAll()
			);
	accelerationTimer().end();
	
	intCorrectTimer().start();
	
		dynPointStruct().correct(ti.dt());
	
		rVelIntegration().correct(ti.dt(), rVelocity(), rAcceleration());

		pFlow::sphereFluidParticlesImplicitKernels::applyImplicitDragLift(
			ti.dt(),
			mass().deviceViewAll(),
			addedMass_.deviceViewAll(),
			dragCoeff_.deviceViewAll(),
			dragFluidVel_.deviceViewAll(),
			liftCrossVec_.deviceViewAll(),
			pStruct().activePointsMaskDevice(),
			dynPointStruct().velocity().deviceViewAll()
			);

		pFlow::sphereFluidParticlesImplicitKernels::applyImplicitRotDrag(
			ti.dt(),
			I().deviceViewAll(),
			addedInertia_.deviceViewAll(),
			rotDragCoeff_.deviceViewAll(),
			pStruct().activePointsMaskDevice(),
			rVelocity().deviceViewAll()
			);
	
	intCorrectTimer().end();
	
	return true;
}

void pFlow::sphereFluidParticlesImplicit::dragCoeffHostUpdatedSync()
{
	copy(dragCoeff_.deviceView(), dragCoeffHost_);
	return;
}

void pFlow::sphereFluidParticlesImplicit::dragFluidVelHostUpdatedSync()
{
	copy(dragFluidVel_.deviceView(), dragFluidVelHost_);
	return;
}

void pFlow::sphereFluidParticlesImplicit::liftCrossVecHostUpdatedSync()
{
	copy(liftCrossVec_.deviceView(), liftCrossVecHost_);
	return;
}

void pFlow::sphereFluidParticlesImplicit::rotDragCoeffHostUpdatedSync()
{
	copy(rotDragCoeff_.deviceView(), rotDragCoeffHost_);
	return;
}

void pFlow::sphereFluidParticlesImplicit::addedMassHostUpdatedSync()
{
	copy(addedMass_.deviceView(), addedMassHost_);
	return;
}

void pFlow::sphereFluidParticlesImplicit::addedInertiaHostUpdatedSync()
{
	copy(addedInertia_.deviceView(), addedInertiaHost_);
	return;
}

void pFlow::sphereFluidParticlesImplicit::displacedMassHostUpdatedSync()
{
	copy(displacedMass_.deviceView(), displacedMassHost_);
	return;
}

void pFlow::sphereFluidParticlesImplicit::fluidSpinHostUpdatedSync()
{
	copy(fluidSpin_.deviceView(), fluidSpinHost_);
	return;
}

void pFlow::sphereFluidParticlesImplicit::setNewParticleSpinFromFluid(int64 lastMaxId)
{
	pFlow::sphereFluidParticlesImplicitKernels::setVelocityFromFluidForNew(
		lastMaxId,
		particleId().deviceViewAll(),
		fluidSpin_.deviceViewAll(),
		pStruct().activePointsMaskDevice(),
		rVelocity().deviceViewAll()
		);
}

void pFlow::sphereFluidParticlesImplicit::setNewParticleVelocityFromFluid(int64 lastMaxId)
{
	pFlow::sphereFluidParticlesImplicitKernels::setVelocityFromFluidForNew(
		lastMaxId,
		particleId().deviceViewAll(),
		dragFluidVel_.deviceViewAll(),
		pStruct().activePointsMaskDevice(),
		dynPointStruct().velocity().deviceViewAll()
		);
}

pFlow::int64 pFlow::sphereFluidParticlesImplicit::maxActiveId()
{
	if(pStruct().numActive() == 0u) return -1;
	return pFlow::sphereFluidParticlesImplicitKernels::maxActiveId(
		particleId().deviceViewAll(),
		pStruct().activePointsMaskDevice()
		);
}
void pFlow::sphereFluidParticlesImplicit::setMaxId(uint32 id)
{
	auto* handler = dynamic_cast<regularParticleIdHandler*>(&idHandler());
	if(!handler)
	{
		fatalErrorInFunction<<
		"the semi-implicit coupling needs the regular particle id handler"<<endl;
		fatalExit;
	}
	handler->setMaxId(id);
}

namespace pFlow
{
namespace
{
	AdamsMoulton4PEC& pecOrFatal(integration& intg)
	{
		auto* pec = dynamic_cast<AdamsMoulton4PEC*>(&intg);
		if(!pec)
		{
			fatalErrorInFunction<<
			"the semi-implicit coupling needs integrationMethod AdamsMoulton4PEC, found "<<
			intg.method()<<endl;
			fatalExit;
		}
		return *pec;
	}
}
}

void pFlow::sphereFluidParticlesImplicit::resetIntegrationHistory()
{
	pecOrFatal(dynPointStruct().integrationPos()).resetHistory();
	pecOrFatal(dynPointStruct().integrationVel()).resetHistory();
	pecOrFatal(rVelIntegration()).resetHistory();
}

void pFlow::sphereFluidParticlesImplicit::saveIntegrationHistory()
{
	pecOrFatal(dynPointStruct().integrationPos()).saveHistory();
	pecOrFatal(dynPointStruct().integrationVel()).saveHistory();
	pecOrFatal(rVelIntegration()).saveHistory();
}

void pFlow::sphereFluidParticlesImplicit::restoreIntegrationHistory()
{
	pecOrFatal(dynPointStruct().integrationPos()).restoreHistory();
	pecOrFatal(dynPointStruct().integrationVel()).restoreHistory();
	pecOrFatal(rVelIntegration()).restoreHistory();
}

void pFlow::sphereFluidParticlesImplicit::restoreKinematics(
	const std::vector<realx3>& pos,
	const std::vector<realx3>& vel,
	const std::vector<realx3>& rVel)
{
	if(!pos.empty()) dynPointStruct().pointPosition().assign(pos);
	if(!vel.empty()) dynPointStruct().velocity().field().assign(vel);
	dynPointStruct().markAllInternalActive(static_cast<uint32>(pos.size()));
	if(!rVel.empty()) rVelocity().field().assign(rVel);
}
