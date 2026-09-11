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
#include <utility>

#include "sphereDEMSystemSemiImplicit.hpp"
#include "vocabs.hpp"

pFlow::uniquePtr<pFlow::interaction>
pFlow::sphereDEMSystemSemiImplicit::createLubInteraction()
{
	fileSystem file = Control().caseSetup().path()+interactionFile__;
	dictionary dict(interactionFile__, file);

	auto interactionDict = dict.subDict("model");

	word clType  = dict.getVal<word>("contactListType");
	word cfModel = interactionDict.getVal<word>("contactForceModel");
	word rfModel = interactionDict.getVal<word>("rollingFrictionModel");

	auto interactionModel = angleBracketsNames3(
			Particles().shapeTypeName()+"InteractionLub",
			angleBracketsNames(rfModel,cfModel),
			Geometry().motionModelTypeName(),
			clType);

	REPORT(0)<<"Creating interaction "<<Green_Text(interactionModel)<<" . . ."<<END_REPORT;

	if( interaction::systemControlvCtorSelector_.search(interactionModel) )
	{
		return interaction::systemControlvCtorSelector_[interactionModel]
			(Control(), Particles(), Geometry());
	}

	printKeys
	(
		fatalError << "Ctor Selector "<< interactionModel << " dose not exist. \n"
		<<"Avaiable ones are: \n\n"
		,
		interaction::systemControlvCtorSelector_
	);
	fatalExit;
	return nullptr;
}

pFlow::sphereDEMSystemSemiImplicit::sphereDEMSystemSemiImplicit(
		word  demSystemName,
		const std::vector<box>& domains,
		int argc,
		char* argv[],
		bool requireRVel)
:
	sphereDEMSystem(demSystemName, domains, argc, argv, requireRVel)
{
	particleDistribution_.reset();
	interaction_.reset();
	insertion_.reset();
	particles_.reset();

	REPORT(0)<<"\nReading sphere particles (implicit coupling fields) . . ."<<END_REPORT;
	auto implicitPtr = makeUnique<sphereFluidParticlesImplicit>(Control(), spheres_());
	implicit_ = implicitPtr.get();
	particles_ = std::move(implicitPtr);

	insertion_ = makeUnique<sphereInsertion>(
		particles_(),
		particles_().spheres());

	REPORT(0)<<"\nCreating interaction model for sphere-sphere contact, sphere-wall contact and lubrication . . ."<<END_REPORT;
	interaction_ = createLubInteraction();

	real minD, maxD;
	particles_->boundingSphereMinMax(minD, maxD);
	particleDistribution_ = makeUnique<domainDistribute>(domains, maxD);
}

pFlow::span<pFlow::uint32> pFlow::sphereDEMSystemSemiImplicit::particleId()
{
	return span<uint32>(particleIdHost_.data(), particleIdHost_.size());
}

pFlow::span<pFlow::real> pFlow::sphereDEMSystemSemiImplicit::parDragCoeff()
{
	auto& hVec = implicit_->dragCoeffHost();
	return span<real>(hVec.data(), hVec.size());
}

pFlow::span<pFlow::realx3> pFlow::sphereDEMSystemSemiImplicit::parDragFluidVel()
{
	auto& hVec = implicit_->dragFluidVelHost();
	return span<realx3>(hVec.data(), hVec.size());
}

pFlow::span<pFlow::realx3> pFlow::sphereDEMSystemSemiImplicit::parFluidSpin()
{
	auto& hVec = implicit_->fluidSpinHost();
	return span<realx3>(hVec.data(), hVec.size());
}

bool pFlow::sphereDEMSystemSemiImplicit::sendFluidSpinToDEM()
{
	implicit_->fluidSpinHostUpdatedSync();
	return true;
}

pFlow::span<pFlow::realx3> pFlow::sphereDEMSystemSemiImplicit::parLiftCrossVec()
{
	auto& hVec = implicit_->liftCrossVecHost();
	return span<realx3>(hVec.data(), hVec.size());
}

pFlow::span<pFlow::real> pFlow::sphereDEMSystemSemiImplicit::parRotDragCoeff()
{
	auto& hVec = implicit_->rotDragCoeffHost();
	return span<real>(hVec.data(), hVec.size());
}

pFlow::span<pFlow::real> pFlow::sphereDEMSystemSemiImplicit::parAddedMass()
{
	auto& hVec = implicit_->addedMassHost();
	return span<real>(hVec.data(), hVec.size());
}

pFlow::span<pFlow::real> pFlow::sphereDEMSystemSemiImplicit::parAddedInertia()
{
	auto& hVec = implicit_->addedInertiaHost();
	return span<real>(hVec.data(), hVec.size());
}

pFlow::span<pFlow::real> pFlow::sphereDEMSystemSemiImplicit::parDisplacedMass()
{
	auto& hVec = implicit_->displacedMassHost();
	return span<real>(hVec.data(), hVec.size());
}

bool pFlow::sphereDEMSystemSemiImplicit::sendDragCoeffToDEM()
{
	implicit_->dragCoeffHostUpdatedSync();
	return true;
}

bool pFlow::sphereDEMSystemSemiImplicit::sendDragFluidVelToDEM()
{
	implicit_->dragFluidVelHostUpdatedSync();
	return true;
}

bool pFlow::sphereDEMSystemSemiImplicit::sendLiftCrossVecToDEM()
{
	implicit_->liftCrossVecHostUpdatedSync();
	return true;
}

bool pFlow::sphereDEMSystemSemiImplicit::sendRotDragCoeffToDEM()
{
	implicit_->rotDragCoeffHostUpdatedSync();
	return true;
}

bool pFlow::sphereDEMSystemSemiImplicit::sendAddedMassToDEM()
{
	implicit_->addedMassHostUpdatedSync();
	return true;
}

bool pFlow::sphereDEMSystemSemiImplicit::sendAddedInertiaToDEM()
{
	implicit_->addedInertiaHostUpdatedSync();
	return true;
}

bool pFlow::sphereDEMSystemSemiImplicit::sendDisplacedMassToDEM()
{
	implicit_->displacedMassHostUpdatedSync();
	return true;
}

void pFlow::sphereDEMSystemSemiImplicit::resetIntegrationHistory()
{
	implicit_->resetIntegrationHistory();
}

void pFlow::sphereDEMSystemSemiImplicit::saveIntegrationHistory()
{
	implicit_->saveIntegrationHistory();
}

void pFlow::sphereDEMSystemSemiImplicit::restoreIntegrationHistory()
{
	implicit_->restoreIntegrationHistory();
}

void pFlow::sphereDEMSystemSemiImplicit::restoreParticleKinematics(
	const std::vector<realx3>& pos,
	const std::vector<realx3>& vel,
	const std::vector<realx3>& rVel)
{
	implicit_->restoreKinematics(pos, vel, rVel);
}

void pFlow::sphereDEMSystemSemiImplicit::setNewParticleVelocityFromFluid()
{
	const int64 activeMaxId = implicit_->maxActiveId();

	if(velFromFluidFirst_)
	{
		velFromFluidFirst_   = false;
		velFromFluidLastId_  = activeMaxId;
		velFromFluidPending_ = activeMaxId;
		return;
	}

	if(activeMaxId > velFromFluidLastId_)
	{
		implicit_->setNewParticleVelocityFromFluid(velFromFluidLastId_);
		if(spinFromFluid_)
		{
			implicit_->setNewParticleSpinFromFluid(velFromFluidLastId_);
		}
	}

	velFromFluidPending_ = activeMaxId;
}

void pFlow::sphereDEMSystemSemiImplicit::commitNewParticleVelocityFromFluid()
{
	velFromFluidLastId_ = velFromFluidPending_;
}

std::vector<pFlow::uint32> pFlow::sphereDEMSystemSemiImplicit::insertionNumInserted()
{
	if(insertion_)
	{
		return insertion_().getAllNumInserted();
	}
	return {};
}

void pFlow::sphereDEMSystemSemiImplicit::setInsertionNumInserted(const std::vector<uint32>& vals)
{
	if(insertion_)
	{
		insertion_().setAllNumInserted(vals);
	}
}

pFlow::uint32 pFlow::sphereDEMSystemSemiImplicit::particleMaxId()
{
	return particles_->maxId();
}

void pFlow::sphereDEMSystemSemiImplicit::setParticleMaxId(uint32 id)
{
	implicit_->setMaxId(id);
}
