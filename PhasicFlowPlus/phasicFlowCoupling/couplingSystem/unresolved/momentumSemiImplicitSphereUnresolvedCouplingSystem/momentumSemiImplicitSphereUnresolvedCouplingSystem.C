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

// Xiao & Sun (2011) Commun. Comput. Phys. 9, 297

#include <algorithm>

#include "momentumSemiImplicitSphereUnresolvedCouplingSystem.hpp"
#include "DEMSystemSemiImplicit.hpp"
#include "particleGaussian.hpp"
#include "fvcCurl.H"

pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::momentumSemiImplicitSphereUnresolvedCouplingSystem
(
    word shapeTypeName,
    word couplingSystemType, 
    Foam::fvMesh& mesh,
    int argc, 
    char* argv[]
)
:
    unresolvedCouplingSystem
    (
        "sphereSemiImplicit", 
        couplingSystemType, 
        mesh, 
        argc, 
        argv
    ),
    porosity_
    (
        porosity::create
        (
            *this, 
            this->cMesh(),
            this->particleDiameter()
        )
    ),
    momentumInteraction_
    (
        *this, 
        porosity_()
    ),
    porosityTimer_
    (
        "porosity", 
        &this->couplingTimers()
    ),
    dragCoeff_("dragCoeff", this->centerMass()),
    dragFluidVel_("dragFluidVel", this->centerMass()),
    liftCrossVec_("liftCrossVec", this->centerMass()),
    rotDragCoeff_("rotDragCoeff", this->centerMass()),
    addedMass_("addedMass", this->centerMass()),
    addedInertia_("addedInertia", this->centerMass()),
    displacedMass_("displacedMass", this->centerMass()),
    fluidSpin_("fluidSpin", this->centerMass()),
    weightSums_("kernelWeightSum", real(0), this->centerMass())
{
    requiresDistribution_ = 
        porosity_().requireCellDistribution()|| momentumInteraction_.requireCellDistribution();

    if(auto* dem = this->pDEMSystem().demSystemPtr())
    {
        couplingDEM_ = dynamic_cast<DEMSystemSemiImplicit*>(dem);
        if(!couplingDEM_)
        {
            fatalErrorInFunction
                << "the semi-implicit coupling needs the sphereSemiImplicitDEMSystem, found "
                << dem->typeName() << endl;
            Plus::processor::abort(0);
        }
    }

    setInsertedVelFromFluid_ =
        this->getOrDefault<Foam::Switch>(
            "setInsertedVelocityFromFluid", Foam::Switch(false));

    if(setInsertedVelFromFluid_)
    {
        Foam::Info<< Blue_Text(
            "Inserted particles will be initialised with the local fluid "
            "velocity (setInsertedVelocityFromFluid = yes).") << Foam::endl;
    }

    setInsertedSpinFromFluid_ =
        this->getOrDefault<Foam::Switch>(
            "setInsertedSpinFromFluid", Foam::Switch(false));

    if(setInsertedSpinFromFluid_)
    {
        if(!setInsertedVelFromFluid_)
        {
            fatalErrorInFunction
                << "setInsertedSpinFromFluid needs setInsertedVelocityFromFluid yes"
                << endl;
            Plus::processor::abort(0);
        }

        spinSampler_ = makeUnique<fluidFieldSampler>(*this, "insertedSpin");
        requiresDistribution_ =
            requiresDistribution_ || spinSampler_->requireCellDistribution();

        if(couplingDEM_) couplingDEM_->setSpinFromFluid(true);

        Foam::Info<< Blue_Text(
            "Inserted particles will be initialised with half the local fluid "
            "vorticity as spin (setInsertedSpinFromFluid = yes).") << Foam::endl;
    }
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::updateKernelWeights()
{
    auto& dist = const_cast<distributionBase&>(this->distribution());

    auto* kernel = dynamic_cast<particleGaussian*>(&dist);

    if(!kernel)
    {
        this->updateDistributionWeights();
        return;
    }

    Foam::scalar dMax = 0;
    for(size_t i=0; i<this->particleDiameter().size(); ++i)
    {
        dMax = Foam::max(dMax,
            static_cast<Foam::scalar>(this->particleDiameter()[i]));
    }
    Foam::reduce(dMax, Foam::maxOp<Foam::scalar>());
    kernel->setShippingRadius(
        this->parMapping().domainExpansionRatio()*dMax);

    kernel->updateWeights(this->particleDiameter());

    kernel->localWeightSums(weightSums_);
    this->globalParticleSum(weightSums_);
    kernel->scaleWeights(weightSums_);
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::calculatePorosity()
{
    // update coupling mesh and map particles 
    this->cMesh().update();

    porosityTimer_.start();

    // update weights for distribution
    if(requiresDistribution_)
        this->updateKernelWeights();

    // calculate porosity 
    porosity_->calculatePorosity();

    porosityTimer_.end();

    Foam::Info<<Blue_Text("Porosity time: ")<< 
                Yellow_Text(porosityTimer_.lastTime())<<
                Yellow_Text(" s")<<Foam::endl;
}


void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::calculateMomentumCoupling()
{
    const auto& U = this->cMesh().mesh().template lookupObject<Foam::volVectorField>("U");

    const auto& vp = this->particleVelocity();

    const auto& wp = this->particleRVelocity();

    auto& fluidForce = this->fluidForce();

    auto& fluidTorque = this->fluidTorque();

    std::fill(fluidForce.begin(), fluidForce.end(), 0.0);
    std::fill(fluidTorque.begin(), fluidTorque.end(), 0.0);
    std::fill(dragCoeff_.begin(), dragCoeff_.end(), 0.0);
    std::fill(dragFluidVel_.begin(), dragFluidVel_.end(), realx3(0.0));
    std::fill(liftCrossVec_.begin(), liftCrossVec_.end(), realx3(0.0));
    std::fill(rotDragCoeff_.begin(), rotDragCoeff_.end(), 0.0);
    std::fill(addedMass_.begin(), addedMass_.end(), 0.0);
    std::fill(addedInertia_.begin(), addedInertia_.end(), 0.0);
    std::fill(displacedMass_.begin(), displacedMass_.end(), 0.0);

    momentumInteraction_.calculateCoupling(
        U, vp, wp, fluidForce, fluidTorque,
        dragCoeff_, dragFluidVel_, liftCrossVec_, rotDragCoeff_,
        addedMass_, addedInertia_, displacedMass_);

    if(setInsertedSpinFromFluid_)
    {
        const Foam::tmp<Foam::volVectorField> tCurlU = Foam::fvc::curl(U);
        const auto& curlU = spinSampler_->sample(tCurlU());

        for(size_t i=0; i<fluidSpin_.size(); ++i)
        {
            fluidSpin_[i] = realx3(
                0.5*curlU[i].x(), 0.5*curlU[i].y(), 0.5*curlU[i].z());
        }
    }
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::calculateHeatCoupling()
{
    notImplementedFunction;
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::calculateMassCoupling()
{
    notImplementedFunction;
}

Foam::tmp<Foam::volScalarField> 
    pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::Sp()const
{
    return Foam::tmp<Foam::volScalarField>(momentumInteraction_.Sp());
}

Foam::tmp<Foam::volVectorField> 
    pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::Su()const
{

    const auto& SU = momentumInteraction_.Su();
    auto tmpLift = momentumInteraction_.liftForce();
    const auto& lift = tmpLift.ref();
    
    auto tmpVirtualMass = momentumInteraction_.virtualMassForce();
    const auto& virtualMass = tmpVirtualMass.ref();
    
    auto SUall = Foam::tmp<Foam::volVectorField>::New
    (
        Foam::IOobject
        (
            "SUall",
            Foam::timeName(this->cMesh().mesh().time()),
            this->cMesh().mesh(),
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE,
            false
        ),
        SU + lift + virtualMass
    );
    return SUall;
}

bool pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::sendDataToDEM(real, real)
{
	this->sendDataTimer().start();
		this->sendFluidForceToDEM();
		this->sendFluidTorqueToDEM();
		sendDragCoeffToDEM();
		sendDragFluidVelToDEM();
		sendLiftCrossVecToDEM();
		sendRotDragCoeffToDEM();
		sendAddedMassToDEM();
		sendAddedInertiaToDEM();
		sendDisplacedMassToDEM();
		if(setInsertedSpinFromFluid_) sendFluidSpinToDEM();
	this->sendDataTimer().end();
	return true;
}

bool pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::collectReal
(
	Plus::realProcCMField& field,
	span<real> all,
	const char* what
)
{
	for(uint32 i=0; i<all.size(); i++)
		all[i] = 0;
	auto mine = makeSpan(field);
	if(!this->parMapping().realScatteredComm().collectSum(mine, all))
	{
		fatalErrorInFunction<<
		"Failed to perform collective sum over processors for "<<what<<endl;
		Plus::processor::abort(0);
	}
	return true;
}

bool pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::collectRealx3
(
	Plus::realx3ProcCMField& field,
	span<realx3> all,
	const char* what
)
{
	for(uint32 i=0; i<all.size(); i++)
		all[i] = zero3;
	auto mine = makeSpan(field);
	if(!this->parMapping().realx3ScatteredComm().collectSum(mine, all))
	{
		fatalErrorInFunction<<
		"Failed to perform collective sum over processors for "<<what<<endl;
		Plus::processor::abort(0);
	}
	return true;
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::sendDragCoeffToDEM()
{
	collectReal(dragCoeff_, couplingDEM_ ? couplingDEM_->parDragCoeff() : span<real>(), "drag coefficient");

	if(couplingDEM_ && !couplingDEM_->sendDragCoeffToDEM())
	{
		fatalErrorInFunction<< "could not perform sendDragCoeffToDEM"<<endl;
		Plus::processor::abort(0);
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::sendDragFluidVelToDEM()
{
	collectRealx3(dragFluidVel_, couplingDEM_ ? couplingDEM_->parDragFluidVel() : span<realx3>(), "drag fluid velocity");

	if(couplingDEM_ && !couplingDEM_->sendDragFluidVelToDEM())
	{
		fatalErrorInFunction<< "could not perform sendDragFluidVelToDEM"<<endl;
		Plus::processor::abort(0);
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::sendLiftCrossVecToDEM()
{
	collectRealx3(liftCrossVec_, couplingDEM_ ? couplingDEM_->parLiftCrossVec() : span<realx3>(), "lift cross-product vector");

	if(couplingDEM_ && !couplingDEM_->sendLiftCrossVecToDEM())
	{
		fatalErrorInFunction<< "could not perform sendLiftCrossVecToDEM"<<endl;
		Plus::processor::abort(0);
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::sendRotDragCoeffToDEM()
{
	collectReal(rotDragCoeff_, couplingDEM_ ? couplingDEM_->parRotDragCoeff() : span<real>(), "rotational drag coefficient");

	if(couplingDEM_ && !couplingDEM_->sendRotDragCoeffToDEM())
	{
		fatalErrorInFunction<< "could not perform sendRotDragCoeffToDEM"<<endl;
		Plus::processor::abort(0);
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::sendAddedMassToDEM()
{
	collectReal(addedMass_, couplingDEM_ ? couplingDEM_->parAddedMass() : span<real>(), "added mass");

	if(couplingDEM_ && !couplingDEM_->sendAddedMassToDEM())
	{
		fatalErrorInFunction<< "could not perform sendAddedMassToDEM"<<endl;
		Plus::processor::abort(0);
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::sendAddedInertiaToDEM()
{
	collectReal(addedInertia_, couplingDEM_ ? couplingDEM_->parAddedInertia() : span<real>(), "added inertia");

	if(couplingDEM_ && !couplingDEM_->sendAddedInertiaToDEM())
	{
		fatalErrorInFunction<< "could not perform sendAddedInertiaToDEM"<<endl;
		Plus::processor::abort(0);
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::sendFluidSpinToDEM()
{
	collectRealx3(fluidSpin_, couplingDEM_ ? couplingDEM_->parFluidSpin() : span<realx3>(), "fluid spin");

	if(couplingDEM_ && !couplingDEM_->sendFluidSpinToDEM())
	{
		fatalErrorInFunction<< "could not perform sendFluidSpinToDEM"<<endl;
		Plus::processor::abort(0);
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::sendDisplacedMassToDEM()
{
	collectReal(displacedMass_, couplingDEM_ ? couplingDEM_->parDisplacedMass() : span<real>(), "displaced mass");

	if(couplingDEM_ && !couplingDEM_->sendDisplacedMassToDEM())
	{
		fatalErrorInFunction<< "could not perform sendDisplacedMassToDEM"<<endl;
		Plus::processor::abort(0);
	}
}

bool pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::globalParticleSum
(
	Plus::realProcCMField& field
)const
{
	if(!Plus::processor::isParallel()) return true;

	auto* self = const_cast<momentumSemiImplicitSphereUnresolvedCouplingSystem*>(this);

	const size_t nAll = self->pDEMSystem().particlesCenterMassAllMaster().size();
	realScratch_.assign(nAll, 0);
	span<real> all(realScratch_.data(), realScratch_.size());

	auto mine = makeSpan(field);

	if(!self->parMapping().realScatteredComm().collectSum(mine, all))
	{
		fatalErrorInFunction<<"globalParticleSum: collectSum failed"<<endl;
		Plus::processor::abort(0);
		return false;
	}

	if(!self->parMapping().realScatteredComm().distribute(all, mine))
	{
		fatalErrorInFunction<<"globalParticleSum: distribute failed"<<endl;
		Plus::processor::abort(0);
		return false;
	}

	return true;
}

bool pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::globalParticleSum
(
	Plus::realx3ProcCMField& field
)const
{
	if(!Plus::processor::isParallel()) return true;

	auto* self = const_cast<momentumSemiImplicitSphereUnresolvedCouplingSystem*>(this);

	const size_t nAll = self->pDEMSystem().particlesCenterMassAllMaster().size();
	realx3Scratch_.assign(nAll, zero3);
	span<realx3> all(realx3Scratch_.data(), realx3Scratch_.size());

	auto mine = makeSpan(field);

	if(!self->parMapping().realx3ScatteredComm().collectSum(mine, all))
	{
		fatalErrorInFunction<<"globalParticleSum: collectSum failed"<<endl;
		Plus::processor::abort(0);
		return false;
	}

	if(!self->parMapping().realx3ScatteredComm().distribute(all, mine))
	{
		fatalErrorInFunction<<"globalParticleSum: distribute failed"<<endl;
		Plus::processor::abort(0);
		return false;
	}

	return true;
}

bool pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::globalParticleSum
(
	Plus::procCMField<Foam::vector>& field
)const
{
	if(!Plus::processor::isParallel()) return true;

	auto* self = const_cast<momentumSemiImplicitSphereUnresolvedCouplingSystem*>(this);

	const size_t n = field.size();
	stagingVec_.resize(n);
	for(size_t i=0; i<n; ++i)
	{
		stagingVec_[i] = realx3(field[i].x(), field[i].y(), field[i].z());
	}

	const size_t nAll = self->pDEMSystem().particlesCenterMassAllMaster().size();
	realx3Scratch_.assign(nAll, zero3);
	span<realx3> all(realx3Scratch_.data(), realx3Scratch_.size());
	span<realx3> mine(stagingVec_.data(), stagingVec_.size());

	auto& comm = self->parMapping().realx3ScatteredComm();

	if(!comm.collectSum(mine, all) || !comm.distribute(all, mine))
	{
		fatalErrorInFunction<<"globalParticleSum(vector) failed"<<endl;
		Plus::processor::abort(0);
		return false;
	}

	for(size_t i=0; i<n; ++i)
	{
		field[i] = Foam::vector
		(
			stagingVec_[i].x(), stagingVec_[i].y(), stagingVec_[i].z()
		);
	}

	return true;
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::applyFluidVelocityToNewParticles()
{
	if(!setInsertedVelFromFluid_) return;

	if(couplingDEM_) couplingDEM_->setNewParticleVelocityFromFluid();
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::commitNewParticleVelocity()
{
	if(!setInsertedVelFromFluid_) return;

	if(couplingDEM_) couplingDEM_->commitNewParticleVelocityFromFluid();
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::saveDEMState()
{
	auto& pds = this->pDEMSystem();

	pds.getDataFromDEM();

	auto pSpan  = pds.particlesCenterMassAllMaster();
	auto vSpan  = pds.particlesVelocityAllMaster();
	auto rvSpan = pds.particlesRVelocityAllMaster();
	savedPos_.assign (pSpan.begin(),  pSpan.end());
	savedVel_.assign (vSpan.begin(),  vSpan.end());
	savedRVel_.assign(rvSpan.begin(), rvSpan.end());

	if(couplingDEM_)
	{
		savedNumInserted_ = couplingDEM_->insertionNumInserted();
		savedMaxId_       = couplingDEM_->particleMaxId();
		couplingDEM_->saveIntegrationHistory();
	}

	savedDEMTime_ = -1.0;
	if(auto* dem = pds.demSystemPtr())
	{
		savedDEMTime_ = dem->Control().time().currentTime();
	}
	else
	{
		savedNumInserted_.clear();
		savedMaxId_ = static_cast<uint32>(-1);
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::restoreDEMState(real tStart)
{
	if(!couplingDEM_) return;

	couplingDEM_->restoreParticleKinematics(savedPos_, savedVel_, savedRVel_);
	couplingDEM_->setInsertionNumInserted(savedNumInserted_);
	couplingDEM_->setParticleMaxId(savedMaxId_);
	// The DEM steps with its own dt and stops within half a step of the fluid time,
	// so its clock at the start of the fluid step differs from tStart by up to half
	// a DEM step. Restarting from tStart dropped that remainder on every repeated
	// exchange: with a fluid step of 1.46 DEM steps a particle received 1/1.46 of
	// the exchanged momentum (5 Oct 2026). Restart from the saved DEM clock.
	couplingDEM_->setControlTime(savedDEMTime_ >= 0 ? savedDEMTime_ : tStart);
	couplingDEM_->restoreIntegrationHistory();
}

size_t pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::numParticlesMaster()
{
	return this->pDEMSystem().particlesCenterMassAllMaster().size();
}

const pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem&
pFlow::coupling::semiImplicitCoupling(const unresolvedCouplingSystem& uCS)
{
	auto* sc = dynamic_cast<const momentumSemiImplicitSphereUnresolvedCouplingSystem*>(&uCS);

	if(!sc)
	{
		fatalErrorInFunction
			<< "this model needs couplingSystemType momentumSemiImplicit, found "
			<< uCS.couplingSystemType() << endl;
		Plus::processor::abort(0);
	}

	return *sc;
}
