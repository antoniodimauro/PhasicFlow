/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Original file: PhasicFlowPlus/phasicFlowCoupling/couplingSystem/unresolved/momentumSphereUnresolvedCouplingSystem/momentumSphereUnresolvedCouplingSystem.C,
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

// Xiao & Sun (2011) Commun. Comput. Phys. 9, 297

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

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

        // J1: at a restart, the largest id already given the fluid velocity, written at that time by
        // commitNewParticleVelocity(); without the file (t = 0, older saves) every particle present counts as set
        if(couplingDEM_)
        {
            const Foam::Time& rt = this->cMesh().mesh().time();
            const Foam::fileName f = rt.globalPath()/rt.timeName()/"uniform"/"velFromFluidLastId";
            std::ifstream in(f.c_str());
            long long id = -1;
            if(in && (in >> id))
            {
                restartVelFromFluidId_    = static_cast<int64>(id);
                hasRestartVelFromFluidId_ = true;
                Foam::Info<< "Inserted-particle velocity: restart with the fluid velocity set up to id "
                    << static_cast<Foam::label>(id)
                    << " (" << f << ")" << Foam::endl;
            }
        }
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

// A particle centre on a processor face lies in a cell of both ranks (couplingMesh::pointInCell counts the face as
// inside), and every per-particle sum to the master (force, drag coefficient, sampled fluid velocity, host-cell
// extent) would count it twice. Each rank that located the particle adds 2^-rank; the sum keeps the lowest of these
// ranks in its leading bit, and only that rank stays the host, the rule of the sphereQuadrature holder. The other
// ranks keep the particle as one hosted elsewhere: the kernel still deposits on their cells. A particle that some
// rank holds but none locates is reported; one inserted in this step is held by no rank until the next mapping
// (16 Oct 2026, validation/pc-tests-2026-10-16/README.md).
void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::holdParticlesOnOneRank()
{
    if(!Plus::processor::isParallel()) return;

    // 2^-rank leaves the normal doubles near rank 1022
    if(Foam::Pstream::nProcs() > 1000)
    {
        fatalErrorInFunction<<"holdParticlesOnOneRank: more than 1000 ranks"<<endl;
        Plus::processor::abort(0);
    }

    auto& cell = const_cast<Plus::procCMField<Foam::label>&>(this->cMesh().parCellIndex());
    const int me = static_cast<int>(Foam::Pstream::myProcNo());

    // x: 2^-rank where this rank located the particle; y: 1 for every particle this rank holds
    Plus::realx3ProcCMField claim("rankClaim", zero3, this->centerMass(), false, false);   // temporary: no resize events
    for(size_t i=0; i<claim.size(); ++i)
    {
        const bool located = i < cell.size() && cell[i] >= 0;
        claim[i] = realx3(located ? static_cast<real>(std::ldexp(1.0, -me)) : real(0), real(1), real(0));
    }

    const size_t nAll = this->pDEMSystem().particlesCenterMassAllMaster().size();
    realx3Scratch_.assign(nAll, zero3);
    span<realx3> all(realx3Scratch_.data(), realx3Scratch_.size());
    auto mine = makeSpan(claim);
    auto& comm = this->parMapping().realx3ScatteredComm();

    if(!comm.collectSum(mine, all))
    {
        fatalErrorInFunction<<"holdParticlesOnOneRank: collectSum failed"<<endl;
        Plus::processor::abort(0);
    }

    if(Plus::processor::isMaster())
    {
        const auto ids = this->pDEMSystem().particleIdAllMaster();
        size_t nShared = 0;
        size_t nLost = 0;
        std::string lostIds;
        for(size_t g=0; g<nAll; ++g)
        {
            if(all[g].x() <= 0)
            {
                if(all[g].y() > 0 && ++nLost <= 10)
                {
                    lostIds += " " + std::to_string(g < ids.size() ? ids[g] : static_cast<uint32>(g));
                }
                continue;
            }
            int e = 0;
            if(std::frexp(static_cast<double>(all[g].x()), &e) != 0.5) ++nShared;
        }

        if(nShared > 0)
        {
            output<<Blue_Text("Particles on a processor face: ")<<Yellow_Text(nShared)<<
            ", each kept by its lowest rank only"<<endl;
        }
        if(nLost > 0)
        {
            WARNING<<nLost<<" particle(s) held by a rank but located in no processor mesh, ids"<<lostIds<<END_WARNING;
        }
    }

    if(!comm.distribute(all, mine))
    {
        fatalErrorInFunction<<"holdParticlesOnOneRank: distribute failed"<<endl;
        Plus::processor::abort(0);
    }

    for(size_t i=0; i<claim.size() && i<cell.size(); ++i)
    {
        if(cell[i] < 0) continue;
        int e = 0;
        std::frexp(static_cast<double>(claim[i].x()), &e);
        if(1 - e != me) cell[i] = -1;
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

    // C11: one kernel width per particle on every rank, from its host cell
    {
        const auto& pci = this->parCellIndex();
        const Foam::scalarField& ext = kernel->cellExtent();
        Plus::realProcCMField hostExt("hostCellExtent", real(0), this->centerMass(), false, false);   // temporary: no resize events
        for(size_t i=0; i<hostExt.size(); ++i)
        {
            if(i < pci.size() && pci[i] >= 0) hostExt[i] = static_cast<real>(ext[pci[i]]);
        }
        this->globalParticleSum(hostExt);
        kernel->setHostCellExtent(hostExt);
    }

    kernel->updateWeights(this->particleDiameter());

    kernel->localWeightSums(weightSums_);
    this->globalParticleSum(weightSums_);
    kernel->scaleWeights(weightSums_);
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::calculatePorosity()
{
    // update coupling mesh and map particles 
    this->cMesh().update();
    holdParticlesOnOneRank();

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

	if(couplingDEM_)
	{
		if(hasRestartVelFromFluidId_)
		{
			couplingDEM_->restoreVelFromFluidId(restartVelFromFluidId_);
			hasRestartVelFromFluidId_ = false;
		}
		couplingDEM_->setNewParticleVelocityFromFluid();
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::commitNewParticleVelocity()
{
	if(!setInsertedVelFromFluid_) return;

	if(couplingDEM_)
	{
		couplingDEM_->commitNewParticleVelocityFromFluid();

		// J1: the written DEM state can hold particles inserted in this step, still with the placeholder velocity;
		// the stored id tells the restart which ones still need the fluid velocity
		const Foam::Time& rt = this->cMesh().mesh().time();
		if(rt.writeTime())
		{
			const Foam::fileName f = rt.globalPath()/rt.timeName()/"uniform"/"velFromFluidLastId";
			Foam::mkDir(f.path());
			std::ofstream out(f.c_str());
			if(out)
			{
				out << static_cast<long long>(couplingDEM_->committedVelFromFluidId()) << "\n";
			}
			else
			{
				Foam::Warning<< "cannot write " << f << Foam::endl;
			}
		}
	}
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

	// C12: acceleration of the last completed step, by particle id, from the velocities of two step starts. A
	// particle without a previous record (first step after a start or a restart, or just inserted) gets zero: the
	// DEM acceleration field holds only the explicit forces, not the implicit drag and lift.
	// H1: with setInsertedVelocityFromFluid, a particle inserted in the last step still has its insertion velocity
	// here and is given the fluid velocity in this step's exchange; its velocity now is not recorded, so its next
	// acceleration is zero too instead of carrying the jump to the fluid velocity (QUESTIONS.md H1).
	if(Plus::processor::isMaster() && savedDEMTime_ > prevVelTime_)
	{
		auto ids = pds.particleIdAllMaster();
		const size_t n = std::min(ids.size(), vSpan.size());
		const int64 unmappedFloor = (setInsertedVelFromFluid_ && couplingDEM_)
			? couplingDEM_->unmappedIdFloor() : std::numeric_limits<int64>::max();

		std::unordered_map<uint32, realx3> now;
		now.reserve(2*n);
		accById_.clear();
		accById_.reserve(2*n);

		const real dtv = savedDEMTime_ - prevVelTime_;

		for(size_t i=0; i<n; ++i)
		{
			const uint32 id = ids[i];
			auto f = prevVelById_.find(id);
			if(prevVelTime_ >= 0 && f != prevVelById_.end())
			{
				accById_[id] = (vSpan[i] - f->second)/dtv;
			}
			if(static_cast<int64>(id) <= unmappedFloor) now[id] = vSpan[i];
		}

		prevVelById_.swap(now);
		prevVelTime_ = savedDEMTime_;
	}
}

void pFlow::coupling::momentumSemiImplicitSphereUnresolvedCouplingSystem::particleAcceleration
(
	Plus::realx3ProcCMField& acc
)const
{
	auto* self = const_cast<momentumSemiImplicitSphereUnresolvedCouplingSystem*>(this);
	auto& pds  = self->pDEMSystem();

	const size_t nAll = pds.particlesCenterMassAllMaster().size();
	accAll_.assign(nAll, zero3);

	if(Plus::processor::isMaster())
	{
		auto ids = pds.particleIdAllMaster();
		for(size_t i=0; i<nAll && i<ids.size(); ++i)
		{
			auto f = accById_.find(ids[i]);
			if(f != accById_.end()) accAll_[i] = f->second;
		}
	}

	span<realx3> all(accAll_.data(), accAll_.size());
	auto mine = makeSpan(acc);

	if(!self->parMapping().realx3ScatteredComm().distribute(all, mine))
	{
		fatalErrorInFunction<<"particleAcceleration: distribute failed"<<endl;
		Plus::processor::abort(0);
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
