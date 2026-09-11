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

#ifndef __momentumSemiImplicitSphereUnresolvedCouplingSystem_hpp__
#define __momentumSemiImplicitSphereUnresolvedCouplingSystem_hpp__


#include <vector>

#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "momentumInteractionImplicit.hpp"
#include "fluidFieldSampler.hpp"

namespace pFlow
{
class DEMSystemSemiImplicit;
}

namespace pFlow::coupling
{

class momentumSemiImplicitSphereUnresolvedCouplingSystem
:
	public unresolvedCouplingSystem
{
private:

	uniquePtr<porosity>					porosity_ = nullptr;

    momentumInteractionImplicit         momentumInteraction_;

	Timer 		porosityTimer_;

    bool 		requiresDistribution_ = false;

	DEMSystemSemiImplicit*			couplingDEM_ = nullptr;

	Plus::realProcCMField       dragCoeff_;

	Plus::realx3ProcCMField     dragFluidVel_;

	Plus::realx3ProcCMField     liftCrossVec_;

	Plus::realProcCMField       rotDragCoeff_;

	Plus::realProcCMField       addedMass_;

	Plus::realProcCMField       addedInertia_;

	Plus::realProcCMField       displacedMass_;

	Plus::realx3ProcCMField     fluidSpin_;

	Plus::realProcCMField       weightSums_;

	bool setInsertedVelFromFluid_ = false;

	bool setInsertedSpinFromFluid_ = false;

	uniquePtr<fluidFieldSampler> spinSampler_ = nullptr;

	mutable std::vector<real>	realScratch_;
	mutable std::vector<realx3>	realx3Scratch_;
	mutable std::vector<realx3>	stagingVec_;

	std::vector<realx3>			savedPos_;
	std::vector<realx3>			savedVel_;
	std::vector<realx3>			savedRVel_;
	std::vector<uint32>			savedNumInserted_;
	uint32						savedMaxId_ = static_cast<uint32>(-1);
	// DEM clock at saveDEMState(): a repeated exchange restarts the DEM from here,
	// not from the fluid time, so the sub-step remainder carries over
	real						savedDEMTime_ = -1.0;

	void updateKernelWeights();

	bool collectReal(Plus::realProcCMField& field, span<real> all, const char* what);

	bool collectRealx3(Plus::realx3ProcCMField& field, span<realx3> all, const char* what);

	void sendDragCoeffToDEM();

	void sendDragFluidVelToDEM();

	void sendLiftCrossVecToDEM();

	void sendRotDragCoeffToDEM();

	void sendAddedMassToDEM();

	void sendAddedInertiaToDEM();

	void sendDisplacedMassToDEM();

	void sendFluidSpinToDEM();

public:

    TypeInfo("sphereUnresolvedCouplingSystem<momentumSemiImplicit>");

    momentumSemiImplicitSphereUnresolvedCouplingSystem(
        word shapeTypeName,
        word couplingSystemType, 
        Foam::fvMesh& mesh,
        int argc, 
        char* argv[]);

    momentumSemiImplicitSphereUnresolvedCouplingSystem(const momentumSemiImplicitSphereUnresolvedCouplingSystem&) = delete;

    momentumSemiImplicitSphereUnresolvedCouplingSystem& operator=(const momentumSemiImplicitSphereUnresolvedCouplingSystem&) = delete;

    momentumSemiImplicitSphereUnresolvedCouplingSystem(momentumSemiImplicitSphereUnresolvedCouplingSystem&&) = delete;

    momentumSemiImplicitSphereUnresolvedCouplingSystem& operator=(momentumSemiImplicitSphereUnresolvedCouplingSystem&&) = delete;

    ~momentumSemiImplicitSphereUnresolvedCouplingSystem() override = default;

    add_vCtor
    (
        unresolvedCouplingSystem,
        momentumSemiImplicitSphereUnresolvedCouplingSystem,
        word
    );

    void calculatePorosity() override;

    void calculateMomentumCoupling() override;

    void calculateHeatCoupling() override;

    void calculateMassCoupling() override;

    const Foam::volScalarField& alpha()const override
    {
        return porosity_().alpha();
    }

    Foam::tmp<Foam::volScalarField> Sp()const override;

    Foam::tmp<Foam::volVectorField> Su()const override;
    
    Foam::tmp<Foam::volScalarField> heatSource() const override
    {
        notImplementedFunction;
        return Foam::tmp<Foam::volScalarField>(nullptr);
    }

    Foam::tmp<Foam::volScalarField> massSource(const word& specieName)const  override
    {
        notImplementedFunction;
        return Foam::tmp<Foam::volScalarField>(nullptr);
    }

    word shapeTypeName() const override
    {
        return "sphere";
    }

    word couplingSystemType()const override
    {
        return "momentumSemiImplicit";
    }

    bool requireCellDistribution()const override
    {
        return requiresDistribution_;
    }

	bool sendDataToDEM(real t, real dt) override;

	inline
	DEMSystemSemiImplicit* couplingDEM()const
	{
		return couplingDEM_;
	}

	inline
	Plus::procDEMSystem& procDEM()const
	{
		return const_cast<momentumSemiImplicitSphereUnresolvedCouplingSystem*>(this)->pDEMSystem();
	}

	bool globalParticleSum(Plus::realProcCMField& field)const;

	bool globalParticleSum(Plus::realx3ProcCMField& field)const;

	bool globalParticleSum(Plus::procCMField<Foam::vector>& field)const;

	void commitHistoryForces()
	{
		momentumInteraction_.commitHistoryForces();
	}

	inline
	bool setInsertedVelFromFluid()const
	{
		return setInsertedVelFromFluid_;
	}

	void applyFluidVelocityToNewParticles();

	void commitNewParticleVelocity();

	void saveDEMState();

	void restoreDEMState(real tStart);

	size_t numParticlesMaster();

	inline
	size_t numSavedParticlesMaster()const
	{
		return savedPos_.size();
	}
};

const momentumSemiImplicitSphereUnresolvedCouplingSystem&
semiImplicitCoupling(const unresolvedCouplingSystem& uCS);

}

#endif
