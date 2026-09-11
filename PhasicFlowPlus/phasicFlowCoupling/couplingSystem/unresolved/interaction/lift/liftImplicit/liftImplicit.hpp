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

// Zeng et al. (2009) Phys. Fluids 21, 033302

#ifndef __liftImplicit_hpp__
#define __liftImplicit_hpp__

#include "lift.hpp"
#include "fluidFieldSampler.hpp"
#include "wallLift.hpp"
#include "surfaceRotationTorqueImplicit.hpp"

namespace pFlow::coupling
{

class distributionBase;

class dragImplicit;

class liftImplicit
:
    public lift
{
private:

    fluidFieldSampler           liftUSampler_;
    fluidFieldSampler           curlUSampler_;

    bool                        wallCorr_     = false;
    Foam::scalar                minWallRadii_ = 2.0;

    mutable Foam::scalar        reportedWallLift_ = 0.0;
    mutable bool                notedShearLift_   = false;

    mutable bool                wallLiftApplied_  = false;

    mutable Plus::realx3ProcCMField liftScratch_;

    uniquePtr<surfaceRotationTorqueImplicit> surfTorqueImplicit_;

    const dragImplicit*         drag_ = nullptr;

    // liftReaction no: the lift is not returned to the fluid. The lift closures
    // give the total lift of a sphere, its own disturbance included; returned
    // to the fluid, the lift drives a flow that the particle samples again
    // (strong in a rotating fluid at small Re_p), so the self-induced lift
    // feedback is not applied either
    bool                        liftReaction_ = true;

protected:

    const Plus::procCMField<Foam::vector>&
    sampleFluidVel(const Foam::volVectorField& U)
    {
        return liftUSampler_.sample(U);
    }

    const Plus::procCMField<Foam::vector>&
    sampleVorticity(const Foam::volVectorField& curlU)
    {
        return curlUSampler_.sample(curlU);
    }

    Plus::realx3ProcCMField& liftScratch()const { return liftScratch_; }

    Foam::vector fluidVelForLift(size_t i, const Foam::vector& sampled)const;

public:

    TypeInfo("liftImplicit");

    void setDrag(const dragImplicit* d);

    bool liftReaction()const { return liftReaction_; }

    const Plus::realx3ProcCMField& liftPerParticle()const { return liftScratch_; }

    liftImplicit(
        const unresolvedCouplingSystem& uCS,
        const porosity& 				prsty);

    virtual ~liftImplicit() = default;

    bool requireCellDistribution()const
    {
        return liftUSampler_.requireCellDistribution()
            || curlUSampler_.requireCellDistribution()
            || surfTorqueImplicit_().requireCellDistribution();
    }

    static Foam::vector wallLiftForce
    (
        const Foam::vector& nw,
        const Foam::vector& uSlip,
        const Foam::scalar  dp,
        const Foam::scalar  rho,
        const Foam::scalar  nu,
        const Foam::scalar  hOverA
    )
    {
        const Foam::vector uPar = uSlip - (uSlip & nw)*nw;

        const Foam::scalar V = Foam::mag(uPar);

        if(V <= Foam::SMALL) return Foam::vector::zero;

        const Foam::scalar Ret = V*dp/nu;
        const Foam::scalar L   = 0.5*hOverA;

        const Foam::scalar cL = zengWallLift(L, Ret);

        return (Foam::constant::mathematical::pi/8.0)*rho*dp*dp*V*V*cL*nw;
    }

    void reportWallLift(Foam::scalar fractionObserved)const;

    bool wallCorrection()const { return wallCorr_; }

    Foam::scalar minWallRadii()const { return minWallRadii_; }

    bool wallLiftApplied()const { return wallLiftApplied_; }

    bool wallLiftActive()const { return wallCorr_; }

    void calculateLiftForceTorqueImplicit(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        const distributionBase&         cellDistribution,
        Plus::realx3ProcCMField&        particleForce,
        Plus::realx3ProcCMField&        particleTorque,
        Plus::realx3ProcCMField&        liftCrossVec,
        Plus::realProcCMField&          rotDragCoeff,
        Plus::realProcCMField&          addedInertia);

    void commitRotationalHistory()
    {
        surfTorqueImplicit_->commitHistory();
    }

    bool rotationalHistoryActive()const
    {
        return surfTorqueImplicit_->historyActive();
    }

    bool rotationalHistoryEvaluated()const
    {
        return surfTorqueImplicit_->historyEvaluated();
    }

    virtual
    void calculateLiftForceImplicit(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        const distributionBase&         cellDistribution,
        Plus::realx3ProcCMField&        particleForce,
        Plus::realx3ProcCMField&        particleTorque,
        Plus::realx3ProcCMField&        liftCrossVec) = 0;

    void calculateLiftForce(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleForce,
        Plus::realx3ProcCMField&        particleTorque) override;

    static
    uniquePtr<liftImplicit> create
    (
        const unresolvedCouplingSystem& uCS,
        const porosity& 				prsty
    );
};

}

#endif
