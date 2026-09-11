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

// Ireland & Desjardins (2017) J. Comput. Phys. 338, 405; Balachandar, Liu & Lakhote (2019) J. Comput. Phys. 376, 160; Maxey & Riley (1983) Phys. Fluids 26, 883; van Hinsberg et al. (2011) J. Comput. Phys. 230, 1465; Brenner (1961) Chem. Eng. Sci. 16, 242; Goldman, Cox & Brenner (1967) Chem. Eng. Sci. 22, 637, 653

#ifndef __dragImplicit_hpp__
#define __dragImplicit_hpp__

#include "drag.hpp"
#include "fluidFieldSampler.hpp"
#include "selfInduction.hpp"
#include "basset.hpp"
#include "pairInteraction.hpp"
#include "sphereQuadrature.hpp"
#include "wallCorrection.hpp"

#include <vector>

namespace pFlow::coupling
{

class dragImplicit
:
    public drag
{
private:

    bool                        undisturbedIsMaterialDerivative_ = false;

    fluidFieldSampler           undisturbedSampler_;

protected:


    bool                        selfInduced_ = true;

    Foam::scalar                maxBeta_ = 0.95;

    static constexpr Foam::scalar selfInducedC_ = selfInducedCentreC;

    Foam::scalar                selfInducedSampleFactor_ = 1.0;

    mutable Foam::scalar        reportedBeta_ = 0.0;

    mutable bool                warnedNoWidth_ = false;

    mutable bool                reportedRegime_ = false;

    mutable Foam::scalar        reportedFaxen_ = -1.0;

    mutable Foam::tmp<Foam::volVectorField> tDDtU_;

    mutable sphereQuadrature        quad_;

    bool faxen_ = false;

    bool wallCorr_ = false;

    // Brenner (1961); Goldman, Cox & Brenner (1967); Ryu & Matsudaira (2010); Zeng et al. (2009)
    bool wallFiniteRe_ = false;

    // Brenner (1961); Vasseur & Cox (1977)
    bool wallGapRe_ = false;

    // Vasseur & Cox (1977) Fig. 4
    Foam::scalar wallGapReScale_ = 1.25;

    // wallNormalImplicit; Izard, Bonometti & Lacaze (2014) eq. (2.18)
    bool wallNormalImplicit_ = false;

    mutable bool warnedWallRe_ = false;

    Foam::scalar minWallRadii_ = 2.0;

    Foam::scalar minBlobWallRadii_ = 2.0;

    mutable Foam::scalar reportedWall_ = 1.0;

    mutable Foam::scalar reportedWallBeta_ = 1.0;

    mutable bool         warnedUnsteadyWall_ = false;

    mutable bool         warnedLubrication_  = false;
    mutable bool         warnedBlobWall_     = false;

    mutable Plus::realProcCMField   spScratch_;
    mutable Plus::realx3ProcCMField upScratch_;

    mutable Plus::realx3ProcCMField fbScratch_;

    // Capecelatro & Desjardins (2013) eq. (28)
    bool conservativeReaction_ = false;

    // fluid receives the drag as a regularised point force K(x) F
    // (Maxey & Patel 2001): first moment of Sp (U - u_p) removed
    bool pointForceReaction_ = false;

    mutable bool reportedConservative_ = false;

    // Esmaily & Horwitz (2018) eq. (2); Balachandar et al. (2019)
    bool selfInducedLift_ = false;

    // Gotoh (1990); Candelier, Mehlig & Magnaudet (2019) eq. (4.10): drag of a
    // sphere in the rotation part of the local flow, and the same screening of
    // the self-induced velocity (selfInduction.hpp)
    bool rotationDrag_ = false;

    // rotationDragCombination quadrature: inertia and rotation both thin the
    // boundary layer, 1/delta^2 = 1/delta_I^2 + 1/delta_E^2, so the two
    // increases of the drag over Stokes add in quadrature (Davis 1992, J. Fluid
    // Mech. 237, 13, Table 2: within 6 % of the exact small-Re result)
    bool rotationDragQuadrature_ = false;

    // rotationDragGeostrophic Stewartson1953: the drag of a sphere moving across
    // the axis of a rapidly rotating fluid, 0.49 x 2 Omega rho V U (Stewartson
    // 1953; Mason 1975 measured 0.5), i.e. (0.49*4/9) Ta over the Stokes drag,
    // added to the rotational drag with the weight Ta/(1 + Ta) = Omega/(Omega +
    // nu/a^2): the geostrophic flow forms at the rotation rate and is diffused
    // at the viscous rate over the sphere
    bool rotationDragGeostrophic_ = false;

    // rotationSelfInducedLateral yes: in a rotating fluid the self-induced
    // velocity of the deposited drag also has a part normal to the slip. The
    // resistance of a sphere in a rotating fluid at small Re is
    // I + Ta^1/2 [[0.524, -0.0517], [0.0517, 0.524]] (Gotoh 1990; Candelier,
    // Mehlig & Magnaudet 2019, J. Fluid Mech. 864, eq. 4.10); for the kernel-
    // equivalent sphere (a_s^2 = pi sigma^2) the lateral part of the mobility
    // over its in-line part is 0.0517 T/(1 + 0.524 T), T = Ta_s^1/2, along
    // omega_hat x slip. It is removed from the sampled slip with the in-line part.
    bool rotationSelfInducedLateral_ = false;

    uniquePtr<fluidFieldSampler> rotVorticitySampler_;

    uniquePtr<fluidFieldSampler> rotStrainSampler_;

    mutable Foam::scalar reportedRotationFactor_ = 0.0;

    void reportRotationDrag
    (
        Foam::scalar factorObserved,
        Foam::scalar taObserved,
        Foam::scalar psiRatioObserved,
        Foam::scalar lateralObserved = 0
    )const;

    mutable std::vector<Foam::vector> liftFeedback_;
    mutable std::vector<Foam::vector> liftFeedbackPos_;
    mutable std::vector<Foam::vector> undisturbedVel_;
    mutable std::vector<char>         undisturbedSet_;

    mutable Plus::realx3ProcCMField ueScratch_;
    mutable Plus::realx3ProcCMField wxScratch_;
    mutable Plus::realx3ProcCMField ubScratch_;

    basset                          basset_;

    pairInteraction                 pair_;

    void reportSelfInduction(Foam::scalar betaObserved, Foam::label anyWidth)const;

    void reportFaxen(Foam::scalar ratioObserved)const;

    void reportDisturbanceRegime(
        Foam::scalar dp,
        Foam::scalar sigma)const;

    Plus::realProcCMField&   spScratch()const   { return spScratch_; }
    Plus::realx3ProcCMField& upScratch()const   { return upScratch_; }
    Plus::realx3ProcCMField& fbScratch()const   { return fbScratch_; }
    Plus::realx3ProcCMField& ueScratch()const   { return ueScratch_; }
    Plus::realx3ProcCMField& wxScratch()const   { return wxScratch_; }
    Plus::realx3ProcCMField& ubScratch()const   { return ubScratch_; }

    void reportConservative(Foam::scalar corrOverDrag)const;

    basset& bassetModel()const { return const_cast<basset&>(basset_); }

    pairInteraction& pairModel()const { return const_cast<pairInteraction&>(pair_); }

    inline Foam::scalar selfInducedSampleFactor()const
    {
        return selfInducedSampleFactor_;
    }

    const Plus::procCMField<Foam::vector>&
    sampleUndisturbed(const Foam::volVectorField& field)
    {
        return undisturbedSampler_.sample(field);
    }

    // vorticity and strain rate sqrt(2 S:S) (x component) at the particles, for
    // the rotation part |omega| - sigma of the local flow
    const Plus::procCMField<Foam::vector>&
    sampleRotationVorticity(const Foam::volVectorField& field)
    {
        return rotVorticitySampler_().sample(field);
    }

    const Plus::procCMField<Foam::vector>&
    sampleRotationStrain(const Foam::volVectorField& field)
    {
        return rotStrainSampler_().sample(field);
    }


public:

    TypeInfo("dragImplicit");

    dragImplicit(
        const unresolvedCouplingSystem& uCS,
        const porosity&                 prsty);

    virtual ~dragImplicit() = default;

    const Foam::volVectorField& materialDerivativeU
    (
        const Foam::volVectorField& U
    )const;

    void clearMaterialDerivative()const
    {
        tDDtU_.clear();
    }

    // Maxey & Riley (1983)
    Foam::tmp<Foam::volVectorField>
    pressureGradient(const Foam::volScalarField& rho)const;

    bool wallCorrection()const { return wallCorr_; }

    bool conservativeReaction()const { return conservativeReaction_; }

    bool pointForceReaction()const { return pointForceReaction_; }

    bool selfInducedLift()const { return selfInducedLift_; }

    const pairInteraction* pairModelPtr()const { return &pair_; }

    bool rotationDrag()const { return rotationDrag_; }

    bool rotationDragQuadrature()const { return rotationDragQuadrature_; }

    bool rotationDragGeostrophic()const { return rotationDragGeostrophic_; }

    bool rotationSelfInducedLateral()const { return rotationSelfInducedLateral_; }

    void setLiftFeedback
    (
        const Plus::realx3ProcCMField& lift,
        const Plus::centerMassField&   centres
    )const;

    bool liftFeedbackFor
    (
        size_t              i,
        const Foam::vector& centre,
        Foam::scalar        radius,
        Foam::vector&       fL
    )const;

    bool undisturbedVelocity(size_t i, Foam::vector& u)const;

    void resizeUndisturbed(size_t n)const
    {
        undisturbedVel_.assign(n, Foam::vector::zero);
        undisturbedSet_.assign(n, 0);
    }

    void setUndisturbed(size_t i, const Foam::vector& u)const
    {
        undisturbedVel_[i] = u;
        undisturbedSet_[i] = 1;
    }

    bool wallFiniteRe()const { return wallFiniteRe_; }

    bool wallGapRe()const { return wallGapRe_; }

    Foam::scalar wallGapReScale()const { return wallGapReScale_; }

    bool wallNormalImplicit()const { return wallNormalImplicit_; }

    void reportWallRe(Foam::scalar ReMaxNearWall)const;

    Foam::scalar minWallRadii()const { return minWallRadii_; }

    Foam::scalar minBlobWallRadii()const { return minBlobWallRadii_; }

    void reportWallSelfInduction(
        Foam::scalar gObserved,
        Foam::scalar minHOverAEff)const;

    void reportWall(Foam::scalar factorObserved, Foam::scalar minHoverA)const;

    bool faxen()const
    {
        return faxen_;
    }

    sphereQuadrature& quadrature()const { return quad_; }

    bool bassetActive()const
    {
        return basset_.active();
    }

    bool bassetEvaluated()const
    {
        return basset_.everEvaluated();
    }

    void commitBassetHistory()const { bassetModel().commitNow(); }

    void commitPairForces(const Plus::realx3ProcCMField& lift)const;

    void addBassetAddedMass(Plus::realProcCMField& addedMass)const;

    bool requireCellDistribution()const
    {
        return undisturbedSampler_.requireCellDistribution();
    }

    inline
    bool undisturbedIsMaterialDerivative()const
    {
        return undisturbedIsMaterialDerivative_;
    }

    virtual
    void calculateDragForceImplicit(
        const fluidAveraging& 	        fluidVelocity,
        const solidAveraging& 	        parVelocity,
        const Plus::realx3ProcCMField&  parRotVelocity,
        const Plus::realProcCMField& 	diameter,
        const distributionBase&         cellDistribution,
        Plus::realx3ProcCMField& 		particleForce,
        Plus::realProcCMField&          dragCoeff,
        Plus::realx3ProcCMField&        dragFluidVel) = 0;

    void calculateDragForce(
        const fluidAveraging& 	        fluidVelocity,
        const solidAveraging& 	        parVelocity,
        const Plus::realProcCMField& 	diameter,
        const distributionBase&         cellDistribution,
        Plus::realx3ProcCMField& 		particleForce) override;

    static
    uniquePtr<dragImplicit> create
    (
        const unresolvedCouplingSystem& uCS,
        const porosity& 				prsty
    );
};

}

#endif
