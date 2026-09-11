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

// Feuillebois & Lasek (1978) Q. J. Mech. Appl. Math. 31, 435; Goldman et al. (1967) Chem. Eng. Sci. 22, 637 and 653

#ifndef __surfaceRotationTorqueImplicit_hpp__
#define __surfaceRotationTorqueImplicit_hpp__

#include "surfaceRotationTorque.hpp"
#include "fluidFieldSampler.hpp"
#include "rotationalHistory.hpp"
#include "wallCorrection.hpp"

namespace pFlow::coupling
{

class pairInteraction;

class surfaceRotationTorqueImplicit
:
    public surfaceRotationTorque
{
private:

    const pairInteraction* pair_ = nullptr;

    fluidFieldSampler    torqueUSampler_;
    fluidFieldSampler    torqueCurlUSampler_;

    rotationalHistory    history_;

    bool                 wallCorr_      = false;
    Foam::scalar         minWallRadii_  = 2.0;

protected:


    rotationalHistory& history() { return history_; }

    const pairInteraction* pair()const { return pair_; }

    const Plus::procCMField<Foam::vector>&
    sampleFluidVel(const Foam::volVectorField& U)
    {
        return torqueUSampler_.sample(U);
    }

    const Plus::procCMField<Foam::vector>&
    sampleVorticity(const Foam::volVectorField& curlU)
    {
        return torqueCurlUSampler_.sample(curlU);
    }

    static Foam::vector wallCouple
    (
        const Foam::scalar  hOverA,
        const Foam::vector& nw,
        const Foam::vector& wp,
        const Foam::vector& curlU,
        const Foam::vector& up,
        const Foam::scalar  kOmega,
        const Foam::scalar  shearC,
        const Foam::scalar  dp
    )
    {
        const Foam::scalar tR = wallRotTorque(hOverA);
        const Foam::scalar tS = wallShearTorque(hOverA);
        const Foam::scalar tT = wallTransTorque(hOverA);

        const Foam::vector wPar = wp    - (wp    & nw)*nw;
        const Foam::vector cPar = curlU - (curlU & nw)*nw;
        const Foam::vector vPar = up    - (up    & nw)*nw;

        return -kOmega*(tR - 1.0)*wPar
             +  shearC*(tS - 1.0)*cPar
             +  2.0*(kOmega/dp)*tT*(nw ^ vPar);
    }


public:

    void setPair(const pairInteraction* p) { pair_ = p; }

    bool wallCorrection()const { return wallCorr_; }

    Foam::scalar minWallRadii()const { return minWallRadii_; }

    bool requireCellDistribution()const
    {
        return torqueUSampler_.requireCellDistribution()
            || torqueCurlUSampler_.requireCellDistribution();
    }

    TypeInfo("surfaceRotationTorqueImplicit");

    surfaceRotationTorqueImplicit(
        const unresolvedCouplingSystem& uCS,
        const porosity&                 prsty);

    virtual ~surfaceRotationTorqueImplicit() = default;

    virtual
    void calculateSurfaceTorqueImplicit(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque,
        Plus::realProcCMField&          rotDragCoeff) = 0;

    void calculateSurfaceTorque(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque) override;

    void calculateTorque(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque,
        Plus::realProcCMField&          rotDragCoeff,
        Plus::realProcCMField&          addedInertia);

    void commitHistory() { history_.commitNow(); }

    bool historyActive()const { return history_.active(); }

    bool historyEvaluated()const { return history_.everEvaluated(); }

    static
    uniquePtr<surfaceRotationTorqueImplicit> create
    (
        const unresolvedCouplingSystem& uCS,
        const porosity&                 prsty
    );
};

}

#endif
