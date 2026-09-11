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

// Feuillebois & Lasek (1978) Q. J. Mech. Appl. Math. 31, 435; Zhang & Stone (1998) J. Fluid Mech. 367, 329

#ifndef __rotationalHistory_hpp__
#define __rotationalHistory_hpp__

#include <unordered_map>
#include <vector>

#include "OFCompatibleHeader.hpp"

#include "processorPlus.hpp"
#include "procCMFields.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class rotationalHistory
{
public:

    static constexpr int nModes = 32;

private:

    const unresolvedCouplingSystem& uCS_;

    bool            active_ = false;

    int             nWindow_ = 20;

    Foam::scalar    lamQ_[nModes];
    Foam::scalar    ampQ_[nModes];

    mutable std::vector<Foam::scalar> wWin_;

    struct state
    {
        std::vector<Foam::vector> w;
        std::vector<Foam::scalar> dt;
        Foam::vector              F[nModes];
        Foam::vector              wAmb;
    };

    std::unordered_map<uint32, state> hist_;

    bool                        hasPending_ = false;
    Foam::scalar                tPending_   = 0.0;
    Foam::scalar                dtPending_  = 0.0;
    std::vector<uint32>         idPend_;
    std::vector<realx3>         wAmbPend_;
    std::vector<real>           tnuPend_;
    std::vector<real>           crPend_;
    std::vector<real>           kwPend_;

    Foam::scalar                tCommitted_ = 0.0;
    bool                        started_    = false;

    bool                        everEvaluated_ = false;
    mutable bool                reportedQuad_  = false;

    std::vector<realx3>         wAmbLoc_;
    std::vector<real>           tnuLoc_;
    std::vector<real>           crLoc_;
    std::vector<real>           kwLoc_;

    std::vector<realx3>         tExplLoc_;
    std::vector<real>           iaLoc_;

    std::vector<realx3>         wAmbAll_;
    std::vector<real>           tnuAll_;
    std::vector<real>           crAll_;
    std::vector<real>           kwAll_;
    std::vector<realx3>         tExplAll_;
    std::vector<real>           iaAll_;

    void setQuadrature();

    Foam::scalar setWindowWeights(
        const std::vector<Foam::scalar>& dt,
        Foam::scalar                     tnu)const;

    Foam::fileName statePath(const Foam::word& timeName)const;
    void writeState()const;
    void readState();

    void commitPending(
        const span<uint32>&  ids,
        const span<realx3>&  rVel);

public:

    rotationalHistory(const unresolvedCouplingSystem& uCS,
                      const Foam::dictionary& dict);

    ~rotationalHistory() = default;

    inline bool active()const { return active_; }

    inline bool everEvaluated()const { return everEvaluated_; }

    void checkQuadrature(Foam::scalar xw)const;

    void beginCoupling(size_t nParLocal);

    inline
    void publish(
        size_t              i,
        const Foam::vector& wAmb,
        Foam::scalar        tnu,
        Foam::scalar        cR,
        Foam::scalar        kOmega)
    {
        if(!active_) return;
        wAmbLoc_[i] = realx3(wAmb.x(), wAmb.y(), wAmb.z());
        tnuLoc_[i]  = static_cast<real>(tnu);
        crLoc_[i]   = static_cast<real>(cR);
        kwLoc_[i]   = static_cast<real>(kOmega);
    }

    void evaluate(Foam::scalar tNow, Foam::scalar dtNow);

    void commitNow();

    inline Foam::vector torque(size_t i)const
    {
        if(!active_ || i >= tExplLoc_.size()) return Foam::vector::zero;
        const realx3& t = tExplLoc_[i];
        return Foam::vector(t.x(), t.y(), t.z());
    }

    inline Foam::scalar addedInertia(size_t i)const
    {
        if(!active_ || i >= iaLoc_.size()) return 0.0;
        return iaLoc_[i];
    }
};

}

#endif
