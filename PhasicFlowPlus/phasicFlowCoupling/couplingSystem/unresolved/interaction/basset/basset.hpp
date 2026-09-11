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

// Maxey & Riley (1983) Phys. Fluids 26, 883; van Hinsberg, ten Thije Boonkkamp & Clercx (2011) J. Comput. Phys. 230, 1465

#ifndef __basset_hpp__
#define __basset_hpp__

#include <unordered_map>
#include <vector>

#include "OFCompatibleHeader.hpp"

#include "processorPlus.hpp"
#include "procCMFields.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class basset
{
public:

    static constexpr int nTail = 10;

    // Basset-Boussinesq; Mei & Adrian (1992) eq. (81)
    enum class historyKernel { Stokes, MeiAdrian };

    static constexpr int nFit = 201;

private:

    const unresolvedCouplingSystem& uCS_;

    bool            active_ = false;

    int             nWindow_ = 20;

    historyKernel   kernel_ = historyKernel::Stokes;

    // Dorgan & Loth (2007) eq. (5); Mei & Adrian (1992)
    Foam::scalar    c1_ = 2.0;
    Foam::scalar    c2_ = 0.105;

    std::vector<Foam::scalar> fitX_;
    std::vector<Foam::scalar> fitP_;

    mutable Foam::scalar minCrossover_ = GREAT;

    int             maxWindow_ = 0;

    Foam::scalar    twin_ = 0.0;

    mutable std::vector<Foam::scalar> wWin_;

    Foam::scalar    timeScale_ = 1.0;

    struct state
    {
        std::vector<Foam::vector> f;
        std::vector<Foam::scalar> dt;
        Foam::vector              F[nTail];
        Foam::vector              u;
    };

    std::unordered_map<uint32, state> hist_;

    bool                        hasPending_ = false;
    Foam::scalar                tPending_   = 0.0;
    Foam::scalar                dtPending_  = 0.0;
    std::vector<uint32>         idPend_;
    std::vector<realx3>         uEffPend_;
    std::vector<real>           ampPend_;
    std::vector<real>           cbPend_;
    std::vector<real>           dpPend_;
    std::vector<real>           nuPend_;

    Foam::scalar                tCommitted_ = 0.0;
    bool                        started_    = false;

    bool                        everEvaluated_ = false;

    mutable bool                reportedKernel_ = false;

    mutable Foam::scalar        reportedOverlap_ = 0.0;

    std::vector<realx3>         uEffLoc_;
    std::vector<real>           ampLoc_;
    std::vector<real>           cbLoc_;
    std::vector<real>           dpLoc_;
    std::vector<real>           nuLoc_;

    std::vector<realx3>         fExplLoc_;
    std::vector<real>           mbLoc_;

    std::vector<realx3>         uEffAll_;
    std::vector<real>           ampAll_;
    std::vector<real>           cbAll_;
    std::vector<real>           dpAll_;
    std::vector<real>           nuAll_;
    std::vector<realx3>         fExplAll_;
    std::vector<real>           mbAll_;

    // Mei & Adrian (1992) eq. (81)
    Foam::scalar setWindowWeights
    (
        const std::vector<Foam::scalar>& dt,
        Foam::scalar C = 0
    )const;

    // Dorgan & Loth (2007) eq. (5)
    Foam::scalar kernelC(Foam::scalar slip, Foam::scalar dp, Foam::scalar nu)const;

    Foam::scalar meiAdrianPsi(Foam::scalar s, Foam::scalar C)const;

    void buildTailFit();

    // van Hinsberg et al. (2011)
    void tailWeights(Foam::scalar C, Foam::scalar a[nTail])const;

    void reportOverlap(Foam::scalar overlap)const;

    int windowSpan(const std::vector<Foam::scalar>& dt)const;

    Foam::fileName statePath(const Foam::word& timeName)const;

    void writeState()const;

    void readState();

    void commitPending(
        const span<uint32>& ids,
        const span<realx3>& vel);

public:

    basset(const unresolvedCouplingSystem& uCS, const Foam::dictionary& dict);

    ~basset() = default;

    inline bool active()const { return active_; }

    inline bool everEvaluated()const { return everEvaluated_; }

    void reportKernelRegime(Foam::scalar slip, Foam::scalar dp,
                            Foam::scalar nu, Foam::scalar dt)const;

    void beginCoupling(size_t nParLocal);

    inline
    void publish(
        size_t                i,
        const Foam::vector&   uEff,
        Foam::scalar          amp,
        Foam::scalar          cB,
        Foam::scalar          dp,
        Foam::scalar          nu)
    {
        if(!active_) return;
        uEffLoc_[i] = realx3(uEff.x(), uEff.y(), uEff.z());
        ampLoc_[i]  = static_cast<real>(amp);
        cbLoc_[i]   = static_cast<real>(cB);
        dpLoc_[i]   = static_cast<real>(dp);
        nuLoc_[i]   = static_cast<real>(nu);
    }

    void evaluate(Foam::scalar tNow, Foam::scalar dtNow);

    void commitNow();

    inline
    Foam::vector explicitForce(size_t i)const
    {
        if(!active_) return Foam::vector::zero;
        return Foam::vector(fExplLoc_[i].x(), fExplLoc_[i].y(), fExplLoc_[i].z());
    }

    inline
    Foam::scalar addedMass(size_t i)const
    {
        return active_ ? static_cast<Foam::scalar>(mbLoc_[i]) : 0.0;
    }
};

}

#endif
