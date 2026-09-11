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
// Maxey & Patel (2001) Int. J. Multiphase Flow 27, 1603; Lomholt & Maxey (2003) J. Comput. Phys. 184, 381

#ifndef __pairInteraction_hpp__
#define __pairInteraction_hpp__

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "OFCompatibleHeader.hpp"
#include "processorPlus.hpp"
#include "procCMFields.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

// Maxey & Patel (2001) eqs. (33), (37), (38); Lomholt & Maxey (2003)
class pairInteraction
{
    const unresolvedCouplingSystem& uCS_;

    bool            active_ = false;

    Foam::scalar    cutWidths_ = 5.0;

    std::vector<realx3> fLoc_;
    std::vector<real>   sigLoc_;
    std::vector<real>   muLoc_;
    std::vector<real>   ownLoc_;

    std::vector<realx3> duLoc_;

    // pairVorticity: the vorticity of each neighbour's Stokeslet that the
    // kernel smears out, added to the vorticity that sets the spin target
    bool                pairVorticity_ = false;
    std::vector<realx3> dwLoc_;
    std::vector<realx3> dwAll_;

    std::unordered_map<uint32, realx3> fById_;

    std::vector<realx3> fAll_;
    std::vector<realx3> duAll_;
    std::vector<real>   sigAll_;
    std::vector<real>   muAll_;
    std::vector<real>   ownAll_;

    mutable Foam::scalar maxDuObs_ = 0;
    mutable Foam::label  nPairsObs_ = 0;
    mutable Foam::label  reportCount_ = 0;

public:

    pairInteraction(const unresolvedCouplingSystem& uCS, const Foam::dictionary& dict);

    bool active()const { return active_; }

    // Maxey & Patel (2001) eqs. (37)-(38)
    static void gaussianStokeslet
    (
        Foam::scalar r,
        Foam::scalar sigma,
        Foam::scalar mu,
        Foam::scalar& A,
        Foam::scalar& B
    );

    void beginCoupling(size_t nParLocal);

    static Foam::vector gaussianFlow
    (
        const Foam::vector& x,
        const Foam::vector& F,
        Foam::scalar sigma,
        Foam::scalar mu
    );

    void publishWidth(size_t i, Foam::scalar sigma, Foam::scalar mu);

    void evaluate(bool faxen);

    Foam::vector du(size_t i)const
    {
        if(!active_ || i >= duLoc_.size()) return Foam::vector::zero;
        const realx3& d = duLoc_[i];
        return Foam::vector(d.x(), d.y(), d.z());
    }

    // vorticity correction at particle i (pairVorticity yes)
    Foam::vector dVort(size_t i)const
    {
        if(!active_ || !pairVorticity_ || i >= dwLoc_.size()) return Foam::vector::zero;
        const realx3& d = dwLoc_[i];
        return Foam::vector(d.x(), d.y(), d.z());
    }

    bool pairVorticity()const { return active_ && pairVorticity_; }

    // vorticity of the Gaussian Stokeslet of width sigma at x:
    // (F x x)/(4 pi mu r^3) [erf(r/(sigma sqrt2)) - sqrt(2/pi) (r/sigma) exp(-r^2/(2 sigma^2))]
    static Foam::vector gaussianVorticity
    (
        const Foam::vector& x,
        const Foam::vector& F,
        Foam::scalar sigma,
        Foam::scalar mu
    );

    void addForce(size_t i, const Foam::vector& F);

    void commit();
};

}

#endif
