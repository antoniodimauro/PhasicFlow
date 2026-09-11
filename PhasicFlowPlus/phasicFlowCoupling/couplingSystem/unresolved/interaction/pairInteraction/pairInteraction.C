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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "pairInteraction.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "momentumSemiImplicitSphereUnresolvedCouplingSystem.hpp"

namespace pFlow::coupling
{

pairInteraction::pairInteraction
(
    const unresolvedCouplingSystem& uCS,
    const Foam::dictionary&         dict
)
:
    uCS_(uCS)
{
    active_ = dict.getOrDefault<Foam::Switch>("pairInteraction", Foam::Switch(false));

    if(!active_) return;

    cutWidths_ = dict.getOrDefault<Foam::scalar>("pairCutWidths", 5.0);
    pairVorticity_ = dict.getOrDefault<Foam::Switch>("pairVorticity", Foam::Switch(false));
    if(pairVorticity_)
    {
        Foam::Info
            << "    Pair vorticity: " << Green_Text("on")
            << ", dw_i = sum_j [W(sigma_seen) - W(sigma_true)](F_j), W the vorticity\n"
               "      of the Gaussian Stokeslet; added to the vorticity that sets the"
               " spin target, so a\n      neighbour turns the particle at half its"
               " Stokeslet vorticity (Jeffrey & Onishi 1984).\n";
    }

    Foam::Info
        << "    Pair interaction: " << Green_Text("near-field correction")
        << " du_i = sum_j [G(sigma_seen) - G(sigma_true)].F_j, G the Stokes flow\n"
           "      of a Gaussian force (Maxey & Patel 2001), sigma_true the force-coupling\n"
           "      pair width (Lomholt & Maxey 2003); neighbours within "
        << cutWidths_ << " widths.\n"
           "      Stokes flow: exact for Re_p << 1; the finite-Re near field (wake,\n"
           "      side-by-side repulsion) is not modelled.\n";
}

Foam::vector pairInteraction::gaussianFlow
(
    const Foam::vector& x,
    const Foam::vector& F,
    const Foam::scalar  sigma,
    const Foam::scalar  mu
)
{
    // Maxey & Patel (2001)
    Foam::scalar A, B;
    gaussianStokeslet(Foam::mag(x), sigma, mu, A, B);
    return A*F + B*(x & F)*x;
}

Foam::vector pairInteraction::gaussianVorticity
(
    const Foam::vector& x,
    const Foam::vector& F,
    const Foam::scalar  sigma,
    const Foam::scalar  mu
)
{
    const Foam::scalar pi = Foam::constant::mathematical::pi;
    const Foam::scalar r  = Foam::mag(x);
    if(r < 1e-3*sigma) return Foam::vector::zero;
    const Foam::scalar rho = r/sigma;
    // fraction of the Gaussian force inside r (the curl of a Gaussian
    // Stokeslet is that of a point force with this fraction of the force)
    const Foam::scalar M = std::erf(rho/Foam::sqrt(2.0))
                         - Foam::sqrt(2.0/pi)*rho*std::exp(-0.5*rho*rho);
    return (M/(4.0*pi*mu*r*r*r))*(F ^ x);
}

void pairInteraction::gaussianStokeslet
(
    Foam::scalar r,
    Foam::scalar sigma,
    Foam::scalar mu,
    Foam::scalar& A,
    Foam::scalar& B
)
{
    const Foam::scalar pi = Foam::constant::mathematical::pi;
    const Foam::scalar c  = 1.0/Foam::sqrt(2.0*pi);

    if(r < 0.05*sigma)
    {
        // Maxey & Patel (2001) eq. (39)
        A = c/(3.0*pi*mu*sigma);
        B = c/(30.0*pi*mu*sigma*sigma*sigma);
        return;
    }

    const Foam::scalar s2 = sigma*sigma/(r*r);
    const Foam::scalar e  = std::erf(r/(sigma*Foam::sqrt(2.0)));
    const Foam::scalar g  = std::exp(-0.5*r*r/(sigma*sigma));

    // Maxey & Patel (2001) eqs. (37)-(38)
    A = ((1.0 + s2)*e - 2.0*(sigma/r)*c*g)/(8.0*pi*mu*r);
    B = ((1.0 - 3.0*s2)*e + 6.0*(sigma/r)*c*g)/(8.0*pi*mu*r*r*r);
}

void pairInteraction::beginCoupling(const size_t nParLocal)
{
    if(!active_) return;

    fLoc_.assign  (nParLocal, realx3(0.0));
    sigLoc_.assign(nParLocal, 0.0);
    muLoc_.assign (nParLocal, 0.0);
    ownLoc_.assign(nParLocal, 0.0);
    duLoc_.assign (nParLocal, realx3(0.0));
    if(pairVorticity_) dwLoc_.assign(nParLocal, realx3(0.0));
}

void pairInteraction::publishWidth
(
    const size_t        i,
    const Foam::scalar  sigma,
    const Foam::scalar  mu
)
{
    if(!active_ || i >= sigLoc_.size()) return;

    sigLoc_[i] = static_cast<real>(sigma);
    muLoc_[i]  = static_cast<real>(mu);
    ownLoc_[i] = 1;
}

void pairInteraction::addForce(const size_t i, const Foam::vector& F)
{
    if(!active_ || i >= fLoc_.size()) return;

    fLoc_[i] += realx3(F.x(), F.y(), F.z());
}

void pairInteraction::evaluate(const bool faxen)
{
    if(!active_) return;

    auto& pMap = const_cast<particleMapping&>(uCS_.parMapping());
    auto& pds  = semiImplicitCoupling(uCS_).procDEM();

    const size_t nAll = pds.particlesCenterMassAllMaster().size();

    sigAll_.assign(nAll, 0.0);
    muAll_.assign (nAll, 0.0);
    ownAll_.assign(nAll, 0.0);

    {
        span<real> mineS(sigLoc_.data(), sigLoc_.size());
        span<real> allS (sigAll_.data(), sigAll_.size());
        span<real> mineM(muLoc_.data(),  muLoc_.size());
        span<real> allM (muAll_.data(),  muAll_.size());
        span<real> mineO(ownLoc_.data(), ownLoc_.size());
        span<real> allO (ownAll_.data(), ownAll_.size());

        if(!pMap.realScatteredComm().collectSum(mineS, allS) ||
           !pMap.realScatteredComm().collectSum(mineM, allM) ||
           !pMap.realScatteredComm().collectSum(mineO, allO))
        {
            fatalErrorInFunction<<"pairInteraction: collect to master failed"<<endl;
            Plus::processor::abort(0);
        }
    }

    duAll_.assign(nAll, realx3(0.0));
    if(pairVorticity_) dwAll_.assign(nAll, realx3(0.0));

    if(Plus::processor::isMaster())
    {
        auto ids = pds.particleIdAllMaster();
        auto pos = pds.particlesCenterMassAllMaster();
        auto dia = pds.particlesDiameterAllMaster();

        if(ids.size() < nAll || dia.size() < nAll)
        {
            fatalErrorInFunction
                << "pairInteraction: the DEM system does not publish a usable "
                   "particle id and diameter per particle." << endl;
            Plus::processor::abort(0);
        }

        fAll_.assign(nAll, realx3(0.0));
        for(size_t i=0; i<nAll; ++i)
        {
            auto f = fById_.find(ids[i]);
            if(f != fById_.end()) fAll_[i] = f->second;
        }

        Foam::scalar hMax = 0;
        for(size_t i=0; i<nAll; ++i)
        {
            if(ownAll_[i] > 0.5)
            {
                sigAll_[i] /= ownAll_[i];
                muAll_[i]  /= ownAll_[i];
            }
            const Foam::scalar a = 0.5*dia[i];
            const Foam::scalar sMax = Foam::max
            (
                Foam::sqrt(Foam::sqr(sigAll_[i]) + (faxen ? a*a/3.0 : 0.0)),
                Foam::sqrt(2.0/Foam::constant::mathematical::pi)*a
            );
            hMax = Foam::max(hMax, sMax*cutWidths_);
        }

        const Foam::scalar cell = hMax;

        Foam::scalar maxDu = 0;
        Foam::label  nPairs = 0;

        if(cell > 0)
        {
            auto key = [cell](const realx3& x, int dx, int dy, int dz)
            {
                const int64_t off = int64_t(1) << 20;
                const int64_t ix = static_cast<int64_t>(std::floor(x.x()/cell)) + dx + off;
                const int64_t iy = static_cast<int64_t>(std::floor(x.y()/cell)) + dy + off;
                const int64_t iz = static_cast<int64_t>(std::floor(x.z()/cell)) + dz + off;
                return (ix << 42) | (iy << 21) | iz;
            };

            std::unordered_map<int64_t, std::vector<size_t>> bins;
            bins.reserve(2*nAll);
            for(size_t j=0; j<nAll; ++j)
            {
                if(ownAll_[j] < 0.5) continue;
                bins[key(pos[j], 0, 0, 0)].push_back(j);
            }

            for(size_t i=0; i<nAll; ++i)
            {
                if(ownAll_[i] < 0.5) continue;

                const Foam::scalar ai  = 0.5*dia[i];
                const Foam::scalar mui = muAll_[i];
                if(mui <= 0) continue;

                const Foam::vector xi(pos[i].x(), pos[i].y(), pos[i].z());
                Foam::vector du = Foam::vector::zero;
                Foam::vector dw = Foam::vector::zero;

                for(int dx=-1; dx<=1; ++dx)
                for(int dy=-1; dy<=1; ++dy)
                for(int dz=-1; dz<=1; ++dz)
                {
                    auto b = bins.find(key(pos[i], dx, dy, dz));
                    if(b == bins.end()) continue;

                    for(const size_t j : b->second)
                    {
                        if(j == i) continue;

                        const Foam::vector Fj(fAll_[j].x(), fAll_[j].y(), fAll_[j].z());
                        if(Foam::magSqr(Fj) <= 0) continue;

                        const Foam::vector x = xi - Foam::vector(pos[j].x(), pos[j].y(), pos[j].z());
                        const Foam::scalar r = Foam::mag(x);
                        const Foam::scalar aj = 0.5*dia[j];

                        const Foam::scalar sSeen = Foam::sqrt
                        (
                            Foam::sqr(sigAll_[j]) + (faxen ? ai*ai/3.0 : 0.0)
                        );
                        const Foam::scalar sTrue = Foam::sqrt
                        (
                            (ai*ai + aj*aj)/Foam::constant::mathematical::pi
                        );

                        if(sigAll_[j] <= 0) continue;

                        const Foam::scalar sMax = Foam::max(sSeen, sTrue);
                        if(r > sMax*cutWidths_) continue;

                        // Lomholt & Maxey (2003)
                        du += gaussianFlow(x, Fj, sSeen, mui)
                            - gaussianFlow(x, Fj, sTrue, mui);
                        if(pairVorticity_)
                        {
                            // the vorticity is sampled at the centre: the kernel width alone
                            dw += gaussianVorticity(x, Fj, sigAll_[j], mui)
                                - gaussianVorticity(x, Fj, sTrue, mui);
                        }
                        ++nPairs;
                    }
                }

                duAll_[i] = realx3(du.x(), du.y(), du.z());
                if(pairVorticity_) dwAll_[i] = realx3(dw.x(), dw.y(), dw.z());
                maxDu = Foam::max(maxDu, Foam::mag(du));
            }
        }

        maxDuObs_  = maxDu;
        nPairsObs_ = nPairs;

        if(reportCount_++ % 200 == 0)
        {
            Foam::Info
                << Blue_Text("Pair interaction: ") << "|du|max = "
                << Yellow_Text(maxDu) << " m/s over " << nPairs
                << " neighbour pairs" << Foam::endl;
        }
    }

    {
        span<realx3> allD (duAll_.data(), duAll_.size());
        span<realx3> mineD(duLoc_.data(), duLoc_.size());

        if(!pMap.realx3ScatteredComm().distribute(allD, mineD))
        {
            fatalErrorInFunction<<"pairInteraction: distribute to ranks failed"<<endl;
            Plus::processor::abort(0);
        }
    }
    if(pairVorticity_)
    {
        span<realx3> allW (dwAll_.data(), dwAll_.size());
        span<realx3> mineW(dwLoc_.data(), dwLoc_.size());
        if(!pMap.realx3ScatteredComm().distribute(allW, mineW))
        {
            fatalErrorInFunction<<"pairInteraction: distribute vorticity to ranks failed"<<endl;
            Plus::processor::abort(0);
        }
    }
}

void pairInteraction::commit()
{
    if(!active_) return;

    auto& pMap = const_cast<particleMapping&>(uCS_.parMapping());
    auto& pds  = semiImplicitCoupling(uCS_).procDEM();

    const size_t nAll = pds.particlesCenterMassAllMaster().size();

    fAll_.assign(nAll, realx3(0.0));

    {
        span<realx3> mineF(fLoc_.data(), fLoc_.size());
        span<realx3> allF (fAll_.data(), fAll_.size());

        if(!pMap.realx3ScatteredComm().collectSum(mineF, allF))
        {
            fatalErrorInFunction<<"pairInteraction: collect forces failed"<<endl;
            Plus::processor::abort(0);
        }
    }

    if(Plus::processor::isMaster())
    {
        auto ids = pds.particleIdAllMaster();

        fById_.clear();
        fById_.reserve(2*nAll);
        for(size_t i=0; i<nAll && i<ids.size(); ++i)
        {
            fById_[ids[i]] = fAll_[i];
        }
    }
}

}
