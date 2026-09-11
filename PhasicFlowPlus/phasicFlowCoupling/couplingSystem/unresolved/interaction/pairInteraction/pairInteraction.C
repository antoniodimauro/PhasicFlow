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
#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

#include "pairInteraction.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "momentumSemiImplicitSphereUnresolvedCouplingSystem.hpp"
#include "DEMSystemSemiImplicit.hpp"

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
        fOld_ = fAll_;

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

        blocks_.clear();

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

                        // the true separation at every gap: the force-coupling mobility is then a Gram matrix,
                        // positive definite for any configuration (solveStrengths); a separation held at a gap of
                        // a/4 left that form and made dense clusters indefinite (validation/pc-tests-2026-10-16)
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

                        // Lomholt & Maxey (2003): C_ij = G(x; sSeen) - G(x; sTrue), kept for the strength solve
                        Foam::scalar As, Bs, At, Bt;
                        gaussianStokeslet(r, sSeen, mui, As, Bs);
                        gaussianStokeslet(r, sTrue, mui, At, Bt);
                        blocks_.push_back({i, j, As - At, Bs - Bt, x});

                        if(Foam::magSqr(Fj) <= 0) continue;

                        du += (As - At)*Fj + (Bs - Bt)*(x & Fj)*x;
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

        // a new step (its start differs from the last one seen; a roll-back repeats the same start): the corrections
        // of the last pass become those under which the DEM produced the step-start velocities
        {
            const auto& tm = semiImplicitCoupling(uCS_).cMesh().mesh().time();
            const Foam::scalar t0 = tm.value() - tm.deltaTValue();
            if(Foam::mag(t0 - stepStart_) > 1e-9*Foam::max(tm.deltaTValue(), SMALL))
            {
                duPrevById_ = duLastById_;
                stepStart_  = t0;
            }
            duPrev_.assign(nAll, realx3(0.0));
            hasPrev_.assign(nAll, 0);
            duLastById_.clear();
            duLastById_.reserve(2*nAll);
            for(size_t i=0; i<nAll; ++i)
            {
                auto f = duPrevById_.find(ids[i]);
                hasPrev_[i] = (f != duPrevById_.end());
                duPrev_[i]  = hasPrev_[i] ? f->second : duAll_[i];
                duLastById_[ids[i]] = duAll_[i];
            }
        }

        if(reportCount_++ % 200 == 0)
        {
            Foam::Info
                << Blue_Text("Pair interaction: ") << "|du|max = "
                << Yellow_Text(maxDu) << " m/s over " << nPairs
                << " neighbour pairs";
            if(maxSweepsObs_ > 0)
            {
                Foam::Info
                    << "; strengths solved in at most " << maxSweepsObs_ << " sweeps (residual "
                    << maxResidualObs_ << ", largest Gershgorin radius " << maxRadiusObs_ << ")";
            }
            Foam::Info<< Foam::endl;
            maxSweepsObs_ = 0;
            maxRadiusObs_ = 0;
            maxResidualObs_ = 0;
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

// The strengths of one pass are, for every particle with neighbours, the Stokes drag on its slip at the end of the
// step: F_i = K0_i (w_i + alpha_i sum_j C_ij F_j + (1 - alpha_i) du_prev_i), K0_i = 3 pi mu d_i, w_i the undisturbed
// slip at the step start without the pair term, du_prev_i the correction under which the DEM produced the step-start
// velocity, and alpha_i the part of a change of the correction that the particle has not followed by the end of the
// step. Each rank returned F_raw_i = K0_i (w_i + du_used_i), so the system is
//     F_i/(alpha_i K0_i) - sum_j C_ij F_j = (F_raw_i - K0_i (du_used_i - (1 - alpha_i) du_prev_i))/(alpha_i K0_i),
// the backward (end-of-step) treatment of the pair interaction in time, with the fixed point of the explicit one.
// Used directly, F_raw is one Jacobi step within the pass and one explicit step in time; both diverge in a dense
// cluster (an eigenvalue of I - K0 C above 2; for overdamped particles a step-to-step gain of (1/K0 - C)^-1 C above
// 1). The matrix diag(1/(alpha K0)) - C is positive definite: 1/K0 = 1/(6 pi mu a) is the self-mobility of the
// force-coupling monopole of width a sqrt(2/pi) (Lomholt & Maxey 2003); with the true separations the pair part is
// the Stokes flow smoothed by sTrue minus that smoothed by the wider sSeen, a positive-weight kernel for equal spheres
// and positive up to O((a/sSeen)^2) for unequal ones; 1/alpha >= 1 only adds to the diagonal. Restarted GMRES on this
// mobility form converges when its symmetric part is positive definite; if it stops short, Jacobi sweeps with each
// update scaled by 1/(1 + R_i), R_i the Gershgorin radius alpha_i K0_i max_k sum_j sum_l |C_ij,kl|, halving the scaling
// when the residual grows over two sweeps, take over. The best iterate is used; if even that is off by more than the
// strengths themselves, F_raw is kept. A particle without a neighbour keeps F_raw (16 Oct 2026,
// validation/pc-tests-2026-10-16/README.md).
void pairInteraction::solveStrengths()
{
    auto& pds = semiImplicitCoupling(uCS_).procDEM();
    auto  dia = pds.particlesDiameterAllMaster();

    const size_t n = fAll_.size();
    if(duAll_.size() < n || muAll_.size() < n || dia.size() < n) return;

    // the unknowns: the particles with at least one neighbour
    std::vector<long>   slot(n, -1);
    std::vector<size_t> idx;
    for(const auto& bk : blocks_)
    {
        if(slot[bk.i] < 0)
        {
            slot[bk.i] = static_cast<long>(idx.size());
            idx.push_back(bk.i);
        }
    }
    const size_t m = idx.size();
    if(m == 0) return;

    std::vector<Foam::scalar> K(m), omega(m);
    std::vector<Foam::vector> rowSum(m, Foam::vector::zero);
    for(const auto& bk : blocks_)
    {
        for(int k=0; k<3; ++k)
        {
            Foam::scalar sum = 0;
            for(int l=0; l<3; ++l)
            {
                sum += Foam::mag((k == l ? bk.dA : 0.0) + bk.dB*bk.x[k]*bk.x[l]);
            }
            rowSum[slot[bk.i]][k] += sum;
        }
    }

    // the particle's response within the step: the DEM's implicit drag, v <- (v + A u_f)/(1 + A) per substep with
    // A = dt K/(m + m_added), leaves the fraction alpha = exp(-K Dt/(m + m_added)) of a change of the correction
    // unfollowed after a fluid step Dt (the DEM step is far below (m + m_added)/K). K, m_added: the values sent to
    // the DEM in the last exchange. A particle with no correction recorded for the previous step (inserted in this
    // step, or the first step after a restart) takes alpha = 1.
    auto* dem = semiImplicitCoupling(uCS_).couplingDEM();
    const Foam::scalar Dt = semiImplicitCoupling(uCS_).cMesh().mesh().time().deltaTValue();
    std::vector<Foam::scalar> K0(m), alpha(m, 1.0);
    {
        span<real> Kd, ma, mp;
        if(dem)
        {
            Kd = dem->parDragCoeff();
            ma = dem->parAddedMass();
            mp = dem->parMass();
        }
        for(size_t k=0; k<m; ++k)
        {
            const size_t i = idx[k];
            const bool known = i < hasPrev_.size() && hasPrev_[i];
            if(known && i < Kd.size() && i < ma.size() && i < mp.size() && mp[i] + ma[i] > 0 && Kd[i] > 0)
            {
                alpha[k] = Foam::max(std::exp(-Kd[i]*Dt/(mp[i] + ma[i])), 1e-6);
            }
        }
    }

    Foam::scalar Rmax = 0;
    for(size_t k=0; k<m; ++k)
    {
        const size_t i = idx[k];
        K0[k] = 3.0*Foam::constant::mathematical::pi*muAll_[i]*dia[i];
        K[k]  = alpha[k]*K0[k];
        const Foam::scalar R =
            K[k]*Foam::max(rowSum[k].x(), Foam::max(rowSum[k].y(), rowSum[k].z()));
        omega[k] = 1.0/(1.0 + R);
        Rmax = Foam::max(Rmax, R);
    }

    std::vector<Foam::vector> Fraw(n);
    Foam::scalar fScale = 0;
    for(size_t i=0; i<n; ++i)
    {
        Fraw[i] = Foam::vector(fAll_[i].x(), fAll_[i].y(), fAll_[i].z());
        fScale = Foam::max(fScale, Foam::mag(Fraw[i]));
    }
    if(fScale <= 0) return;

    const size_t N3 = 3*m;
    auto get = [](const std::vector<double>& v, size_t k)
    {
        return Foam::vector(v[3*k], v[3*k+1], v[3*k+2]);
    };
    auto add = [](std::vector<double>& v, size_t k, const Foam::vector& a)
    {
        v[3*k] += a.x(); v[3*k+1] += a.y(); v[3*k+2] += a.z();
    };

    // right-hand side: F = F_raw + K0 (alpha C F - du_used + (1 - alpha) du_prev), divided by K = alpha K0; with the
    // fixed strengths of neighbours that are not unknowns
    std::vector<double> rhs(N3, 0.0), x(N3, 0.0);
    for(size_t k=0; k<m; ++k)
    {
        const size_t i = idx[k];
        const Foam::vector du(duAll_[i].x(), duAll_[i].y(), duAll_[i].z());
        const Foam::vector dp = i < duPrev_.size()
            ? Foam::vector(duPrev_[i].x(), duPrev_[i].y(), duPrev_[i].z()) : du;
        add(rhs, k, (Fraw[i] - K0[k]*(du - (1.0 - alpha[k])*dp))/K[k]);
        // warm start: the strengths of the previous pass
        const Foam::vector x0 = i < fOld_.size()
            ? Foam::vector(fOld_[i].x(), fOld_[i].y(), fOld_[i].z()) : Fraw[i];
        add(x, k, x0);
    }
    for(const auto& bk : blocks_)
    {
        if(slot[bk.j] >= 0) continue;
        const Foam::vector& Fj = Fraw[bk.j];
        add(rhs, slot[bk.i], bk.dA*Fj + bk.dB*(bk.x & Fj)*bk.x);
    }

    // y = (1/K - C) v over the unknowns
    auto apply = [&](const std::vector<double>& v, std::vector<double>& y)
    {
        for(size_t k=0; k<m; ++k)
        {
            y[3*k] = v[3*k]/K[k]; y[3*k+1] = v[3*k+1]/K[k]; y[3*k+2] = v[3*k+2]/K[k];
        }
        for(const auto& bk : blocks_)
        {
            const long kj = slot[bk.j];
            if(kj < 0) continue;
            const Foam::vector Fj = get(v, static_cast<size_t>(kj));
            add(y, slot[bk.i], -(bk.dA*Fj + bk.dB*(bk.x & Fj)*bk.x));
        }
    };
    auto dot = [](const std::vector<double>& a, const std::vector<double>& b)
    {
        double s = 0;
        for(size_t q=0; q<a.size(); ++q) s += a[q]*b[q];
        return s;
    };

    std::vector<double> Ax(N3), r(N3), w(N3);

    // the residual of the force equation, max_i K_i |r_i| / max|F_raw|
    auto forceResidual = [&]()
    {
        apply(x, Ax);
        Foam::scalar res = 0;
        for(size_t k=0; k<m; ++k)
        {
            const Foam::vector rk(rhs[3*k] - Ax[3*k], rhs[3*k+1] - Ax[3*k+1], rhs[3*k+2] - Ax[3*k+2]);
            res = Foam::max(res, K[k]*Foam::mag(rk));
        }
        return res/fScale;
    };

    const Foam::scalar tol = 1e-6;
    Foam::label iterations = 0;
    Foam::scalar res = forceResidual();

    // the best iterate so far (smallest force residual): the one used if the solve stops short
    std::vector<double> xBest = x;
    Foam::scalar resBest = res;
    auto keepBest = [&]()
    {
        if(res < resBest)
        {
            resBest = res;
            xBest = x;
        }
    };

    // restarted GMRES(40)
    if(res >= tol)
    {
        const int mR = 40;
        const int maxRestarts = 50;
        std::vector<std::vector<double>> V(mR + 1, std::vector<double>(N3, 0.0));
        std::vector<double> H((mR + 1)*mR, 0.0), cs(mR, 0.0), sn(mR, 0.0), g(mR + 1, 0.0), yv(mR, 0.0);

        for(int restart=0; restart<maxRestarts && res >= tol; ++restart)
        {
            apply(x, Ax);
            for(size_t q=0; q<N3; ++q) r[q] = rhs[q] - Ax[q];
            const double beta = std::sqrt(dot(r, r));
            if(beta <= 0) break;
            for(size_t q=0; q<N3; ++q) V[0][q] = r[q]/beta;
            std::fill(g.begin(), g.end(), 0.0);
            g[0] = beta;

            int kEnd = 0;
            for(int k=0; k<mR; ++k)
            {
                apply(V[k], w);
                ++iterations;
                for(int j=0; j<=k; ++j)
                {
                    const double h = dot(w, V[j]);
                    H[j*mR + k] = h;
                    for(size_t q=0; q<N3; ++q) w[q] -= h*V[j][q];
                }
                const double hk = std::sqrt(dot(w, w));
                H[(k + 1)*mR + k] = hk;
                if(hk > 0) for(size_t q=0; q<N3; ++q) V[k + 1][q] = w[q]/hk;

                for(int j=0; j<k; ++j)
                {
                    const double t = cs[j]*H[j*mR + k] + sn[j]*H[(j + 1)*mR + k];
                    H[(j + 1)*mR + k] = -sn[j]*H[j*mR + k] + cs[j]*H[(j + 1)*mR + k];
                    H[j*mR + k] = t;
                }
                const double ha = H[k*mR + k], hb = H[(k + 1)*mR + k];
                const double den = std::sqrt(ha*ha + hb*hb);
                cs[k] = den > 0 ? ha/den : 1.0;
                sn[k] = den > 0 ? hb/den : 0.0;
                H[k*mR + k] = cs[k]*ha + sn[k]*hb;
                H[(k + 1)*mR + k] = 0.0;
                g[k + 1] = -sn[k]*g[k];
                g[k]     =  cs[k]*g[k];
                kEnd = k + 1;
                if(std::abs(g[k + 1]) <= 1e-12*beta || hk <= 0) break;
            }

            for(int j=kEnd-1; j>=0; --j)
            {
                double sum = g[j];
                for(int l=j+1; l<kEnd; ++l) sum -= H[j*mR + l]*yv[l];
                yv[j] = H[j*mR + j] != 0 ? sum/H[j*mR + j] : 0.0;
            }
            for(int j=0; j<kEnd; ++j)
            {
                for(size_t q=0; q<N3; ++q) x[q] += yv[j]*V[j][q];
            }
            res = forceResidual();
            keepBest();
        }
    }

    // fallback: damped Jacobi sweeps from the GMRES iterate
    Foam::scalar scale = 1.0;
    bool fallback = false;
    if(res >= tol)
    {
        fallback = true;
        const Foam::label maxSweeps = 5000;
        x = xBest;
        res = resBest;
        Foam::scalar resPrev = res, resPrev2 = res;
        for(Foam::label sweep=1; sweep<=maxSweeps && res >= tol; ++sweep)
        {
            apply(x, Ax);
            for(size_t k=0; k<m; ++k)
            {
                for(int c=0; c<3; ++c)
                {
                    x[3*k + c] += scale*omega[k]*K[k]*(rhs[3*k + c] - Ax[3*k + c]);
                }
            }
            ++iterations;
            res = forceResidual();
            keepBest();
            // growth over two consecutive sweeps halves the scaling
            if(sweep > 2 && res > resPrev && resPrev > resPrev2 && scale > 1.0/64.0) scale *= 0.5;
            resPrev2 = resPrev;
            resPrev = res;
        }
    }

    x = xBest;
    res = resBest;
    // even the best iterate is off by more than the strengths themselves: keep the strengths of the ranks (F_raw)
    const bool useRaw = res > 1.0;

    maxSweepsObs_   = Foam::max(maxSweepsObs_, iterations);
    maxRadiusObs_   = Foam::max(maxRadiusObs_, Rmax);
    maxResidualObs_ = Foam::max(maxResidualObs_, res);

    if(iterations > 200 || fallback)
    {
        Foam::Info<< "Pair strengths: " << iterations << " iterations" << (fallback ? " (Jacobi fallback)" : "")
            << ", residual " << res << ", largest Gershgorin radius " << Rmax << ", scaling " << scale << Foam::endl;
    }

    if(res >= tol)
    {
        // the particles with the largest remaining residual, by id
        auto ids = pds.particleIdAllMaster();
        apply(x, Ax);
        std::vector<std::pair<Foam::scalar, size_t>> order;
        for(size_t k=0; k<m; ++k)
        {
            const Foam::vector rk(rhs[3*k] - Ax[3*k], rhs[3*k+1] - Ax[3*k+1], rhs[3*k+2] - Ax[3*k+2]);
            order.push_back({K[k]*Foam::mag(rk), idx[k]});
        }
        std::sort(order.begin(), order.end(), [](const auto& p, const auto& q){ return p.first > q.first; });
        std::string list;
        for(size_t k=0; k<order.size() && k<5; ++k)
        {
            const size_t i = order[k].second;
            list += " " + std::to_string(i < ids.size() ? ids[i] : static_cast<uint32>(i));
        }
        WARNING<<"pairInteraction: the strengths of this pass stopped at residual "<<res<<" after "<<iterations
            <<" iterations (largest Gershgorin radius "<<Rmax<<"), ids"<<list
            <<(useRaw ? "; the strengths of the ranks (F_raw) are kept" : "; the best iterate is used")
            <<END_WARNING;
    }

    if(useRaw) return;

    for(size_t k=0; k<m; ++k)
    {
        const Foam::vector Fk = get(x, k);
        fAll_[idx[k]] = realx3(Fk.x(), Fk.y(), Fk.z());
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
        // no pair within the cutoff: no block, nothing to solve, and every strength stays F_raw as the ranks sent it
        if(!blocks_.empty()) solveStrengths();

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
