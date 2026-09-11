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

#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include <cmath>

#include <fstream>
#include <iomanip>

#include "basset.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "momentumSemiImplicitSphereUnresolvedCouplingSystem.hpp"

namespace pFlow::coupling
{

static constexpr int     nvH = 10;
static constexpr double  tvH[nvH] =
    { 0.1, 0.3, 1.0, 3.0, 10.0, 40.0, 190.0, 1000.0, 6500.0, 50000.0 };
static constexpr double  avH[nvH] =
{
    0.23477481312586, 0.28549576238194, 0.28479416718255, 0.26149775537574,
    0.32056200511938, 0.35354490689146, 0.39635904496921, 0.42253908596514,
    0.48317384225265, 0.63661146557001
};

basset::basset
(
    const unresolvedCouplingSystem& uCS,
    const Foam::dictionary&         dict
)
:
    uCS_(uCS)
{
    active_ = dict.getOrDefault<Foam::Switch>("bassetHistory", Foam::Switch(false));

    if(!active_) return;

    nWindow_ = dict.getOrDefault<Foam::label>("bassetWindowSteps", 20);

    if(nWindow_ < 2)
    {
        fatalErrorInFunction
            << "bassetWindowSteps must be at least 2, got " << nWindow_
            << " in " << dict.name() << endl;
        Plus::processor::abort(0);
    }

    maxWindow_ = 8*nWindow_;

    const Foam::word kName =
        dict.getOrDefault<Foam::word>("bassetKernel", Foam::word("Stokes"));

    if(kName == "Stokes")
    {
        kernel_ = historyKernel::Stokes;
    }
    else if(kName == "MeiAdrian")
    {
        kernel_ = historyKernel::MeiAdrian;
        c1_ = dict.getOrDefault<Foam::scalar>("bassetKernelC1", 2.0);
        c2_ = dict.getOrDefault<Foam::scalar>("bassetKernelC2", 0.105);
        if(c1_ <= 0 || c2_ < 0)
        {
            fatalErrorInFunction
                << "bassetKernelC1 must be > 0 and bassetKernelC2 >= 0, got "
                << c1_ << " and " << c2_ << " in " << dict.name() << endl;
            Plus::processor::abort(0);
        }
        buildTailFit();
    }
    else
    {
        fatalErrorInFunction
            << "bassetKernel must be Stokes or MeiAdrian, got " << kName
            << " in " << dict.name() << endl;
        Plus::processor::abort(0);
    }

    readState();

    Foam::Info
        << "    Basset-Boussinesq history force is " << Green_Text("active")
        << ": van Hinsberg et al., JCP 230 (2011) 1465.\n"
           "      Exact kernel inside t_win = " << nWindow_
        << " initial steps, " << nTail << " exponentials beyond it.\n"
           "      t_win is a fixed time and the window is sub-divided on each "
           "particle's\n      own step lags, so adjustTimeStep moves neither "
           "the split point nor\n      the tail time constants.\n";

    if(kernel_ == historyKernel::MeiAdrian)
    {
        Foam::Info
            << "      Kernel: " << Green_Text("finite-Re") << ", Dorgan & Loth (2007) eq (5) "
               "with c1 = " << c1_ << ", c2 = " << c2_ << "\n"
               "      (Mei & Adrian 1992: 2, 0.105): the Stokes kernel times\n"
               "      psi(s) = (1 + (C s^(3/2))^(1/c1))^-c1, C from the instantaneous "
               "slip,\n      f_H = (0.75 + c2 Re)^3. In the "
               "window\n      each interval is weighted by psi at its mid-lag; the "
               "tail weights are\n      refitted to psi/sqrt(s) at every step "
               "(least squares over 1 to 1e4 window\n      lengths, the same "
               << nTail << " time constants).\n\n";
    }
    else
    {
        Foam::Info<< "      Kernel: Stokes (Basset-Boussinesq).\n\n";
    }
}

Foam::scalar basset::kernelC
(
    Foam::scalar slip,
    Foam::scalar dp,
    Foam::scalar nu
)const
{
    if(kernel_ != historyKernel::MeiAdrian) return 0.0;
    if(slip <= 0 || dp <= 0 || nu <= 0) return 0.0;

    const Foam::scalar pi = Foam::constant::mathematical::pi;
    const Foam::scalar Re = slip*dp/nu;
    const Foam::scalar fH = Foam::pow(0.75 + c2_*Re, 3);

    // Dorgan & Loth (2007) eq. (5); Mei & Adrian (1992) eq. (81)
    return Foam::sqrt(pi)*Foam::pow(slip, 3)/(2.0*Foam::pow(nu, 1.5)*fH);
}

Foam::scalar basset::meiAdrianPsi(Foam::scalar s, Foam::scalar C)const
{
    if(C <= 0 || s <= 0) return 1.0;
    return Foam::pow(1.0 + Foam::pow(C*Foam::pow(s, 1.5), 1.0/c1_), -c1_);
}

void basset::buildTailFit()
{
    // van Hinsberg et al. (2011) tail exponentials, least-squares fit
    fitX_.resize(nFit);
    std::vector<Foam::scalar> A(nFit*nTail);

    for(int j=0; j<nFit; ++j)
    {
        const Foam::scalar x = Foam::pow(10.0, 4.0*j/(nFit - 1));
        fitX_[j] = x;
        for(int k=0; k<nTail; ++k)
        {
            A[j*nTail + k] = Foam::sqrt(Foam::constant::mathematical::e/tvH[k])
                           * Foam::exp(-x/(2.0*tvH[k]))*Foam::sqrt(x);
        }
    }

    Foam::scalar M[nTail][nTail];
    for(int r=0; r<nTail; ++r)
    {
        for(int c=0; c<nTail; ++c)
        {
            Foam::scalar sum = 0;
            for(int j=0; j<nFit; ++j) sum += A[j*nTail + r]*A[j*nTail + c];
            M[r][c] = sum;
        }
    }

    std::vector<Foam::scalar> B(nTail*nFit);
    for(int k=0; k<nTail; ++k)
        for(int j=0; j<nFit; ++j)
            B[k*nFit + j] = A[j*nTail + k];

    // Gauss-Jordan, partial pivoting
    for(int c=0; c<nTail; ++c)
    {
        int piv = c;
        for(int r=c+1; r<nTail; ++r)
            if(Foam::mag(M[r][c]) > Foam::mag(M[piv][c])) piv = r;
        if(piv != c)
        {
            for(int q=0; q<nTail; ++q) std::swap(M[c][q], M[piv][q]);
            for(int j=0; j<nFit; ++j) std::swap(B[c*nFit + j], B[piv*nFit + j]);
        }
        const Foam::scalar d = M[c][c];
        for(int q=0; q<nTail; ++q) M[c][q] /= d;
        for(int j=0; j<nFit; ++j) B[c*nFit + j] /= d;
        for(int r=0; r<nTail; ++r)
        {
            if(r == c) continue;
            const Foam::scalar f = M[r][c];
            if(f == 0) continue;
            for(int q=0; q<nTail; ++q) M[r][q] -= f*M[c][q];
            for(int j=0; j<nFit; ++j) B[r*nFit + j] -= f*B[c*nFit + j];
        }
    }

    fitP_ = B;
}

void basset::tailWeights(Foam::scalar C, Foam::scalar a[nTail])const
{
    if(kernel_ != historyKernel::MeiAdrian)
    {
        for(int k=0; k<nTail; ++k) a[k] = avH[k];
        return;
    }

    for(int k=0; k<nTail; ++k) a[k] = 0;

    for(int j=0; j<nFit; ++j)
    {
        const Foam::scalar psi = meiAdrianPsi(fitX_[j]*twin_, C);
        for(int k=0; k<nTail; ++k) a[k] += fitP_[k*nFit + j]*psi;
    }
}

Foam::fileName basset::statePath(const Foam::word& timeName)const
{
    const Foam::Time& rt = uCS_.cMesh().mesh().time();

    return rt.globalPath()/timeName/"uniform"/"bassetHistory";
}

void basset::writeState()const
{
    if(!active_ || !Plus::processor::isMaster()) return;

    const Foam::Time& rt = uCS_.cMesh().mesh().time();
    const Foam::fileName f = statePath(rt.timeName());

    Foam::mkDir(f.path());

    std::ofstream out(f.c_str());
    if(!out) { Foam::Warning<<"basset: cannot write "<<f<<Foam::endl; return; }

    out << std::setprecision(17);
    out << "nTail "     << nTail      << "\n";
    out << "nWindow "   << nWindow_   << "\n";
    out << "tWin "      << twin_      << "\n";
    out << "tCommitted "<< tCommitted_<< "\n";
    out << "nParticles "<< hist_.size() << "\n";

    for(const auto& kv : hist_)
    {
        out << kv.first;

        for(int k=0; k<nTail; ++k)
        {
            out << ' ' << kv.second.F[k].x() << ' ' << kv.second.F[k].y()
                << ' ' << kv.second.F[k].z();
        }

        out << ' ' << kv.second.dt.size();

        for(const Foam::vector& fv : kv.second.f)
        {
            out << ' ' << fv.x() << ' ' << fv.y() << ' ' << fv.z();
        }

        for(const Foam::scalar sv : kv.second.dt) out << ' ' << sv;

        out << ' ' << kv.second.u.x() << ' ' << kv.second.u.y()
            << ' ' << kv.second.u.z();
        out << "\n";
    }
}

void basset::readState()
{
    if(!active_ || !Plus::processor::isMaster()) return;

    const Foam::Time& rt = uCS_.cMesh().mesh().time();
    const Foam::fileName f = statePath(rt.timeName());

    std::ifstream in(f.c_str());
    if(!in) return;

    std::string key;
    int    nt = 0, nw = 0;
    double tw = 0, tc = 0;
    long   np = 0;

    in >> key >> nt >> key >> nw >> key >> tw >> key >> tc >> key >> np;

    if(nt != nTail || nw != nWindow_)
    {
        Foam::Warning
            << "basset: " << f << " was written with nTail = " << nt
            << ", nWindow = " << nw << " but this run uses " << nTail
            << ", " << nWindow_ << ". Ignoring it and starting the history "
            << "from zero." << Foam::endl;
        return;
    }

    hist_.clear();
    hist_.reserve(static_cast<size_t>(np)*2);

    for(long n=0; n<np; ++n)
    {
        uint32 id = 0;
        if(!(in >> id)) break;

        state st;
        double x,y,z;

        for(int k=0; k<nTail; ++k)
        {
            in >> x >> y >> z; st.F[k] = Foam::vector(x,y,z);
        }

        size_t ndt = 0;
        in >> ndt;

        if(ndt == 0 || ndt > size_t(maxWindow_) + 1) break;

        st.f.resize(ndt + 1);
        for(size_t m=0; m<ndt+1; ++m)
        {
            in >> x >> y >> z; st.f[m] = Foam::vector(x,y,z);
        }

        st.dt.resize(ndt);
        for(size_t m=0; m<ndt; ++m) in >> st.dt[m];

        in >> x >> y >> z; st.u = Foam::vector(x,y,z);

        hist_.emplace(id, st);
    }

    twin_       = tw;
    tCommitted_ = tc;
    started_    = true;

    Foam::Info
        << "    Basset history resumed from " << f.name() << ": "
        << hist_.size() << " particles, committed at t = " << tCommitted_
        << " s.\n";
}

void basset::reportKernelRegime
(
    Foam::scalar slip,
    Foam::scalar dp,
    Foam::scalar nu,
    Foam::scalar dt
)const
{
    if(reportedKernel_ || slip <= 0 || dp <= 0 || nu <= 0) return;
    reportedKernel_ = true;

    if(kernel_ == historyKernel::MeiAdrian)
    {
        Foam::Info
            << Blue_Text("Basset kernel: ")
            << "finite-Re (c1 = " << c1_ << ", c2 = " << c2_ << "), in use; the "
               "per-step crossover time is reported with the history force.\n";
        return;
    }

    const Foam::scalar a   = 0.5*dp;
    const Foam::scalar ReP = slip*dp/nu;
    const Foam::scalar fH  = 0.75 + 0.105*ReP;
    const Foam::scalar sc  = Foam::pow(4.0/Foam::constant::mathematical::pi, 1.0/3.0)
                           * nu*fH*fH/(slip*slip);

    const Foam::scalar s    = twin_ > 0 ? twin_ : nWindow_*dt;
    const Foam::scalar t1   = Foam::pow(Foam::constant::mathematical::pi*s*nu/(a*a), 0.25);
    const Foam::scalar t2   = Foam::sqrt(0.5*Foam::constant::mathematical::pi
                                        *Foam::pow(slip,3)*s*s/(a*nu*Foam::pow(fH,3)));
    const Foam::scalar ratio = (a/Foam::sqrt(Foam::constant::mathematical::pi*nu*s))
                             > 0
                             ? Foam::pow(t1+t2,-2)
                               / (a/Foam::sqrt(Foam::constant::mathematical::pi*nu*s))
                             : 1.0;

    Foam::Info
        << Blue_Text("Basset kernel: ")
        << "Stokes. At Re_p = " << Yellow_Text(ReP)
        << " the finite-Re kernel of Mei & Adrian (1992) eq (81)\n"
           "      is " << Yellow_Text(ratio)
        << " of it at the window edge (" << s << " s); crossover at "
        << Yellow_Text(sc) << " s.\n";

    if(ratio < 0.8)
    {
        if(Plus::processor::isMaster())
        {
            WARNING
                << "The Stokes history kernel over-predicts the memory by about "
                << 1.0/Foam::max(ratio, Foam::SMALL) << "x at the window edge for "
                << "this slip. It is accurate at small slip, which is the trapping "
                << "regime, but a particle in transit at this Re_p carries too much "
                << "history force." << END_WARNING;
        }
    }
}

int basset::windowSpan(const std::vector<Foam::scalar>& dt)const
{
    const Foam::scalar reach = twin_*(1.0 - 1.0e-12);

    Foam::scalar x = 0.0;

    for(size_t n=0; n<dt.size(); ++n)
    {
        x += dt[n];

        if(x >= reach) return static_cast<int>(n) + 1;
    }

    return static_cast<int>(dt.size());
}

void basset::reportOverlap(Foam::scalar overlap)const
{
    if(overlap <= reportedOverlap_ + 0.05) return;

    reportedOverlap_ = overlap;

    Foam::Info
        << Blue_Text("Basset window: ")
        << "edge overruns t_win by " << Yellow_Text(100.0*overlap)
        << Yellow_Text(" %") << ", double-counted between window and tail."
        << Foam::endl;

    if(overlap > 0.5)
    {
        if(Plus::processor::isMaster())
        {
            WARNING
                << "The Basset window now overruns t_win by "
                << 100.0*overlap << " %. dt has grown enough that one step is a "
                << "large part of the window, so the split point is poorly "
                << "resolved. Raise bassetWindowSteps or cap maxDeltaT."
                << END_WARNING;
        }
    }
}

Foam::scalar basset::setWindowWeights
(
    const std::vector<Foam::scalar>& dt,
    Foam::scalar C
)const
{
    wWin_.assign(dt.size() + 1, 0.0);

    Foam::scalar xa = 0.0;

    for(size_t n=1; n<=dt.size(); ++n)
    {
        const Foam::scalar step = dt[n-1];

        if(step <= 0) break;

        const Foam::scalar xb = xa + step;
        Foam::scalar M0 = 2.0*(Foam::sqrt(xb) - Foam::sqrt(xa));
        Foam::scalar M1 = (2.0/3.0)*(Foam::pow(xb, 1.5) - Foam::pow(xa, 1.5));

        if(C > 0)
        {
            const Foam::scalar psi = meiAdrianPsi(0.5*(xa + xb), C);
            M0 *= psi;
            M1 *= psi;
        }

        wWin_[n]     += (M1 - xa*M0)/step;
        wWin_[n - 1] += (xb*M0 - M1)/step;

        xa = xb;
    }

    return xa;
}

void basset::beginCoupling(size_t nParLocal)
{
    if(!active_) return;

    uEffLoc_.assign(nParLocal, realx3(0.0));
    ampLoc_.assign (nParLocal, 0.0);
    cbLoc_.assign  (nParLocal, 0.0);
    dpLoc_.assign  (nParLocal, 0.0);
    nuLoc_.assign  (nParLocal, 0.0);
    fExplLoc_.assign(nParLocal, realx3(0.0));
    mbLoc_.assign  (nParLocal, 0.0);
}

void basset::commitPending
(
    const span<uint32>& ids,
    const span<realx3>& vel
)
{
    const size_t nPaired = std::min(ids.size(), vel.size());

    std::unordered_map<uint32, size_t> nowIndex;
    nowIndex.reserve(nPaired*2);

    for(size_t i=0; i<nPaired; ++i) nowIndex[ids[i]] = i;

    const Foam::scalar dt = dtPending_;

    Foam::scalar decay[nTail], tvi[nTail], amp[nTail];

    for(int i=0; i<nTail; ++i)
    {
        tvi[i]   = tvH[i]*twin_;
        decay[i] = Foam::exp(-dt/(2.0*tvi[i]));
        amp[i]   = 2.0*Foam::sqrt(Foam::constant::mathematical::e*tvi[i]);
    }

    Foam::scalar maxOverlap = 0.0;

    std::unordered_map<uint32, state> next;
    next.reserve(idPend_.size()*2);

    Foam::scalar maxFB = 0.0;

    const size_t nPend = std::min(idPend_.size(), cbPend_.size());

    for(size_t j=0; j<nPend; ++j)
    {
        const uint32 id = idPend_[j];

        auto found = nowIndex.find(id);
        if(found == nowIndex.end()) continue;

        auto old = hist_.find(id);

        if(cbPend_[j] <= 0)
        {
            if(old != hist_.end()) next.emplace(id, old->second);
            continue;
        }

        const realx3& v = vel[found->second];

        const Foam::vector uEff
        (
            uEffPend_[j].x(), uEffPend_[j].y(), uEffPend_[j].z()
        );

        const Foam::vector fNew =
            ampPend_[j]*(uEff - Foam::vector(v.x(), v.y(), v.z()));

        state st;

        if(old == hist_.end())
        {
            const std::vector<Foam::scalar> even(maxWindow_, dt);

            const int nw = Foam::max(1, windowSpan(even));

            st.f.assign (nw + 2, fNew);
            st.dt.assign(nw + 1, dt);
            for(int k=0; k<nTail; ++k) st.F[k] = Foam::vector::zero;
        }
        else
        {
            st = old->second;
            st.f.insert (st.f.begin(),  fNew);
            st.dt.insert(st.dt.begin(), dt);

            for(int k=0; k<nTail; ++k) st.F[k] = decay[k]*old->second.F[k];

            const int nwPrev = static_cast<int>(old->second.dt.size()) - 1;

            const int nw = Foam::min
            (
                Foam::min(windowSpan(st.dt), maxWindow_),
                nwPrev + 1
            );

            Foam::scalar xa = 0.0;
            for(int m=0; m<nw; ++m) xa += st.dt[m];

            for(int m=nw; m<=nwPrev; ++m)
            {
                const Foam::scalar step = st.dt[m];

                if(step <= 0) break;

                const Foam::scalar xb = xa + step;
                const Foam::vector g  = (st.f[m] - st.f[m+1])/step;

                for(int k=0; k<nTail; ++k)
                {
                    st.F[k] += amp[k]*g
                             *( Foam::exp(-xa/(2.0*tvi[k]))
                              - Foam::exp(-xb/(2.0*tvi[k])) );
                }

                xa = xb;
            }

            const Foam::vector fOldest  = st.f.back();
            const Foam::scalar dtOldest = st.dt.back();

            st.f.resize (nw + 2, fOldest);
            st.dt.resize(nw + 1, dtOldest);
        }

        st.u = uEff;

        const Foam::scalar Ck = kernelC
        (
            Foam::mag(fNew),
            j < dpPend_.size() ? dpPend_[j] : 0.0,
            j < nuPend_.size() ? nuPend_[j] : 0.0
        );

        if(Ck > 0)
        {
            minCrossover_ = Foam::min(minCrossover_, Foam::pow(Ck, -2.0/3.0));
        }

        Foam::scalar aT[nTail];
        tailWeights(Ck, aT);

        Foam::vector fb = Foam::vector::zero;
        for(int k=0; k<nTail; ++k) fb += aT[k]*st.F[k];

        const int nw = static_cast<int>(st.dt.size()) - 1;

        const Foam::scalar edge = setWindowWeights
        (
            std::vector<Foam::scalar>(st.dt.begin(), st.dt.begin() + nw),
            Ck
        );

        maxOverlap = Foam::max(maxOverlap, edge/twin_ - 1.0);

        for(int m=0; m<=nw; ++m)
        {
            if(!(st.dt[m] > 0)) continue;

            fb += wWin_[m]*((st.f[m] - st.f[m+1])/st.dt[m]);
        }

        maxFB = Foam::max(maxFB, Foam::mag(cbPend_[j]*fb));

        next.emplace(id, st);
    }

    hist_.swap(next);

    Foam::Info
        << Blue_Text("Basset history: ")
        << "|F_B|max = " << Yellow_Text(maxFB) << " N";

    if(kernel_ == historyKernel::MeiAdrian && minCrossover_ < GREAT)
    {
        Foam::Info
            << ", finite-Re crossover t_c min = " << Yellow_Text(minCrossover_)
            << " s (" << minCrossover_/Foam::max(twin_, SMALL) << " windows)";
    }

    Foam::Info<< Foam::endl;

    minCrossover_ = GREAT;

    reportOverlap(maxOverlap);

    tCommitted_ = tPending_;
    started_    = true;
    hasPending_ = false;
}

void basset::commitNow()
{
    if(!active_) return;

    auto& pds = semiImplicitCoupling(uCS_).procDEM();

    pds.getDataFromDEM();

    if(!Plus::processor::isMaster()) return;
    if(!hasPending_) return;

    auto ids = pds.particleIdAllMaster();
    auto vel = pds.particlesVelocityAllMaster();

    if(ids.size() < vel.size()) return;

    commitPending(ids, vel);

    if(uCS_.cMesh().mesh().time().writeTime()) writeState();
}

void basset::evaluate(Foam::scalar tNow, Foam::scalar dtNow)
{
    if(!active_) return;

    everEvaluated_ = true;

    auto& pMap = const_cast<particleMapping&>(uCS_.parMapping());
    auto& pds  = semiImplicitCoupling(uCS_).procDEM();

    const size_t nAll = pds.particlesCenterMassAllMaster().size();

    uEffAll_.assign(nAll, realx3(0.0));
    ampAll_.assign (nAll, 0.0);
    cbAll_.assign  (nAll, 0.0);
    dpAll_.assign  (nAll, 0.0);
    nuAll_.assign  (nAll, 0.0);

    {
        span<realx3> mineU(uEffLoc_.data(), uEffLoc_.size());
        span<realx3> allU (uEffAll_.data(), uEffAll_.size());
        span<real>   mineA(ampLoc_.data(),  ampLoc_.size());
        span<real>   allA (ampAll_.data(),  ampAll_.size());
        span<real>   mineC(cbLoc_.data(),   cbLoc_.size());
        span<real>   allC (cbAll_.data(),   cbAll_.size());

        span<real>   mineD(dpLoc_.data(),   dpLoc_.size());
        span<real>   allD (dpAll_.data(),   dpAll_.size());
        span<real>   mineN(nuLoc_.data(),   nuLoc_.size());
        span<real>   allN (nuAll_.data(),   nuAll_.size());

        if(!pMap.realx3ScatteredComm().collectSum(mineU, allU) ||
           !pMap.realScatteredComm().collectSum(mineA, allA)   ||
           !pMap.realScatteredComm().collectSum(mineC, allC)   ||
           !pMap.realScatteredComm().collectSum(mineD, allD)   ||
           !pMap.realScatteredComm().collectSum(mineN, allN))
        {
            fatalErrorInFunction<<"basset: collect to master failed"<<endl;
            Plus::processor::abort(0);
        }
    }

    fExplAll_.assign(nAll, realx3(0.0));
    mbAll_.assign  (nAll, 0.0);

    if(Plus::processor::isMaster())
    {
        auto ids = pds.particleIdAllMaster();
        auto vel = pds.particlesVelocityAllMaster();

        if(ids.size() < nAll || vel.size() != nAll)
        {
            fatalErrorInFunction
                << "basset: the DEM system does not publish a usable particle id "
                   "per particle (" << ids.size() << " ids, " << vel.size()
                << " velocities, " << nAll << " particles)."
                << endl;
            Plus::processor::abort(0);
        }

        if(hasPending_ && tNow > tPending_)
        {
            commitPending(ids, vel);
        }

        if(!started_)
        {
            tCommitted_ = tNow - dtNow;
        }

        const Foam::scalar dt = tNow - tCommitted_;

        if(dt > 0)
        {
            if(twin_ <= 0) twin_ = nWindow_*dt;

            for(size_t i=0; i<nAll; ++i)
            {
                const Foam::scalar cB  = cbAll_[i];
                const Foam::scalar amp = ampAll_[i];

                if(cB <= 0) continue;

                auto old = hist_.find(ids[i]);

                std::vector<Foam::scalar> steps;
                steps.reserve(maxWindow_ + 2);
                steps.push_back(dt);

                if(old != hist_.end())
                {
                    steps.insert
                    (
                        steps.end(),
                        old->second.dt.begin(),
                        old->second.dt.end()
                    );
                }

                int nw = Foam::min(windowSpan(steps), maxWindow_);

                if(old != hist_.end())
                {
                    nw = Foam::min(nw, static_cast<int>(old->second.f.size()) - 1);
                }

                nw = Foam::max(nw, 1);

                while(static_cast<int>(steps.size()) < nw + 1)
                {
                    steps.push_back(steps.back());
                }

                const Foam::vector uE
                (
                    uEffAll_[i].x(), uEffAll_[i].y(), uEffAll_[i].z()
                );
                const Foam::vector vP(vel[i].x(), vel[i].y(), vel[i].z());
                const Foam::scalar C = kernelC
                (
                    amp*Foam::mag(uE - vP), dpAll_[i], nuAll_[i]
                );

                setWindowWeights
                (
                    std::vector<Foam::scalar>(steps.begin(), steps.begin() + nw),
                    C
                );

                mbAll_[i] = static_cast<real>(cB*amp*wWin_[0]);

                if(old == hist_.end()) continue;

                const state& st = old->second;

                const Foam::vector uNow
                (
                    uEffAll_[i].x(), uEffAll_[i].y(), uEffAll_[i].z()
                );

                Foam::vector win = wWin_[0]*(amp*(uNow - st.u)/dt);

                for(int m=1; m<=nw; ++m)
                {
                    win += wWin_[m]*((st.f[m-1] - st.f[m])/steps[m]);
                }

                Foam::scalar aT[nTail];
                tailWeights(C, aT);

                Foam::vector tail = Foam::vector::zero;
                for(int k=0; k<nTail; ++k) tail += aT[k]*st.F[k];

                const Foam::vector fe = cB*(win + tail);

                fExplAll_[i] = realx3(fe.x(), fe.y(), fe.z());
            }
        }

        idPend_.assign(ids.begin(), ids.begin() + nAll);
        uEffPend_ = uEffAll_;
        ampPend_  = ampAll_;
        cbPend_   = cbAll_;
        dpPend_   = dpAll_;
        nuPend_   = nuAll_;
        tPending_ = tNow;
        dtPending_= dt;
        hasPending_ = (dt > 0);
    }

    {
        span<realx3> allF (fExplAll_.data(), fExplAll_.size());
        span<realx3> mineF(fExplLoc_.data(), fExplLoc_.size());
        span<real>   allM (mbAll_.data(),    mbAll_.size());
        span<real>   mineM(mbLoc_.data(),    mbLoc_.size());

        if(!pMap.realx3ScatteredComm().distribute(allF, mineF) ||
           !pMap.realScatteredComm().distribute(allM, mineM))
        {
            fatalErrorInFunction<<"basset: distribute to ranks failed"<<endl;
            Plus::processor::abort(0);
        }
    }
}

}
