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

#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <cmath>

#include <fstream>
#include <iomanip>

#include "rotationalHistory.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "momentumSemiImplicitSphereUnresolvedCouplingSystem.hpp"

namespace pFlow::coupling
{

static inline Foam::scalar erfcxSqrt(const Foam::scalar x)
{
    if(x < 25.0)
    {
        return Foam::exp(x)*Foam::erfc(Foam::sqrt(x));
    }

    const Foam::scalar t = 0.5/x;
    Foam::scalar ser = 1.0, term = 1.0;

    for(int n=1; n<=11; ++n)
    {
        term *= -(2.0*n - 1.0)*t;
        ser  += term;
    }

    return ser/Foam::sqrt(Foam::constant::mathematical::pi*x);
}

static inline Foam::scalar kernelK(const Foam::scalar x)
{
    if(x < 25.0)
    {
        return 1.0/Foam::sqrt(Foam::constant::mathematical::pi*x) - erfcxSqrt(x);
    }

    const Foam::scalar t = 0.5/x;
    Foam::scalar ser = 0.0, term = 1.0;

    for(int n=1; n<=11; ++n)
    {
        term *= -(2.0*n - 1.0)*t;
        ser  -= term;
    }

    return ser/Foam::sqrt(Foam::constant::mathematical::pi*x);
}

static inline Foam::scalar kernelI0(const Foam::scalar X)
{
    return 1.0 - erfcxSqrt(X);
}

static inline Foam::scalar kernelI1(const Foam::scalar X)
{
    const Foam::scalar pi = Foam::constant::mathematical::pi;

    if(X < 1.0e-4)
    {
        const Foam::scalar s = Foam::sqrt(X);
        return (2.0/(3.0*Foam::sqrt(pi)))*X*s - X*X;
    }

    return (1.0 - X)*erfcxSqrt(X) - 1.0 + 2.0*Foam::sqrt(X/pi);
}

static inline Foam::scalar expM0(const Foam::scalar z)
{
    if(Foam::mag(z) < 1.0e-3)
    {
        return 1.0 - 0.5*z + z*z/6.0 - z*z*z/24.0;
    }
    return (1.0 - Foam::exp(-z))/z;
}

rotationalHistory::rotationalHistory
(
    const unresolvedCouplingSystem& uCS,
    const Foam::dictionary&         dict
)
:
    uCS_(uCS)
{
    active_ = dict.getOrDefault<Foam::Switch>
        ("rotationalHistory", Foam::Switch(false));

    if(!active_) return;

    nWindow_ = dict.getOrDefault<Foam::label>("rotationalWindowSteps", 20);

    if(nWindow_ < 2)
    {
        fatalErrorInFunction
            << "rotationalWindowSteps must be at least 2, got " << nWindow_
            << " in " << dict.name() << endl;
        Plus::processor::abort(0);
    }

    setQuadrature();
    readState();

    Foam::Info
        << "    Rotational history couple is " << Green_Text("active")
        << ": Feuillebois & Lasek, QJMAM 31 (1978) 435.\n"
           "      Exact kernel over a window of " << nWindow_
        << " steps, " << nModes << " quadrature exponentials beyond it.\n"
           "      Window nodes follow each particle's own step history, so a "
           "changing\n      dt moves the window edge rather than stretching "
           "the memory.\n\n";
}

void rotationalHistory::setQuadrature()
{
    const Foam::scalar pi = Foam::constant::mathematical::pi;
    const int n = nModes;

    for(int i=1; i<=n; ++i)
    {
        Foam::scalar x = Foam::cos(pi*(i - 0.25)/(n + 0.5));
        Foam::scalar dp = 0.0;

        for(int it=0; it<100; ++it)
        {
            Foam::scalar p0 = 1.0, p1 = 0.0;

            for(int j=1; j<=n; ++j)
            {
                const Foam::scalar p2 = p1;
                p1 = p0;
                p0 = ((2.0*j - 1.0)*x*p1 - (j - 1.0)*p2)/j;
            }

            dp = n*(x*p0 - p1)/(x*x - 1.0);

            const Foam::scalar dx = -p0/dp;
            x += dx;

            if(Foam::mag(dx) < 1.0e-15) break;
        }

        const Foam::scalar w  = 2.0/((1.0 - x*x)*dp*dp);
        const Foam::scalar th = 0.25*pi*(x + 1.0);
        const Foam::scalar wq = 0.25*pi*w;
        const Foam::scalar t2 = Foam::sqr(Foam::tan(th));

        lamQ_[i-1] = t2;
        ampQ_[i-1] = (2.0/pi)*wq*t2;
    }
}

void rotationalHistory::checkQuadrature(const Foam::scalar xw)const
{
    if(!active_ || reportedQuad_) return;
    reportedQuad_ = true;

    const Foam::scalar xMax = 1.0e3;

    Foam::scalar worst = 0.0;
    Foam::scalar atX   = xw;

    for(int m=0; m<60; ++m)
    {
        const Foam::scalar x = xw*Foam::pow(1.4, m);
        if(x > xMax) break;

        Foam::scalar q = 0.0;
        for(int k=0; k<nModes; ++k) q += ampQ_[k]*Foam::exp(-x*lamQ_[k]);

        const Foam::scalar e = Foam::mag(q/kernelK(x) - 1.0);

        if(e > worst) { worst = e; atX = x; }
    }

    Foam::scalar unit = 0.0;
    for(int k=0; k<nModes; ++k) unit += ampQ_[k]/lamQ_[k];

    Foam::Info
        << Blue_Text("Rotational history: ")
        << "tail quadrature worst error " << Yellow_Text(worst)
        << " at x = " << Yellow_Text(atX)
        << " over lags x = " << xw << " to " << xMax
        << ", kernel mass " << Yellow_Text(unit) << " (exact 1)" << Foam::endl;

    if(worst > 0.05)
    {
        if(Plus::processor::isMaster())
        {
            WARNING
                << "The " << nModes << "-node quadrature of the rotational kernel "
                << "is " << 100.0*worst << " % off at x = " << atX << ". The "
                << "window edge sits at x = " << xw << " = " << nWindow_
                << " steps of dt/(a^2/nu); a window that short leaves the tail to "
                << "cover lags the nodes do not reach. Raise "
                << "rotationalWindowSteps." << END_WARNING;
        }
    }
}

Foam::scalar rotationalHistory::setWindowWeights
(
    const std::vector<Foam::scalar>& dt,
    const Foam::scalar               tnu
)const
{
    wWin_.assign(nWindow_ + 1, 0.0);

    if(!(tnu > 0)) return 0.0;

    Foam::scalar xa  = 0.0;
    Foam::scalar I0a = 0.0, I1a = 0.0;

    for(int n=1; n<=nWindow_; ++n)
    {
        const Foam::scalar step =
            (size_t(n-1) < dt.size() && dt[n-1] > 0) ? dt[n-1] : 0.0;

        if(step <= 0) break;

        const Foam::scalar xb  = xa + step/tnu;
        const Foam::scalar I0b = kernelI0(xb);
        const Foam::scalar I1b = kernelI1(xb);

        const Foam::scalar M0 = I0b - I0a;
        const Foam::scalar M1 = I1b - I1a;
        const Foam::scalar dx = xb - xa;

        if(!(dx > 0)) break;

        wWin_[n]     += tnu*(M1 - xa*M0)/dx;
        wWin_[n - 1] += tnu*(xb*M0 - M1)/dx;

        xa  = xb;
        I0a = I0b;
        I1a = I1b;
    }

    return xa;
}

Foam::fileName rotationalHistory::statePath(const Foam::word& timeName)const
{
    const Foam::Time& rt = uCS_.cMesh().mesh().time();
    return rt.globalPath()/timeName/"uniform"/"rotationalHistory";
}

void rotationalHistory::writeState()const
{
    if(!active_ || !Plus::processor::isMaster()) return;

    const Foam::Time& rt = uCS_.cMesh().mesh().time();
    const Foam::fileName f = statePath(rt.timeName());

    Foam::mkDir(f.path());

    std::ofstream out(f.c_str());
    if(!out)
    {
        Foam::Warning<<"rotationalHistory: cannot write "<<f<<Foam::endl;
        return;
    }

    out << std::setprecision(17);
    out << "nModes "    << nModes     << "\n";
    out << "nWindow "   << nWindow_   << "\n";
    out << "tCommitted "<< tCommitted_<< "\n";
    out << "nParticles "<< hist_.size() << "\n";

    for(const auto& kv : hist_)
    {
        out << kv.first;

        for(int k=0; k<nModes; ++k)
        {
            out << ' ' << kv.second.F[k].x() << ' ' << kv.second.F[k].y()
                << ' ' << kv.second.F[k].z();
        }

        for(int n=0; n<nWindow_+2; ++n)
        {
            const Foam::vector& wv = kv.second.w[n];
            out << ' ' << wv.x() << ' ' << wv.y() << ' ' << wv.z();
        }

        for(int n=0; n<nWindow_+1; ++n) out << ' ' << kv.second.dt[n];

        out << ' ' << kv.second.wAmb.x() << ' ' << kv.second.wAmb.y()
            << ' ' << kv.second.wAmb.z();
        out << "\n";
    }
}

void rotationalHistory::readState()
{
    if(!active_ || !Plus::processor::isMaster()) return;

    const Foam::Time& rt = uCS_.cMesh().mesh().time();
    const Foam::fileName f = statePath(rt.timeName());

    std::ifstream in(f.c_str());
    if(!in) return;

    std::string key;
    int    nm = 0, nw = 0;
    double tc = 0;
    long   np = 0;

    in >> key >> nm >> key >> nw >> key >> tc >> key >> np;

    if(nm != nModes || nw != nWindow_)
    {
        Foam::Warning
            << "rotationalHistory: " << f << " was written with nModes = "
            << nm << ", nWindow = " << nw << " but this run uses " << nModes
            << ", " << nWindow_ << ". Ignoring it and starting the memory "
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

        for(int k=0; k<nModes; ++k)
        {
            in >> x >> y >> z; st.F[k] = Foam::vector(x,y,z);
        }

        st.w.resize(nWindow_ + 2);
        for(int m=0; m<nWindow_+2; ++m)
        {
            in >> x >> y >> z; st.w[m] = Foam::vector(x,y,z);
        }

        st.dt.resize(nWindow_ + 1);
        for(int m=0; m<nWindow_+1; ++m) in >> st.dt[m];

        in >> x >> y >> z; st.wAmb = Foam::vector(x,y,z);

        hist_.emplace(id, st);
    }

    tCommitted_ = tc;
    started_    = true;

    Foam::Info
        << "    Rotational history resumed from " << f.name() << ": "
        << hist_.size() << " particles, committed at t = " << tCommitted_
        << " s.\n";
}

void rotationalHistory::beginCoupling(const size_t nParLocal)
{
    if(!active_) return;

    wAmbLoc_.assign (nParLocal, realx3(0.0));
    tnuLoc_.assign  (nParLocal, 0.0);
    crLoc_.assign   (nParLocal, 0.0);
    kwLoc_.assign   (nParLocal, 0.0);
    tExplLoc_.assign(nParLocal, realx3(0.0));
    iaLoc_.assign   (nParLocal, 0.0);
}

void rotationalHistory::commitPending
(
    const span<uint32>& ids,
    const span<realx3>& rVel
)
{
    const size_t nPaired = std::min(ids.size(), rVel.size());

    std::unordered_map<uint32, size_t> nowIndex;
    nowIndex.reserve(nPaired*2);

    for(size_t i=0; i<nPaired; ++i) nowIndex[ids[i]] = i;

    const Foam::scalar dt = dtPending_;

    std::unordered_map<uint32, state> next;
    next.reserve(idPend_.size()*2);

    Foam::scalar maxT     = 0.0;
    Foam::scalar maxRatio = 0.0;

    const size_t nPend = std::min(idPend_.size(), kwPend_.size());

    for(size_t j=0; j<nPend; ++j)
    {
        const uint32 id = idPend_[j];

        auto found = nowIndex.find(id);
        if(found == nowIndex.end()) continue;

        auto old = hist_.find(id);

        if(crPend_[j] <= 0 || tnuPend_[j] <= 0)
        {
            if(old != hist_.end()) next.emplace(id, old->second);
            continue;
        }

        const realx3& wp = rVel[found->second];

        const Foam::vector wNew =
            Foam::vector(wp.x(), wp.y(), wp.z())
          - Foam::vector(wAmbPend_[j].x(), wAmbPend_[j].y(), wAmbPend_[j].z());

        state st;

        const Foam::scalar tnu = tnuPend_[j];

        if(old == hist_.end())
        {
            st.w.assign (nWindow_ + 2, wNew);
            st.dt.assign(nWindow_ + 1, dt);
            for(int k=0; k<nModes; ++k) st.F[k] = Foam::vector::zero;
        }
        else
        {
            st = old->second;
            st.w.insert(st.w.begin(), wNew);
            st.w.resize(nWindow_ + 2);
            st.dt.insert(st.dt.begin(), dt);
            st.dt.resize(nWindow_ + 1);

            const Foam::scalar stepN =
                st.dt[nWindow_] > 0 ? st.dt[nWindow_] : dt;

            Foam::scalar xw = 0.0;
            for(int n=0; n<nWindow_; ++n) xw += st.dt[n]/tnu;

            const Foam::scalar h = stepN/tnu;

            const Foam::vector gN =
                (st.w[nWindow_] - st.w[nWindow_ + 1])/stepN;

            for(int k=0; k<nModes; ++k)
            {
                const Foam::scalar z     = lamQ_[k]*h;
                const Foam::scalar decay = Foam::exp(-lamQ_[k]*dt/tnu);
                const Foam::scalar shift = Foam::exp(-lamQ_[k]*xw);

                const Foam::vector di = tnu*h*shift*expM0(z)*gN;

                st.F[k] = decay*old->second.F[k] + di;
            }
        }

        st.wAmb = Foam::vector(wAmbPend_[j].x(), wAmbPend_[j].y(),
                               wAmbPend_[j].z());

        setWindowWeights(st.dt, tnu);

        Foam::vector th = Foam::vector::zero;
        for(int k=0; k<nModes; ++k) th += ampQ_[k]*st.F[k];

        for(int n=0; n<=nWindow_; ++n)
        {
            const Foam::scalar sn = st.dt[n] > 0 ? st.dt[n] : dt;

            if(!(sn > 0)) continue;

            th += wWin_[n]*((st.w[n] - st.w[n+1])/sn);
        }

        const Foam::scalar magT = Foam::mag(crPend_[j]*th);

        maxT = Foam::max(maxT, magT);

        const Foam::scalar magS = kwPend_[j]*Foam::mag(st.w[0]);

        if(magS > 0)
        {
            maxRatio = Foam::max(maxRatio, magT/magS);
        }

        next.emplace(id, st);
    }

    hist_.swap(next);

    Foam::Info
        << Blue_Text("Rotational history: ")
        << "|T_R|max = " << Yellow_Text(maxT) << " N m, "
        << "largest fraction of the steady couple " << Yellow_Text(maxRatio)
        << Foam::endl;

    tCommitted_ = tPending_;
    started_    = true;
    hasPending_ = false;
}

void rotationalHistory::commitNow()
{
    if(!active_) return;

    auto& pds = semiImplicitCoupling(uCS_).procDEM();

    pds.getDataFromDEM();

    if(!Plus::processor::isMaster()) return;
    if(!hasPending_) return;

    auto ids  = pds.particleIdAllMaster();
    auto rVel = pds.particlesRVelocityAllMaster();

    if(ids.size() < rVel.size()) return;

    commitPending(ids, rVel);

    if(uCS_.cMesh().mesh().time().writeTime()) writeState();
}

void rotationalHistory::evaluate
(
    const Foam::scalar tNow,
    const Foam::scalar dtNow
)
{
    if(!active_) return;

    everEvaluated_ = true;

    auto& pMap = const_cast<particleMapping&>(uCS_.parMapping());
    auto& pds  = semiImplicitCoupling(uCS_).procDEM();

    const size_t nAll = pds.particlesCenterMassAllMaster().size();

    wAmbAll_.assign(nAll, realx3(0.0));
    tnuAll_.assign (nAll, 0.0);
    crAll_.assign  (nAll, 0.0);
    kwAll_.assign  (nAll, 0.0);

    {
        span<realx3> mineW(wAmbLoc_.data(), wAmbLoc_.size());
        span<realx3> allW (wAmbAll_.data(), wAmbAll_.size());
        span<real>   mineT(tnuLoc_.data(),  tnuLoc_.size());
        span<real>   allT (tnuAll_.data(),  tnuAll_.size());
        span<real>   mineC(crLoc_.data(),   crLoc_.size());
        span<real>   allC (crAll_.data(),   crAll_.size());
        span<real>   mineK(kwLoc_.data(),   kwLoc_.size());
        span<real>   allK (kwAll_.data(),   kwAll_.size());

        if(!pMap.realx3ScatteredComm().collectSum(mineW, allW) ||
           !pMap.realScatteredComm().collectSum(mineT, allT)   ||
           !pMap.realScatteredComm().collectSum(mineC, allC)   ||
           !pMap.realScatteredComm().collectSum(mineK, allK))
        {
            fatalErrorInFunction
                <<"rotationalHistory: collect to master failed"<<endl;
            Plus::processor::abort(0);
        }
    }

    tExplAll_.assign(nAll, realx3(0.0));
    iaAll_.assign   (nAll, 0.0);

    if(Plus::processor::isMaster())
    {
        auto ids  = pds.particleIdAllMaster();
        auto rVel = pds.particlesRVelocityAllMaster();

        if(ids.size() < nAll || rVel.size() != nAll)
        {
            fatalErrorInFunction
                << "rotationalHistory: the DEM system does not publish a usable "
                   "particle id and rotational velocity per particle ("
                << ids.size() << " ids, " << rVel.size()
                << " rotational velocities, " << nAll << " particles)."
                << endl;
            Plus::processor::abort(0);
        }

        if(hasPending_ && tNow > tPending_)
        {
            commitPending(ids, rVel);
        }

        if(!started_)
        {
            tCommitted_ = tNow - dtNow;
        }

        const Foam::scalar dt = tNow - tCommitted_;

        if(dt > 0)
        {
            for(size_t i=0; i<nAll; ++i)
            {
                const Foam::scalar cR  = crAll_[i];
                const Foam::scalar tnu = tnuAll_[i];

                if(cR <= 0 || tnu <= 0) continue;

                auto old = hist_.find(ids[i]);

                std::vector<Foam::scalar> steps(nWindow_ + 1, dt);

                if(old != hist_.end())
                {
                    for(int n=1; n<=nWindow_; ++n)
                    {
                        steps[n] = (size_t(n-1) < old->second.dt.size()
                                    && old->second.dt[n-1] > 0)
                                 ? old->second.dt[n-1] : dt;
                    }
                }

                const Foam::scalar xw = setWindowWeights(steps, tnu);

                checkQuadrature(xw);

                iaAll_[i] = static_cast<real>(cR*wWin_[0]);

                if(old == hist_.end()) continue;

                const state& st = old->second;

                const Foam::vector wAmbNow
                (
                    wAmbAll_[i].x(), wAmbAll_[i].y(), wAmbAll_[i].z()
                );

                Foam::vector win = wWin_[0]*(-(wAmbNow - st.wAmb)/dt);

                for(int n=1; n<=nWindow_; ++n)
                {
                    win += wWin_[n]*((st.w[n-1] - st.w[n])/steps[n]);
                }

                Foam::vector tail = Foam::vector::zero;
                for(int k=0; k<nModes; ++k) tail += ampQ_[k]*st.F[k];

                const Foam::vector te = -cR*(win + tail);

                tExplAll_[i] = realx3(te.x(), te.y(), te.z());
            }
        }

        idPend_.assign(ids.begin(), ids.begin() + nAll);
        wAmbPend_ = wAmbAll_;
        tnuPend_  = tnuAll_;
        crPend_   = crAll_;
        kwPend_   = kwAll_;
        tPending_ = tNow;
        dtPending_= dt;
        hasPending_ = (dt > 0);
    }

    {
        span<realx3> allT (tExplAll_.data(), tExplAll_.size());
        span<realx3> mineT(tExplLoc_.data(), tExplLoc_.size());
        span<real>   allI (iaAll_.data(),    iaAll_.size());
        span<real>   mineI(iaLoc_.data(),    iaLoc_.size());

        if(!pMap.realx3ScatteredComm().distribute(allT, mineT) ||
           !pMap.realScatteredComm().distribute(allI, mineI))
        {
            fatalErrorInFunction
                <<"rotationalHistory: distribute to ranks failed"<<endl;
            Plus::processor::abort(0);
        }
    }
}

}
