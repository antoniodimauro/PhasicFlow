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

// Costa et al. (2015) Phys. Rev. E 92, 053012

#ifndef __lubrication_hpp__
#define __lubrication_hpp__

#include "types.hpp"
#include "dictionary.hpp"

namespace pFlow::lubricationModels
{

INLINE_FUNCTION_HD
real lubricationLambda(const real eps)
{
    const real lne = log(eps);

    return static_cast<real>(1.0)/eps
         - static_cast<real>(1.0/5.0)*lne
         - static_cast<real>(1.0/21.0)*eps*lne;
}

class lubrication
{
private:

    bool active_ = false;

    real mu_ = 0.0;

    real epsOut_ = 1.0;

    real roughness_ = 0.0;

    real epsMin_ = 1.0e-4;

    real cflLub_ = 0.5;

    static constexpr real aMaxDt = static_cast<real>(0.158);

public:

    INLINE_FUNCTION_HD
    lubrication() = default;

    INLINE_FUNCTION_HD
    lubrication(const lubrication&) = default;

    INLINE_FUNCTION_HD
    lubrication& operator=(const lubrication&) = default;

    INLINE_FUNCTION_HD
    ~lubrication() = default;

    bool read(const dictionary& dict)
    {
        if( !dict.containsDictionay("lubrication") )
        {
            active_ = false;
            return true;
        }

        const auto& lDict = dict.subDict("lubrication");

        active_ = lDict.getValOrSet<Logical>("active", Logical(true))();

        if(!active_) return true;

        mu_        = lDict.getVal<real>("fluidViscosity");
        roughness_ = lDict.getVal<real>("roughness");
        epsOut_    = lDict.getValOrSet<real>("epsOut", static_cast<real>(1.0));
        epsMin_    = lDict.getValOrSet<real>("epsMin", static_cast<real>(1.0e-4));
        cflLub_    = lDict.getValOrSet<real>("cflLub", static_cast<real>(0.5));

        if(mu_ <= 0)
        {
            fatalErrorInFunction<<
            "fluidViscosity in the lubrication dictionary must be positive, "
            "but it is "<<mu_<<".\n";
            return false;
        }

        if(roughness_ <= 0)
        {
            fatalErrorInFunction<<
            "roughness in the lubrication dictionary must be positive: it is "
            "the physical asperity height in metres, and it is what stops the "
            "1/eps resistance from diverging. Costa et al. (2015) use 1e-5 of "
            "a radius for smooth steel and 1e-3 for rough surfaces.\n";
            return false;
        }

        if(epsOut_ <= 0)
        {
            fatalErrorInFunction<<
            "epsOut in the lubrication dictionary must be positive.\n";
            return false;
        }

        return true;
    }

    INLINE_FUNCTION_HD
    bool active()const
    {
        return active_;
    }

    INLINE_FUNCTION_HD
    real epsOut()const
    {
        return epsOut_;
    }

    INLINE_FUNCTION_HD
    real cflLub()const
    {
        return cflLub_;
    }

    static INLINE_FUNCTION_HD
    real reducedRadius(const real Ri, const real Rj)
    {
        return Ri*Rj/max(Ri+Rj, smallValue);
    }

    INLINE_FUNCTION_HD
    real range(const real Ri, const real Rj)const
    {
        return epsOut_*reducedRadius(Ri, Rj);
    }

    INLINE_FUNCTION_HD
    real maxCoeffSum(const real mass, const real dt)const
    {
        if(mass <= 0 || dt <= 0) return largeValue;

        return cflLub_*aMaxDt*mass/(static_cast<real>(2.0)*dt);
    }

    INLINE_FUNCTION_HD
    real minResolvedEps
    (
        const real Rl,
        const real mMin,
        const real dt
    )const
    {
        const real cMax = maxCoeffSum(mMin, dt);

        if(cMax <= 0 || cMax >= largeValue) return 0.0;

        return static_cast<real>(6.0)*Pi*mu_*Rl/cMax;
    }

    INLINE_FUNCTION_HD
    real coefficient(const real Ri, const real Rj, const real gap)const
    {
        if(!active_) return 0.0;

        if(gap < 0) return 0.0;

        const real Rl = reducedRadius(Ri, Rj);

        if(Rl <= 0) return 0.0;

        const real eps = gap/Rl;

        if(eps >= epsOut_) return 0.0;

        const real epsSig = max(roughness_/Rl, epsMin_);

        if(epsOut_ <= epsSig) return 0.0;

        const real e   = max(eps, epsSig);
        const real lam = lubricationLambda(e) - lubricationLambda(epsOut_);

        if(lam <= 0) return 0.0;

        return static_cast<real>(6.0)*Pi*mu_*Rl*lam;
    }

    INLINE_FUNCTION_HD
    realx3 force
    (
        const real    Ri,
        const real    Rj,
        const real    gap,
        const realx3& Vr,
        const realx3& Nij,
        const real    scale
    )const
    {
        const real c = scale*coefficient(Ri, Rj, gap);

        if(c <= 0) return realx3(0.0);

        const real vrn = dot(Vr, Nij);

        return (-c*vrn)*Nij;
    }
};

}

#endif
