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

// Brenner (1961) Chem. Eng. Sci. 16, 242; Goldman, Cox & Brenner (1967) Chem. Eng. Sci. 22, 637, 653; Ryu & Matsudaira (2010) Chem. Eng. Sci. 65, 4913; Zeng, Najjar, Balachandar & Fischer (2009) Phys. Fluids 21, 033302

#ifndef __wallCorrection_hpp__
#define __wallCorrection_hpp__

#include "OFCompatibleHeader.hpp"

namespace pFlow::coupling
{

inline Foam::scalar wallInterp
(
    const Foam::scalar  e,
    const Foam::scalar* eTab,
    const Foam::scalar* red,
    const int           n
)
{
    if (e <= 0.0)        return red[0];
    if (e >= eTab[n-1])  return red[n-1];

    int k = 1;
    while (k < n-1 && eTab[k] < e) ++k;

    const Foam::scalar t = (e - eTab[k-1])/(eTab[k] - eTab[k-1]);

    return red[k-1] + t*(red[k] - red[k-1]);
}

inline constexpr Foam::scalar wallTabE[9] =
{
    0.0, 0.0993279, 0.2658022, 0.4250960, 0.6480543,
    0.8868189, 0.9566279, 0.9950207, 0.9968085
};

inline Foam::scalar wallShearForce(const Foam::scalar hOverA)
{
    static constexpr Foam::scalar red[9] =
    {
        0.5625,
        0.590972, 0.628663, 0.653970, 0.677567,
        0.694618, 0.698495, 0.700387, 0.700435
    };

    const Foam::scalar e = 1.0/Foam::max(hOverA, Foam::scalar(1));

    return 1.0 + e*wallInterp(e, wallTabE, red, 9);
}

inline Foam::scalar wallShearTorque(const Foam::scalar hOverA)
{
    static constexpr Foam::scalar red[9] =
    {
        0.1875,
        0.193883, 0.153894, 0.128877, 0.094832,
        0.066329, 0.059752, 0.056419, 0.056267
    };

    const Foam::scalar e = 1.0/Foam::max(hOverA, Foam::scalar(1));

    return 1.0 - Foam::pow(e,3)*wallInterp(e, wallTabE, red, 9);
}

inline Foam::scalar wallTransTorque(const Foam::scalar hOverA)
{
    static constexpr Foam::scalar red[9] =
    {
        0.09375,
        0.090143, 0.084463, 0.080916, 0.083054,
        0.119188, 0.173760, 0.348765, 0.389894
    };

    const Foam::scalar e = 1.0/Foam::max(hOverA, Foam::scalar(1));

    return Foam::pow(e,4)*wallInterp(e, wallTabE, red, 9);
}

inline Foam::scalar wallRotForce(const Foam::scalar hOverA)
{
    return (4.0/3.0)*wallTransTorque(hOverA);
}

inline Foam::scalar wallRotTorque(const Foam::scalar hOverA)
{
    static constexpr Foam::scalar red[9] =
    {
        0.3125,
        0.306131, 0.314178, 0.325446, 0.366688,
        0.555893, 0.799136, 1.528316, 1.695482
    };

    const Foam::scalar e = 1.0/Foam::max(hOverA, Foam::scalar(1));

    return 1.0 + Foam::pow(e,3)*wallInterp(e, wallTabE, red, 9);
}

inline Foam::scalar wallDragNormal(const Foam::scalar hOverA, const int nMax = 200)
{
    if (hOverA <= 1.0) return Foam::GREAT;

    const Foam::scalar al = Foam::acosh(hOverA);
    const Foam::scalar sh = Foam::sinh(al);
    const Foam::scalar s2 = Foam::sinh(2.0*al);

    Foam::scalar sum = 0.0;

    for (int n = 1; n <= nMax; ++n)
    {
        const Foam::scalar c  = Foam::scalar(n*(n+1))
                              / Foam::scalar((2*n-1)*(2*n+3));
        const Foam::scalar x  = (2*n + 1)*al;
        const Foam::scalar A  = (2*n + 1)*s2;
        const Foam::scalar B  = Foam::sqr(Foam::scalar(2*n + 1))*sh*sh;

        const Foam::scalar e1 = Foam::exp(-x);
        const Foam::scalar e2 = e1*e1;

        const Foam::scalar nS = -2.0*e2 + (A + 2.0 + B)*e1;
        const Foam::scalar dS =  1.0 + e2 - (2.0 + B)*e1;

        sum += c*(nS/dS);
    }

    return (4.0/3.0)*sh*sum;
}

inline Foam::scalar wallTransForce(const Foam::scalar hOverA)
{
    static constexpr Foam::scalar red[9] =
    {
        0.5625,
        0.594999, 0.653870, 0.724307, 0.875698,
        1.298349, 1.722195, 2.800243, 3.031977
    };

    const Foam::scalar e = 1.0/Foam::max(hOverA, Foam::scalar(1));

    return 1.0 + e*wallInterp(e, wallTabE, red, 9);
}

// Schiller & Naumann (1933)
inline Foam::scalar wallStdDrag(const Foam::scalar Re)
{
    return 1.0 + 0.15*Foam::pow(Foam::max(Re, Foam::scalar(0)), 0.687);
}

// Ryu & Matsudaira (2010) eq. (6)
inline Foam::scalar wallDragNormalRe(const Foam::scalar hOverA, const Foam::scalar Re)
{
    const Foam::scalar sRe = 1.0 + Foam::sqrt(Foam::max(Re, Foam::scalar(0)));

    if (hOverA >= 1.4)
    {
        const Foam::scalar ex = 0.39*Foam::pow(hOverA, -0.85);
        return wallDragNormal(hOverA)/Foam::pow(sRe, ex);
    }

    // Brenner (1961); Cox & Brenner (1967)
    const Foam::scalar lam14 = wallDragNormal(1.4);
    const Foam::scalar ex14  = 0.39*Foam::pow(1.4, -0.85);
    return wallDragNormal(hOverA) - (lam14 - lam14/Foam::pow(sRe, ex14));
}

// Vasseur & Cox (1977) Fig. 4
inline Foam::scalar wallDragNormalGapRe(const Foam::scalar hOverA, const Foam::scalar ReGap)
{
    return 1.0 + (wallDragNormal(hOverA) - 1.0)/(1.0 + Foam::max(ReGap, Foam::scalar(0)));
}

// Zeng et al. (2009) eqs. (16)-(18)
inline Foam::scalar wallTransForceRe(const Foam::scalar hOverA, const Foam::scalar Re)
{
    const Foam::scalar d   = Foam::max(0.5*(hOverA - 1.0), Foam::scalar(1.0e-3));
    const Foam::scalar sd  = Foam::sqrt(d);
    const Foam::scalar c0  = 1.028 - 0.07/(1.0 + 4.0*d*d)
                           - (8.0/15.0)*Foam::log(270.0*d/(135.0 + 256.0*d));
    const Foam::scalar al  = 0.15*(1.0 - Foam::exp(-sd));
    const Foam::scalar be  = 0.687 + 0.313*Foam::exp(-2.0*sd);
    const Foam::scalar R   = Foam::max(Re, Foam::scalar(0));
    return c0*(1.0 + al*Foam::pow(R, be))/wallStdDrag(R);
}

// Zeng et al. (2009) eqs. (11), (13), (14)
inline Foam::scalar wallShearForceRe(const Foam::scalar hOverA, const Foam::scalar Re)
{
    const Foam::scalar d   = Foam::max(0.5*(hOverA - 1.0), Foam::scalar(1.0e-3));
    const Foam::scalar c0  = 1.0 + 0.138*Foam::exp(-2.0*d) + 9.0/(16.0*(1.0 + 2.0*d));
    const Foam::scalar al  = 0.15 - 0.046*(1.0 - 0.16*d*d)*Foam::exp(-0.7*d);
    const Foam::scalar be  = 0.687 + 0.066*(1.0 - 0.76*d*d)*Foam::exp(-Foam::pow(d, 0.9));
    const Foam::scalar R   = Foam::max(Re, Foam::scalar(0));
    return c0*(1.0 + al*Foam::pow(R, be))/wallStdDrag(R);
}

}

#endif
