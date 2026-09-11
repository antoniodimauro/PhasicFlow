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

// Ireland & Desjardins (2017) J. Comput. Phys. 338, 405; Balachandar, Liu & Lakhote (2019) J. Comput. Phys. 376, 160; Pakseresht, Esmaily & Apte (2020) J. Comput. Phys. 420, 109711

#ifndef __selfInduction_hpp__
#define __selfInduction_hpp__

#include "OFCompatibleHeader.hpp"
#include "wallCorrection.hpp"

namespace pFlow::coupling
{

inline constexpr Foam::scalar selfInducedCentreC = 0.79788456080286536;

inline constexpr Foam::scalar selfInducedAverageFactor = 0.70710678118654752;

inline Foam::scalar selfInducedSampleFactor(const Foam::word& fluidVelocityMethod)
{
    if (fluidVelocityMethod == "distribution")
    {
        return selfInducedAverageFactor;
    }

    return 1.0;
}

inline Foam::word fluidVelocityMethod(const Foam::dictionary& unresolvedDict)
{
    if (!unresolvedDict.found("momentumInteraction"))
    {
        return "cell";
    }

    return unresolvedDict.subDict("momentumInteraction")
        .getOrDefault<Foam::word>("fluidVelocity", "cell");
}

inline Foam::scalar erfcScaled(const Foam::scalar z)
{
    if (z < 5.0)
    {
        return Foam::exp(z*z)*Foam::erfc(z);
    }

    Foam::scalar sum = 0.0, term = 1.0;

    for (int k = 1; k <= 12; ++k)
    {
        sum  += term;
        term *= -Foam::scalar(2*k - 1)/(2.0*z*z);
    }

    return sum/(z*Foam::sqrt(Foam::constant::mathematical::pi));
}

inline Foam::scalar psiOseen(const Foam::scalar ReSigma)
{
    const Foam::scalar pi = Foam::constant::mathematical::pi;

    if (ReSigma < 0.2)
    {
        const Foam::scalar r = ReSigma;
        const Foam::scalar k = Foam::sqrt(0.5*pi);

        return 1.0
             - (3.0/8.0)*k*r
             + (1.0/5.0)*r*r
             - (1.0/16.0)*k*r*r*r
             + (1.0/35.0)*r*r*r*r
             - (1.0/128.0)*k*r*r*r*r*r;
    }

    const Foam::scalar z  = ReSigma/Foam::sqrt(2.0);
    const Foam::scalar br = pi
                          - Foam::sqrt(2.0*pi)*ReSigma
                          + 0.5*pi*ReSigma*ReSigma
                          - pi*erfcScaled(z);

    return 3.0/Foam::sqrt(2.0*pi) * br / (ReSigma*ReSigma*ReSigma);
}

inline Foam::scalar selfInducedSigmaHat
(
    const Foam::scalar dp,
    const Foam::scalar sigma
)
{
    return 2.0 * sigma / Foam::max(dp, Foam::SMALL);
}

inline Foam::scalar selfInducedAlphaDeficit(const Foam::scalar sigmaHat)
{
    if (sigmaHat <= 0) return 0.0;

    const Foam::scalar A = selfInducedCentreC / sigmaHat;
    const Foam::scalar E = Foam::exp(-1.0 / (2.0 * sigmaHat * sigmaHat));

    return Foam::erf(1.0 / (sigmaHat * Foam::sqrt(2.0))) - A * E;
}

inline Foam::scalar selfInducedZetaU
(
    const Foam::scalar sampleFactor,
    const Foam::scalar dp,
    const Foam::scalar sigma,
    const Foam::scalar ReP,
    const Foam::scalar phiDrag
)
{
    if (sigma <= 0) return 0.0;

    const Foam::scalar sh = selfInducedSigmaHat(dp, sigma);
    const Foam::scalar A  = selfInducedCentreC / sh;
    const Foam::scalar E  = Foam::exp(-1.0 / (2.0 * sh * sh));

    const Foam::scalar denom = 1.0 - selfInducedAlphaDeficit(sh);

    const Foam::scalar stokesZeta =
        sampleFactor * A * E / Foam::max(denom, Foam::SMALL);

    const Foam::scalar ReSigma = ReP * sigma / Foam::max(dp, Foam::SMALL);

    return stokesZeta * phiDrag * psiOseen(ReSigma);
}

// Gotoh (1990) J. Stat. Phys. 59, 371; Candelier, Mehlig & Magnaudet (2019)
// J. Fluid Mech. 864, 554, eq. (4.10): a sphere at rest in the laboratory, in a
// solid-body rotation of angular velocity Omega, has the in-plane drag
// 6 pi mu a U (1 + c Ta^1/2), Ta = Omega a^2/nu, c = 3 sqrt(2)(19 + 9 sqrt(3))/280,
// derived for Re_p << Ta^1/2 << 1.
inline constexpr Foam::scalar rotationDragC = 0.52409427182290080;

inline Foam::scalar rotationDragFactorStokes(const Foam::scalar Ta)
{
    return 1.0 + rotationDragC*Foam::sqrt(Foam::max(Ta, Foam::scalar(0)));
}

// The rotation screens the far field of the particle's own regularised point
// force as it raises the drag: its mobility at the Stokes radius of the Gaussian,
// a_s = sqrt(pi) sigma (Maxey & Patel 2001), falls by 1/(1 + c Ta_s^1/2),
// Ta_s = Omega a_s^2/nu. The analogue of psiOseen for the inertial screening.
inline Foam::scalar psiRotation
(
    const Foam::scalar omegaRot,
    const Foam::scalar sigma,
    const Foam::scalar nu
)
{
    const Foam::scalar aS2 = Foam::constant::mathematical::pi*sigma*sigma;

    return 1.0/rotationDragFactorStokes(omegaRot*aS2/Foam::max(nu, Foam::SMALL));
}

inline Foam::scalar selfInducedBlobRadius(const Foam::scalar sigma)
{
    return Foam::sqrt(0.5*Foam::constant::mathematical::pi)*sigma;
}

inline Foam::scalar selfInducedWallFactor
(
    const Foam::scalar hOverAEff,
    const Foam::scalar cosSlipWall
)
{
    const Foam::scalar gPar  = 1.0/wallTransForce(hOverAEff);
    const Foam::scalar gPerp = 1.0/wallDragNormal(hOverAEff);

    const Foam::scalar c2 = Foam::min(cosSlipWall*cosSlipWall, Foam::scalar(1));

    return c2*gPerp + (1.0 - c2)*gPar;
}

}

#endif
