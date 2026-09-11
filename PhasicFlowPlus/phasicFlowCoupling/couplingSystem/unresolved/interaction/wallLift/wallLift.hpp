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

// Zeng et al. (2009) Phys. Fluids 21, 033302

#ifndef __wallLift_hpp__
#define __wallLift_hpp__

#include "OFCompatibleHeader.hpp"

namespace pFlow::coupling
{

inline Foam::scalar zengLiftLowRe(const Foam::scalar Lstar)
{
    if (Lstar <= 0) return 1.125;

    if (Lstar < 10.0)
    {
        return (1.125 + 5.78e-6*Foam::pow(Lstar, 4.58))
              *Foam::exp(-0.292*Lstar);
    }

    return 8.94*Foam::pow(Lstar, -2.09);
}

inline Foam::scalar zengLiftContact(const Foam::scalar Ret)
{
    return 0.313 + 0.812*Foam::exp(-0.125*Foam::pow(Ret, 0.77));
}

inline Foam::scalar zengLiftOuter
(
    const Foam::scalar L,
    const Foam::scalar Ret
)
{
    const Foam::scalar A = 1.0 + 0.329*Ret + 0.00485*Ret*Ret;
    const Foam::scalar g = -0.9*Foam::tanh(0.022*Ret);

    return A*zengLiftLowRe(L*Ret)*Foam::pow(L, g);
}

inline Foam::scalar zengWallLift(const Foam::scalar L, const Foam::scalar Ret)
{
    const Foam::scalar delta  = Foam::max(L - 0.5, Foam::scalar(0));
    const Foam::scalar deltaC = 3.0*Foam::exp(-0.17*Foam::pow(Ret, 0.7));

    const Foam::scalar near =
        zengLiftContact(Ret) - zengLiftOuter(0.5, Ret);

    return zengLiftOuter(L, Ret)
         + near*Foam::exp(-11.0*Foam::pow(delta/deltaC, 1.2));
}

}

#endif
