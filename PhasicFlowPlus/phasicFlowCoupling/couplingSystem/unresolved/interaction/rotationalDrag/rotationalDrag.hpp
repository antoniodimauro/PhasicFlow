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

// Dennis et al. (1980) J. Fluid Mech. 101, 257

#ifndef __rotationalDrag_hpp__
#define __rotationalDrag_hpp__

#include "OFCompatibleHeader.hpp"

namespace pFlow::coupling
{

inline Foam::scalar rotationalDragFactor(const Foam::scalar ReOmega)
{
    const Foam::scalar R = 0.25*ReOmega;

    auto highR = [](const Foam::scalar r)
    {
        return (6.45*Foam::sqrt(r) + 32.1)
             / (16.0*Foam::constant::mathematical::pi);
    };

    auto lowR = [](const Foam::scalar r)
    {
        const Foam::scalar u2 = Foam::sqr(0.1*r);
        Foam::scalar       p  = u2;
        Foam::scalar       f  = 1.0 + (1.0/12.0)*p;
        p *= u2; f -= 0.007542671*p;
        p *= u2; f -= 0.005353489*p;
        p *= u2; f += 0.005824484*p;
        p *= u2; f -= 0.003398210*p;
        p *= u2; f += 0.000896000*p;
        p *= u2; f += 0.000740000*p;
        return f;
    };

    if(R <= 10.0) return lowR(R);
    if(R >= 20.0) return highR(R);

    const Foam::scalar w = 0.1*(R - 10.0);
    return (1.0 - w)*lowR(10.0) + w*highR(20.0);
}

}

#endif
