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

// Shi & Rzehak (2019) Chem. Eng. Sci. 208, 115145, eq. (20)

#ifndef __shearLiftJ_hpp__
#define __shearLiftJ_hpp__

#include "OFCompatibleHeader.hpp"

namespace pFlow::coupling
{

inline Foam::scalar shearLiftJ(const Foam::scalar eps)
{
    if (eps <= 0) return 0.0;

    return 2.255/Foam::pow(1.0 + 0.2/(eps*eps), 1.5);
}

}

#endif
