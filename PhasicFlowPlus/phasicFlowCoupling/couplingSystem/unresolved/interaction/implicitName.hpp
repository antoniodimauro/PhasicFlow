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

#ifndef __implicitName_hpp__
#define __implicitName_hpp__

#include "OFCompatibleHeader.hpp"

namespace pFlow::coupling
{

inline Foam::word implicitName(const Foam::word& model)
{
    static const Foam::word suffix("Implicit");

    if(model.size() >= suffix.size() &&
       model.compare(model.size() - suffix.size(), suffix.size(), suffix) == 0)
    {
        return model;
    }
    return model + suffix;
}

}

#endif
