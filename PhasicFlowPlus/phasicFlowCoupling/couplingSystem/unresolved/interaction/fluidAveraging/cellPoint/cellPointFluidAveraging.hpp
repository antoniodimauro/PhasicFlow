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

#ifndef __cellPointFluidAveraging_hpp__
#define __cellPointFluidAveraging_hpp__

#include "fluidAveraging.hpp"

namespace pFlow::coupling
{

class cellPointFluidAveraging
:
    public fluidAveraging
{
public:

    TypeInfo("cellPoint");

    cellPointFluidAveraging(
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const word&                     name);

    virtual ~cellPointFluidAveraging() = default;

    add_vCtor
    (
        fluidAveraging,
        cellPointFluidAveraging,
        word
    );

    void calculate(const Foam::volVectorField& orgField);

    bool requireCellDistribution()const override
    {
        return false;
    }
};

}

#endif //__cellPointFluidAveraging_hpp__
