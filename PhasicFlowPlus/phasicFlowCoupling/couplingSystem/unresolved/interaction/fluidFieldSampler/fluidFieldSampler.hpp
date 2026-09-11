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

// Anderson & Jackson (1967) Ind. Eng. Chem. Fundam. 6, 527

#ifndef __fluidFieldSampler_hpp__
#define __fluidFieldSampler_hpp__

#include "OFCompatibleHeader.hpp"
#include "procCMFields.hpp"
#include "uniquePtr.hpp"
#include "fluidAveraging.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class fluidFieldSampler
{
private:

    word                        method_;

    uniquePtr<fluidAveraging>   averaging_ = nullptr;

public:

    fluidFieldSampler(
        const unresolvedCouplingSystem& uCS,
        const word&                     name);

    ~fluidFieldSampler() = default;

    const Plus::procCMField<Foam::vector>&
    sample(const Foam::volVectorField& src);

    bool requireCellDistribution()const;

    const word& method()const
    {
        return method_;
    }

    static word readMethod(const unresolvedCouplingSystem& uCS);
};

}

#endif
