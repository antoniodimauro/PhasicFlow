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

#include "fluidFieldSampler.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "streams.hpp"
#include "processorPlus.hpp"

pFlow::word
pFlow::coupling::fluidFieldSampler::readMethod
(
    const unresolvedCouplingSystem& uCS
)
{
    const Foam::dictionary& dict =
        uCS.unresolvedDict().subDict("momentumInteraction");

    const Foam::word m =
        dict.getOrDefault<Foam::word>("fluidSampling", "cell");

    if(m != "cell" && m != "interpolate" && m != "cellPoint" && m != "distribution")
    {
        if(Plus::processor::isMaster())
        {
            fatalErrorInFunction
                << "Unknown fluidSampling \"" << m << "\" in "
                << dict.name()
                << "\nAvailable ones are: cell, interpolate, cellPoint, distribution"
                << endl;
        }
        Plus::processor::abort(0);
    }

    return m;
}

pFlow::coupling::fluidFieldSampler::fluidFieldSampler
(
    const unresolvedCouplingSystem& uCS,
    const word&                     name
)
:
    method_(readMethod(uCS)),
    averaging_(fluidAveraging::create(method_, uCS, name))
{}

const pFlow::Plus::procCMField<Foam::vector>&
pFlow::coupling::fluidFieldSampler::sample
(
    const Foam::volVectorField& src
)
{
    averaging_().calculate(src);
    return averaging_().field();
}

bool pFlow::coupling::fluidFieldSampler::requireCellDistribution()const
{
    return averaging_().requireCellDistribution();
}
