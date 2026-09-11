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

#include "cellPointFluidAveraging.hpp"
#include "unresolvedCouplingSystem.hpp"

#include "interpolationCellPoint.H"

pFlow::coupling::cellPointFluidAveraging::cellPointFluidAveraging
(
    const word&                     type,
    const unresolvedCouplingSystem& uCS,
    const word&                     name
)
:
    fluidAveraging(type, uCS, name)
{
}

void pFlow::coupling::cellPointFluidAveraging::calculate(const Foam::volVectorField& orgField)
{
    const size_t numPar = averagedField_.size();
    const Plus::centerMassField& centerMass = averagedField_.centerMass();
    const Foam::fvMesh& mesh = uCS_.cMesh().mesh();
    const Plus::procCMField<Foam::label>& parCellIndex = uCS_.parCellIndex();

    mesh.cellCentres();
    mesh.faceCentres();
    mesh.tetBasePtIs();

    const Foam::interpolationCellPoint<Foam::vector> interp(orgField);

    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<numPar; i++)
    {
        const auto celli = parCellIndex[i];

        if(celli < 0)
        {
            averagedField_[i] = Foam::vector(0,0,0);
            continue;
        }

        const Foam::vector p{centerMass[i].x(), centerMass[i].y(), centerMass[i].z()};

        averagedField_[i] = interp.interpolate(p, celli);
    }
}
