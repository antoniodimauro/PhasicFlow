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

// Maxey & Riley (1983) Phys. Fluids 26, 883
#include "noneVirtualMassImplicit.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "distributionBase.hpp"

pFlow::coupling::noneVirtualMassImplicit::noneVirtualMassImplicit
(
    const unresolvedCouplingSystem& uCS,
    const porosity& prsty
)
:
    virtualMassImplicit(uCS, prsty)
{
    tmpVirtualMassForce_ = Foam::tmp<Foam::volVectorField>::New
    (
        Foam::IOobject
        (
            "virtualMassForce",
            Foam::timeName(this->mesh().time()),
            this->mesh(),
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE
        ),
        this->mesh(),
        Foam::dimensionedVector
        (
            "virtualMassForce",
            Foam::dimensionSet(1,-2,-2,0,0),
            Foam::vector(0,0,0)
        )
    );
}

void pFlow::coupling::noneVirtualMassImplicit::calculateVirtualMassForceImplicit
(
    const Foam::volVectorField& U,
    const Foam::volVectorField& DDtUField,
    sphereQuadrature* quadPtr,
    const Plus::realProcCMField& diameter,
    const distributionBase& cellDistribution,
    Plus::realx3ProcCMField& particleForce,
    Plus::realProcCMField& addedMass
)
{
}

