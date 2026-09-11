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

#include "noneImplicit.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "distributionBase.hpp"

pFlow::coupling::noneImplicit::noneImplicit
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    liftImplicit(uCS, prsty)
{
    this->setLiftActive(false);

    tmpLiftForce_ = Foam::tmp<Foam::volVectorField>::New
    (
        Foam::IOobject
        (
            "liftForce",
            Foam::timeName(this->mesh().time()),
            this->mesh(),
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE
        ),
        this->mesh(),
        Foam::dimensionedVector
        (
            "liftForce",
            Foam::dimensionSet(1,-2,-2,0,0),
            Foam::vector(0,0,0)
        )
    );
}


void pFlow::coupling::noneImplicit::calculateLiftForceImplicit
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  parVel,
    const Plus::realx3ProcCMField&  parRotVel,
    const Plus::realProcCMField&    diameter,
    const distributionBase&         cellDistribution,
    Plus::realx3ProcCMField&        particleForce,
    Plus::realx3ProcCMField&        particleTorque,
    Plus::realx3ProcCMField&        liftCrossVec
) 
{
    
}