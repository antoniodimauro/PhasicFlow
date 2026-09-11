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

#include "surfaceRotationTorqueImplicit.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "implicitName.hpp"

pFlow::coupling::surfaceRotationTorqueImplicit::surfaceRotationTorqueImplicit
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    surfaceRotationTorque(uCS, prsty),
    torqueUSampler_(uCS, "torqueU"),
    torqueCurlUSampler_(uCS, "torqueCurlU"),
    history_(uCS, getDict(uCS))
{
    const Foam::dictionary& dragDict =
        uCS.unresolvedDict().subDict("momentumInteraction").subDict("drag");

    wallCorr_ = dragDict.getOrDefault<Foam::Switch>("wallCorrection", false);

    if(wallCorr_)
    {
        minWallRadii_ = dragDict.getOrDefault<Foam::scalar>("minWallRadii", 2.0);

        Foam::Info
            << "    Near-wall couple: " << Green_Text("on")
            << " (Goldman, Cox & Brenner 1967 I table 2 and II\n"
               "      table 1: the wall resists spin, weakens the shear couple,"
               " and makes\n      a translating sphere feel one at all)\n";
    }
}

void pFlow::coupling::surfaceRotationTorqueImplicit::calculateTorque
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  parVel,
    const Plus::realx3ProcCMField&  parRotVel,
    const Plus::realProcCMField&    diameter,
    Plus::realx3ProcCMField&        particleTorque,
    Plus::realProcCMField&          rotDragCoeff,
    Plus::realProcCMField&          addedInertia
)
{
    const size_t nPar = diameter.size();

    history_.beginCoupling(nPar);

    calculateSurfaceTorqueImplicit(U, parVel, parRotVel, diameter,
                           particleTorque, rotDragCoeff);

    if(!history_.active()) return;

    const Foam::Time& rt = Porosity().uCS().cMesh().mesh().time();

    history_.evaluate(rt.value(), rt.deltaTValue());

    const auto& parCellInd = this->parCellIndex();

    for(size_t i=0; i<nPar; ++i)
    {
        if(parCellInd[i] < 0) continue;

        const Foam::vector t = history_.torque(i);

        particleTorque[i] += realx3(t.x(), t.y(), t.z());
        addedInertia[i]   += static_cast<real>(history_.addedInertia(i));
    }
}

void pFlow::coupling::surfaceRotationTorqueImplicit::calculateSurfaceTorque
(
    const Foam::volVectorField&,
    const Plus::realx3ProcCMField&,
    const Plus::realx3ProcCMField&,
    const Plus::realProcCMField&,
    Plus::realx3ProcCMField&
)
{
    fatalErrorInFunction
        << "an implicit surfaceRotationTorque model was selected outside the "
           "momentumSemiImplicit coupling system."
        << endl;
    Plus::processor::abort(0);
}

pFlow::uniquePtr<pFlow::coupling::surfaceRotationTorqueImplicit>
    pFlow::coupling::surfaceRotationTorqueImplicit::create
(
    const unresolvedCouplingSystem &uCS,
    const porosity &prsty
)
{
    const auto& liftDict = surfaceRotationTorque::getDict(uCS);
    auto type = implicitName(liftDict.getOrDefault<Foam::word>("surfaceRotationTorque", "none"));

    if( surfaceRotationTorque::couplingSystemvCtorSelector_.search(type))
    {
        Foam::Info<<"    Creating surfaceRotationTorque model "<<Green_Text(type)<<" ...\n\n";
        auto basePtr = surfaceRotationTorque::couplingSystemvCtorSelector_[type] (uCS, prsty);
        return uniquePtr<surfaceRotationTorqueImplicit>(
            dynamic_cast<surfaceRotationTorqueImplicit*>(basePtr.release()));
    }

    if(Plus::processor::isMaster())
    {
        printKeys
        (
            fatalErrorInFunction << "Ctor Selector "<< type << " does not exist"
            " for surfaceRotationTorque method in "<< liftDict.name()
            <<"\nAvailable ones are: \n"
            ,
            surfaceRotationTorque::couplingSystemvCtorSelector_
        )<<endl;
    }
    Plus::processor::abort(0);
    return nullptr;
}
