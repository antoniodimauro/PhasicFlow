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

#include "liftImplicit.hpp"
#include "dragImplicit.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "implicitName.hpp"
#include "distributionBase.hpp"

pFlow::coupling::liftImplicit::liftImplicit
(
    const unresolvedCouplingSystem& uCS,
    const porosity& 				prsty
)
:
    lift(uCS, prsty),
    liftUSampler_(uCS, "liftU"),
    curlUSampler_(uCS, "curlU"),
    liftScratch_("liftForceScratch", realx3(0), uCS.centerMass()),
    surfTorqueImplicit_(surfaceRotationTorqueImplicit::create(uCS, prsty))
{
    const Foam::dictionary& dragDict =
        uCS.unresolvedDict().subDict("momentumInteraction").subDict("drag");

    wallCorr_ = dragDict.getOrDefault<Foam::Switch>("wallCorrection", false);

    liftReaction_ = uCS.unresolvedDict().subDict("momentumInteraction").subDict("lift")
        .getOrDefault<Foam::Switch>("liftReaction", true);

    if(!liftReaction_)
    {
        Foam::Info
            << "    Lift reaction: " << Green_Text("not returned to the fluid")
            << ". The lift closures give the total lift, the\n"
               "      particle's own disturbance included; returned to the fluid,"
               " the lift drives a flow\n      that the particle would sample"
               " again. No self-induced lift feedback either.\n";
    }

    if(wallCorr_)
    {
        minWallRadii_ = dragDict.getOrDefault<Foam::scalar>("minWallRadii", 2.0);

        Foam::Info
            << "    Wall-induced lift: " << Green_Text("on")
            << " (Zeng et al., Phys. Fluids 21 (2009) 033302,\n"
               "      eqs 30-35, the translation-induced half. Their"
               " shear-induced half is a\n      total that already contains the"
               " unbounded lift, so it is not added.)\n";
    }
}

void pFlow::coupling::liftImplicit::setDrag(const dragImplicit* d)
{
    drag_ = d;
    if(d && surfTorqueImplicit_)
    {
        surfTorqueImplicit_->setPair(d->pairModelPtr());
    }
}

void pFlow::coupling::liftImplicit::reportWallLift(Foam::scalar fractionObserved)const
{
    if(!wallCorr_) return;

    wallLiftApplied_ = true;

    Foam::reduce(fractionObserved, Foam::maxOp<Foam::scalar>());

    if(!notedShearLift_ && Plus::processor::isMaster())
    {
        notedShearLift_ = true;

        if(Plus::processor::isMaster())
        {
            INFORMATION
                << "Wall lift: the TRANSLATION-induced half of Zeng et al. (2009) "
                << "is applied. The SHEAR-induced half, their eq (28), is not: it "
                << "is a total that already contains the unbounded shear lift this "
                << "run gets from its lift model, and their parameterisation ties "
                << "the shear rate to the wall distance through one Reynolds "
                << "number, so the wall part cannot be separated out. The wall's "
                << "effect on the shear-induced lift is therefore NOT modelled."
                << END_INFO;
        }
    }

    if(fractionObserved > reportedWallLift_ + 0.05)
    {
        reportedWallLift_ = fractionObserved;

        Foam::Info
            << Blue_Text("Wall lift: ")
            << "largest wall-normal lift is "
            << Yellow_Text(100.0*fractionObserved)
            << Yellow_Text(" %") << " of the drag on that particle"
            << Foam::endl;
    }
}

void pFlow::coupling::liftImplicit::calculateLiftForceTorqueImplicit
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  parVel,
    const Plus::realx3ProcCMField&  parRotVel,
    const Plus::realProcCMField&    diameter,
    const distributionBase&         cellDistribution,
    Plus::realx3ProcCMField&        particleForce,
    Plus::realx3ProcCMField&        particleTorque,
    Plus::realx3ProcCMField&        liftCrossVec,
    Plus::realProcCMField&          rotDragCoeff,
    Plus::realProcCMField&          addedInertia
)
{
    calculateLiftForceImplicit(U, parVel, parRotVel, diameter, cellDistribution,
        particleForce, particleTorque, liftCrossVec);
    surfTorqueImplicit_->calculateTorque(U, parVel, parRotVel, diameter,
        particleTorque, rotDragCoeff, addedInertia);
}

void pFlow::coupling::liftImplicit::calculateLiftForce
(
    const Foam::volVectorField&,
    const Plus::realx3ProcCMField&,
    const Plus::realx3ProcCMField&,
    const Plus::realProcCMField&,
    Plus::realx3ProcCMField&,
    Plus::realx3ProcCMField&
)
{
    fatalErrorInFunction
        << "an implicit lift model was selected outside the momentumSemiImplicit "
           "coupling system."
        << endl;
    Plus::processor::abort(0);
}

pFlow::uniquePtr<pFlow::coupling::liftImplicit> pFlow::coupling::liftImplicit::create
(
    const unresolvedCouplingSystem& uCS,
    const porosity& 				prsty
)
{
    const auto& liftDict = lift::getDict(uCS);
    auto liftType = implicitName(liftDict.getOrDefault<Foam::word>("model", "none"));

    if( lift::couplingSystemvCtorSelector_.search(liftType))
    {
        Foam::Info<<"    Creating lift force "<<Green_Text(liftType)<<" ...\n\n";
        auto basePtr = lift::couplingSystemvCtorSelector_[liftType] (uCS, prsty);
        return uniquePtr<liftImplicit>(dynamic_cast<liftImplicit*>(basePtr.release()));
    }

    if(Plus::processor::isMaster())
    {
        printKeys
        (
            fatalErrorInFunction << "Ctor Selector "<< liftType << " does not exist"
            " for lift method in "<< liftDict.name()
            <<"\nAvailable ones are: \n"
            ,
            lift::couplingSystemvCtorSelector_
        )<<endl;
    }
    Plus::processor::abort(0);
    return nullptr;
}

Foam::vector pFlow::coupling::liftImplicit::fluidVelForLift
(
    size_t              i,
    const Foam::vector& sampled
)const
{
    Foam::vector u;

    if(drag_ && drag_->undisturbedVelocity(i, u)) return u;

    return sampled;
}
