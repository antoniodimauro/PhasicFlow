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

#include "virtualMassImplicit.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "implicitName.hpp"

namespace pFlow::coupling
{

virtualMassImplicit::virtualMassImplicit
(
    const unresolvedCouplingSystem& uCS,
    const porosity& prsty
)
:
    virtualMass(uCS, prsty),
    ddtUSampler_(uCS, "DDtU"),
    vmScratch_("virtualMassScratch", realx3(0), uCS.centerMass())
{}

void virtualMassImplicit::reportAddedMassFaxen(Foam::scalar ratioObserved)const
{
    Foam::reduce(ratioObserved, Foam::maxOp<Foam::scalar>());

    if(ratioObserved > reportedAmFaxen_ + 0.05)
    {
        reportedAmFaxen_ = ratioObserved;

        Foam::Info
            << Blue_Text("Added-mass Faxen: ")
            << "volume mean departs from the point value by "
            << Yellow_Text(100.0*ratioObserved)
            << Yellow_Text(" % at a particle")
            << Foam::endl;
    }
}

void virtualMassImplicit::calculateVirtualMassForce
(
    const Foam::volVectorField& U,
    const Plus::realx3ProcCMField& parAcc,
    const Plus::realProcCMField& diameter,
    Plus::realx3ProcCMField& particleForce
)
{
    fatalErrorInFunction
        << "the explicit virtual mass interface is not available for "
        << typeName() << ", use calculateVirtualMassForceImplicit" << endl;
    fatalExit;
}

uniquePtr<virtualMassImplicit> virtualMassImplicit::create
(
    const unresolvedCouplingSystem& uCS,
    const porosity& prsty
)
{
    const auto& vmDict = virtualMass::getDict(uCS);

    auto vmType = implicitName(vmDict.getOrDefault<Foam::word>("model", "none"));

    if(virtualMass::couplingSystemvCtorSelector_.search(vmType))
    {
        Foam::Info<<"    Creating virtual mass force "<<Green_Text(vmType)<<" ...\n\n";
        uniquePtr<virtualMass> basePtr =
            virtualMass::couplingSystemvCtorSelector_[vmType](uCS, prsty);
        auto* implicitPtr = dynamic_cast<virtualMassImplicit*>(basePtr.release());
        if(!implicitPtr)
        {
            fatalErrorInFunction
                << vmType << " is not an implicit virtual mass model" << endl;
            fatalExit;
        }
        return uniquePtr<virtualMassImplicit>(implicitPtr);
    }

    if(Plus::processor::isMaster())
    {
        printKeys
        (
            fatalErrorInFunction << "Ctor Selector "<< vmType << " does not exist"
            " for virtualMass method in "<< vmDict.name()
            <<"\nAvailable ones are: \n",
            virtualMass::couplingSystemvCtorSelector_
        )<<endl;
    }
    Plus::processor::abort(0);

    return nullptr;
}

} // pFlow::coupling
