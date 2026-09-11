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

#include "momentumInteractionImplicit.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "distributionBase.hpp"

namespace pFlow::coupling
{


momentumInteractionImplicit::momentumInteractionImplicit
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    porosity_(prsty),
    momentumInteractionTimer_
    (
        "momentumInteraction", 
        &uCS.couplingTimers()
    )
{

    auto momExch = dict().get<Foam::word>("momentumExchange");

    if(momExch == "distribution")
    {
        momentumExchangeDistribute_ = true;
    }
    else if(momExch == "cell")
    {
        momentumExchangeDistribute_ = false;
    }
    else
    {
        Foam::Info
            << "Unknown momentum exchange method: " << momExch 
            << " in "<< dict().name() << Foam::endl;
        Plus::processor::abort(0);
    }

    drag_ = dragImplicit::create(uCS, porosity_);

    lift_ = liftImplicit::create(uCS, porosity_);

    lift_->setDrag(drag_.get());

    virtualMass_ = virtualMassImplicit::create(uCS, porosity_);

    auto fluidAveragingType = dict().template get<Foam::word>("fluidVelocity");

    fluidAveraging_ = fluidAveraging::create(fluidAveragingType, uCS, "fluidVelocity");

    auto solidAveragingType = dict().template get<Foam::word>("solidVelocity");
    solidAveraging_ = solidAveraging::create(solidAveragingType, uCS, porosity_, "Us");

    requireCellDistribution_ = momentumExchangeDistribute_ ||
        fluidAveraging_ -> requireCellDistribution() ||
        solidAveraging_ -> requireCellDistribution() ||
        drag_ -> requireCellDistribution() ||
        lift_ -> requireCellDistribution() ||
        virtualMass_ -> requireCellDistribution();

    if(requireCellDistribution_)
    {
        INFORMATION<<"Cell distribution is active."<<END_INFO;
    }
    
    if(!momentumExchangeDistribute_)
    {
        noDistribution_ = makeUnique<PCM>(uCS.cMesh(), uCS.centerMass());
    }

    buoyancy_ = dict().getOrDefault<Foam::Switch>("buoyancy", Foam::Switch(false));

    if(buoyancy_)
    {
        Foam::Info
            << "    Buoyancy on particles is "<< Green_Text("active")
            << ": DEM body force is (m_p - rho_f*V_p)*g,\n"
               "      with g taken from settings/settingsDict.\n\n";
    }
}

void momentumInteractionImplicit::calculateDisplacedMass
(
    const Plus::realProcCMField& diameter,
    Plus::realProcCMField&       displacedMass
)const
{
    if(!buoyancy_) return;

    const auto& rho =
        porosity_.mesh().template lookupObject<Foam::volScalarField>("rho");
    const auto& parCellInd = porosity_.parCellIndex();

    const size_t nPar = diameter.size();

    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<nPar; ++i)
    {
        const Foam::label cellI = parCellInd[i];

        if(cellI < 0) continue;

        const Foam::scalar dp = diameter[i];
        const Foam::scalar Vp =
            Foam::constant::mathematical::pi/6.0 * Foam::pow(dp,3.0);

        displacedMass[i] += static_cast<real>(rho[cellI]*Vp);
    }
}

const unresolvedCouplingSystem& momentumInteractionImplicit::uCS()const
{
    return porosity_.uCS();
}


const Foam::dictionary& momentumInteractionImplicit::dict()const
{
    return  momentumInteractionImplicit::getDict(uCS());

}

void momentumInteractionImplicit::calculateCoupling
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  vp,
    const Plus::realx3ProcCMField&  wp,
    Plus::realx3ProcCMField&        fluidForce,
    Plus::realx3ProcCMField&        fluidTorque,
    Plus::realProcCMField&          dragCoeff,
    Plus::realx3ProcCMField&        dragFluidVel,
    Plus::realx3ProcCMField&        liftCrossVec,
    Plus::realProcCMField&          rotDragCoeff,
    Plus::realProcCMField&          addedMass,
    Plus::realProcCMField&          addedInertia,
    Plus::realProcCMField&          displacedMass
)
{

    momentumInteractionTimer_.start();

    // calculate fluid averaging
    fluidAveraging_->calculate(U);

    // calculate solid averaging
    solidAveraging_->calculate(vp);

    const distributionBase& cellDistribution =
        momentumExchangeDistribute_ ? uCS().distribution() : noDistribution_();

    drag_->clearMaterialDerivative();

    // calculate drag force
    drag_->calculateDragForceImplicit
    (
        *fluidAveraging_,
        *solidAveraging_,
        wp,
        porosity_.particleDiameter(),
        cellDistribution,
        fluidForce,
        dragCoeff,
        dragFluidVel);

    if(drag_->bassetActive() && !drag_->bassetEvaluated())
    {
        if(Plus::processor::isMaster())
        {
            fatalErrorInFunction
                << "bassetHistory yes is set, but the selected drag model does "
                   "not evaluate the history force.\n"
                   "It is implemented for the sphereDrag family only."
                << endl;
        }
        Plus::processor::abort(0);
    }

    drag_->addBassetAddedMass(addedMass);

    lift_->calculateLiftForceTorqueImplicit(
        U,
        vp,
        wp,
        porosity_.particleDiameter(),
        cellDistribution,
        fluidForce,
        fluidTorque,
        liftCrossVec,
        rotDragCoeff,
        addedInertia);

    if(lift_->liftReaction())
    {
        drag_->setLiftFeedback(lift_->liftPerParticle(), uCS().centerMass());
    }

    drag_->commitPairForces(lift_->liftPerParticle());

    if(lift_->rotationalHistoryActive() && !lift_->rotationalHistoryEvaluated())
    {
        if(Plus::processor::isMaster())
        {
            fatalErrorInFunction
                << "rotationalHistory yes is set, but the selected "
                   "surfaceRotationTorque model does not drive it."
                << endl;
        }
        Plus::processor::abort(0);
    }

    if(lift_->wallLiftActive() && !lift_->wallLiftApplied())
    {
        if(Plus::processor::isMaster())
        {
            fatalErrorInFunction
                << "wallCorrection yes is set, but the selected lift model "
                   "does not apply the wall-induced lift of Zeng et al. "
                   "(2009).\nIt is implemented for the Shi2019 lift model."
                << endl;
        }
        Plus::processor::abort(0);
    }

    virtualMass_->calculateVirtualMassForceImplicit(
        U,
        drag_->materialDerivativeU(U),
        drag_->faxen() ? &drag_->quadrature() : nullptr,
        porosity_.particleDiameter(),
        cellDistribution,
        fluidForce,
        addedMass);

    calculateDisplacedMass(porosity_.particleDiameter(), displacedMass);

    momentumInteractionTimer_.end();

    Foam::Info<<Blue_Text("Momentum interaction time: ")<< 
                Yellow_Text(momentumInteractionTimer_.lastTime())<<
                Yellow_Text(" s")<<Foam::endl;
}


const Foam::dictionary& momentumInteractionImplicit::getDict(const unresolvedCouplingSystem& uCS)
{
    return uCS.unresolvedDict().subDict("momentumInteraction");	
}


} // pFlow::coupling
