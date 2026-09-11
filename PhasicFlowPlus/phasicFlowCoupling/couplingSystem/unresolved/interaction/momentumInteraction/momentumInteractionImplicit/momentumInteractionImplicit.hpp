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

#ifndef __momentumInteractionImplicit_hpp__
#define __momentumInteractionImplicit_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "Timer.hpp"

// from phasicFlowPlus
#include "procCMFields.hpp"
#include "dragImplicit.hpp"
#include "liftImplicit.hpp"
#include "virtualMassImplicit.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"
#include "PCM.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;
class porosity;


class momentumInteractionImplicit
{
private:

    const porosity&                 porosity_;

    bool                            momentumExchangeDistribute_;

    bool                            requireCellDistribution_ = false;

    uniquePtr<dragImplicit>         drag_;

    uniquePtr<liftImplicit>         lift_;

    uniquePtr<virtualMassImplicit>  virtualMass_;

    uniquePtr<fluidAveraging>       fluidAveraging_;

    uniquePtr<solidAveraging>       solidAveraging_;

    uniquePtr<PCM>                  noDistribution_ = nullptr;

    Timer                           momentumInteractionTimer_;

    bool                            buoyancy_ = false;

    void calculateDisplacedMass(
        const Plus::realProcCMField& diameter,
        Plus::realProcCMField&       displacedMass)const;

public:

    TypeInfo("momentumInteractionImplicit");

    momentumInteractionImplicit(
        const unresolvedCouplingSystem& uCS,
        const porosity& prsty);

    virtual ~momentumInteractionImplicit() = default;

    inline
    const Foam::volVectorField& Su()const
    {
        return std::as_const<const drag&>(*drag_).Su();
    }

    inline
    const Foam::volScalarField& Sp()const
    {
        return std::as_const<const drag&>(*drag_).Sp();
    }

    inline
    Foam::tmp<Foam::volVectorField> liftForce()const
    {
        return lift_->liftForce();
    }

    inline
    Foam::tmp<Foam::volVectorField> virtualMassForce()const
    {
        return virtualMass_->virtualMassForce();
    }

    inline
    const porosity& Porosity()const
    {
        return porosity_;
    }

    const unresolvedCouplingSystem& uCS()const;

    inline
    bool requireCellDistribution()const
    {
        return requireCellDistribution_;
    }

    void calculateCoupling(
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
        Plus::realProcCMField&          displacedMass);

    inline
    bool buoyancy()const
    {
        return buoyancy_;
    }

    void commitHistoryForces()const
    {
        drag_->commitBassetHistory();
        lift_->commitRotationalHistory();
    }

    const Foam::dictionary& dict()const;

    static
    const Foam::dictionary& getDict(const unresolvedCouplingSystem& uCS);
};

}

#endif //__momentumInteractionImplicit_hpp__
