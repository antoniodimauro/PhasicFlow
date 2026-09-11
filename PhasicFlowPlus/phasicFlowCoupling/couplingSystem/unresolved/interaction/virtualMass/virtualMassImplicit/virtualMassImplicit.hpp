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

#ifndef __virtualMassImplicit_hpp__
#define __virtualMassImplicit_hpp__

#include "virtualMass.hpp"
#include "fluidFieldSampler.hpp"
#include "sphereQuadrature.hpp"

namespace pFlow::coupling
{

class distributionBase;

class unresolvedCouplingSystem;

class virtualMassImplicit
:
    public virtualMass
{
private:

    fluidFieldSampler    ddtUSampler_;

    mutable Plus::realx3ProcCMField vmScratch_;

    mutable Foam::scalar         reportedAmFaxen_ = -1.0;

protected:

    const Plus::procCMField<Foam::vector>&
    sampleDDtU(const Foam::volVectorField& ddtU)
    {
        return ddtUSampler_.sample(ddtU);
    }

    void reportAddedMassFaxen(Foam::scalar ratioObserved)const;

    Plus::realx3ProcCMField& vmScratch()const { return vmScratch_; }

public:

    TypeInfo("virtualMassImplicit");

    virtualMassImplicit(const unresolvedCouplingSystem& uCS, const porosity& prsty);

    virtual ~virtualMassImplicit() = default;

    bool requireCellDistribution()const
    {
        return ddtUSampler_.requireCellDistribution();
    }

    virtual
    void calculateVirtualMassForceImplicit
    (
        const Foam::volVectorField& U,
        const Foam::volVectorField& DDtUField,
        sphereQuadrature* quadPtr,
        const Plus::realProcCMField& diameter,
        const distributionBase& cellDistribution,
        Plus::realx3ProcCMField& particleForce,
        Plus::realProcCMField& addedMass
    ) = 0;

    void calculateVirtualMassForce
    (
        const Foam::volVectorField& U,
        const Plus::realx3ProcCMField& parAcc,
        const Plus::realProcCMField& diameter,
        Plus::realx3ProcCMField& particleForce
    ) override;

    static
    uniquePtr<virtualMassImplicit> create
    (
        const unresolvedCouplingSystem& uCS,
        const porosity& prsty
    );
};

} // pFlow::coupling

#endif
