/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Original file: PhasicFlowPlus/phasicFlowCoupling/couplingSystem/unresolved/interaction/virtualMass/none/noneVirtualMass.hpp,
  by Hamidreza Norouzi and Bahram Haddadi (PhasicFlowPlus as distributed, commit bbf2110).
  Copyright (C): Antonio Di Mauro
  email: antoniodimauro03@gmail.com
  Derived from the original file in 2026 by Antonio Di Mauro
  (T-junction CFD-DEM study).
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
#ifndef __noneVirtualMassImplicit_hpp__
#define __noneVirtualMassImplicit_hpp__

#include "virtualMassImplicit.hpp"

namespace pFlow::coupling
{

class distributionBase;

class unresolvedCouplingSystem;

class noneVirtualMassImplicit : public virtualMassImplicit
{
private:

    tmp<Foam::volVectorField> tmpVirtualMassForce_;

public:

    TypeInfo("noneImplicit");

    noneVirtualMassImplicit(const unresolvedCouplingSystem& uCS, const porosity& prsty);

    virtual ~noneVirtualMassImplicit() = default;

    add_vCtor
    (
        virtualMass,
        noneVirtualMassImplicit,
        couplingSystem
    );

    void calculateVirtualMassForceImplicit
    (
        const Foam::volVectorField& U,
        const Foam::volVectorField& DDtUField,
        sphereQuadrature* quadPtr,
        const Plus::realProcCMField& diameter,
        const distributionBase& cellDistribution,
        Plus::realx3ProcCMField& particleForce,
        Plus::realProcCMField& addedMass
    ) override;

    Foam::tmp<Foam::volVectorField> virtualMassForce()const override
    {
        return tmpVirtualMassForce_;
    }
};

} // pFlow::coupling


#endif
