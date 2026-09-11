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

#ifndef __constantCoeffVirtualMassImplicit_hpp__
#define __constantCoeffVirtualMassImplicit_hpp__

#include "virtualMassImplicit.hpp"

namespace pFlow::coupling
{

class distributionBase;

class constantCoeffVirtualMassImplicit 
: 
    public virtualMassImplicit
{
private:

    tmp<Foam::volVectorField> tmpVirtualMassForce_;
    
    Foam::scalar 	Cvm_;
    
    bool isCompressible_ = false;
        
public:

    TypeInfo("constantCoeffImplicit");

    constantCoeffVirtualMassImplicit(const unresolvedCouplingSystem& uCS, const porosity& prsty);

    virtual ~constantCoeffVirtualMassImplicit() = default;

    add_vCtor
    (
        virtualMass,
        constantCoeffVirtualMassImplicit,
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
    
    inline 
    bool isCompressible()const
    {
        return isCompressible_;
    }
    
};

} // pFlow::coupling

#endif
