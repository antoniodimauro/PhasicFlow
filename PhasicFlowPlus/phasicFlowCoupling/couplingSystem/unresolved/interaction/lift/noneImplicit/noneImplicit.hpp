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

#ifndef __noneImplicit_hpp__
#define __noneImplicit_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"

// from PhasicFlowPlus
#include "procCMFields.hpp"
#include "liftImplicit.hpp"


namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class noneImplicit
:
    public liftImplicit
{
private:

    tmp<Foam::volVectorField>   tmpLiftForce_;
    
public:

    // type info
    TypeInfo("noneImplicit");

    noneImplicit(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);
    

    virtual ~noneImplicit() = default;
    
    add_vCtor
    (
        lift,
        noneImplicit,
        couplingSystem
    );

    void calculateLiftForceImplicit(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        const distributionBase&         cellDistribution,
        Plus::realx3ProcCMField&        particleForce,
        Plus::realx3ProcCMField&        particleTorque,
        Plus::realx3ProcCMField&        liftCrossVec) override;

    Foam::tmp<Foam::volVectorField> liftForce()const override
    {
        return tmpLiftForce_;
    }

};
    
}

#endif // __noneImplicit_hpp__