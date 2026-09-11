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

// Shi & Rzehak (2019) Chem. Eng. Sci. 208, 115145

#ifndef __Shi2019Implicit_hpp__
#define __Shi2019Implicit_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"

// from PhasicFlowPlus
#include "procCMFields.hpp"
#include "liftImplicit.hpp"


namespace pFlow::coupling
{


class Shi2019Implicit
:
    public liftImplicit
{
private:

    tmp<Foam::volVectorField>   tmpLiftForce_;

    Foam::scalar                residualRe_;

    // rotation lift on the vorticity in excess of the strain rate:
    // none, or Bluemink2010 (Bluemink et al. 2010, eq. 2.8)
    bool                        rotationLift_;

    // rotationLiftEkman: the Ekman layer (nu/Omega)^1/2 thins the boundary
    // layer as the inertial one (nu a/|u - v|)^1/2 does; the fit of Bluemink
    // et al. (2010) is taken at Re_a + Ta, i.e. Re_p + 2 Ta, Ta = Omega a^2/nu
    bool                        rotationLiftEkman_ = false;

    // spinLiftHighRotation OesterleDinh1998: at high Rr = |omega_p| d/|u - v|
    // the spin-lift coefficient of Oesterle & Dinh (1998, Exp. Fluids 25, 16),
    // C = 0.45 + (Rr - 0.45) exp(-0.05684 Rr^0.4 Re_p^0.7), fitted on
    // 10 <= Re_p <= 140, 2 <= Rr <= 12 and equal to Rubinow & Keller (1961) as
    // Re_p -> 0. Shi & Rzehak (2019) below Rr = 2, Oesterle & Dinh above
    // Rr = 5, and between them a smoothstep blend over the range where both
    // hold (Shi & Rzehak propose theirs for 0.1 <= Rr <= 10).
    bool                        spinLiftOesterleDinh_ = false;

    uniquePtr<fluidFieldSampler> strainSampler_;
    
public:

    // type info
    TypeInfo("Shi2019Implicit");

    Shi2019Implicit(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);
    

    virtual ~Shi2019Implicit() = default;
    
    add_vCtor
    (
        lift,
        Shi2019Implicit,
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

#endif // __Shi2019Implicit_hpp__