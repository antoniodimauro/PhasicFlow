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

// Anderson & Jackson (1967) Ind. Eng. Chem. Fundam. 6, 527; Ireland & Desjardins (2017) J. Comput. Phys. 338, 405

#ifndef __particleGaussian_hpp__
#define __particleGaussian_hpp__

#include <vector>

#include "distribution.hpp"
#include "procCMFields.hpp"

namespace pFlow::coupling
{

class couplingMesh;

class particleGaussian
:
    public distribution
{
private:

    Foam::scalar    maxSolidFraction_;

    Foam::scalar    minFilterPerDiameter_;

    Foam::scalar    cFloor_ = 0.0;

    Foam::scalar    minSigmaCells_;

    Foam::scalar    macroLength_;

    mutable bool    warnedMacro_ = false;

    Foam::scalar    nSigma_;

    Foam::scalar    shippingRadius_ = 0.0;

    Foam::label     maxLayers_;

    Plus::procCMField<Foam::scalar>  sigma_;

    mutable std::vector<std::vector<Foam::label>>  stamp_;

    mutable std::vector<Foam::label>               stampCount_;

    mutable std::vector<std::vector<Foam::label>>  frontier_;
    mutable std::vector<std::vector<Foam::label>>  nextFrontier_;

    mutable Foam::label  layersUsed_ = 0;

    bool            truncationAbort_ = true;

    mutable std::vector<char>  procBoundaryCell_;

    mutable Foam::scalar  reportedClearance_ = Foam::GREAT;

    mutable Foam::scalar  maxCellSize_ = -1.0;

    mutable Foam::scalarField cellExtent_;

    mutable Foam::scalar  maxSigmaSeen_    = 0.0;
    mutable Foam::scalar  maxDiameterSeen_ = 0.0;

    void prepareScratch()const;

    void markProcessorBoundaryCells()const;

    void checkProcessorTruncation(Foam::scalar procClearance)const;

    static Foam::scalar gaussianTailBeyond(Foam::scalar rOverSigma);

public:

    TypeInfo("particleGaussian");

    particleGaussian(
        const Foam::dictionary&      parrentDict,
        const couplingMesh&          cMesh,
        const Plus::centerMassField& centerMass);

    ~particleGaussian() = default;

    add_vCtor
    (
        distributionBase,
        particleGaussian,
        dictionary
    );

    inline Foam::scalar maxSigmaBound(Foam::scalar dp)const
    {
        return sigma(dp, maxCellSize());
    }

    Foam::scalar maxCellSize()const;

    const Foam::scalarField& cellExtent()const;

    inline Foam::scalar sigma(Foam::scalar dp, Foam::scalar dCell)const
    {
        return Foam::max(cFloor_*dp, minSigmaCells_*dCell);
    }

    Foam::scalar kernelWidth(Foam::label parIndx)const
    {
        if(parIndx < 0 || static_cast<size_t>(parIndx) >= sigma_.size())
        {
            return 0.0;
        }
        return sigma_[parIndx];
    }

    void updateWeights(const Plus::procCMField<real>& parDiameter) override;

    void setShippingRadius(Foam::scalar r) { shippingRadius_ = r; }

    inline Foam::scalar shippingRadius()const { return shippingRadius_; }

    void localWeightSums(Plus::realProcCMField& sums)const;

    void scaleWeights(const Plus::realProcCMField& totals);

    Foam::word distributionMethodName()const override
    {
        return "particleGaussian";
    }
};

}

#endif
