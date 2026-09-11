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

#include <vector>
#include <cmath>

#include "particleGaussian.hpp"
#include "couplingMesh.hpp"
#include "streams.hpp"
#include "processorPlus.hpp"
#include "selfInduction.hpp"

#include "processorPolyPatch.H"

#include <sstream>

#ifdef _OPENMP
    #include <omp.h>
#endif

pFlow::coupling::particleGaussian::particleGaussian
(
    const Foam::dictionary&      parrentDict,
    const couplingMesh&          cMesh,
    const Plus::centerMassField& centerMass
)
:
    distribution(parrentDict, cMesh, centerMass),
    maxSolidFraction_(lookupOrDefaultDict<Foam::scalar>(
        infoDict(parrentDict, "particleGaussianInfo"), "maxSolidFraction", 0.05)),
    minFilterPerDiameter_(lookupOrDefaultDict<Foam::scalar>(
        infoDict(parrentDict, "particleGaussianInfo"), "minFilterPerDiameter", 4.0)),
    minSigmaCells_(lookupOrDefaultDict<Foam::scalar>(
        infoDict(parrentDict, "particleGaussianInfo"), "minSigmaCells", 1.6986)),
    macroLength_(lookupOrDefaultDict<Foam::scalar>(
        infoDict(parrentDict, "particleGaussianInfo"), "macroscopicLength", 0.0)),
    nSigma_(lookupOrDefaultDict<Foam::scalar>(
        infoDict(parrentDict, "particleGaussianInfo"), "nSigma", 3.5)),
    maxLayers_(lookupOrDefaultDict<Foam::label>(
        infoDict(parrentDict, "particleGaussianInfo"), "maxLayers", 60)),
    sigma_("kernelSigma", Foam::scalar(0), centerMass)
{
    const Foam::word truncation = lookupOrDefaultDict<Foam::word>(
        infoDict(parrentDict, "particleGaussianInfo"),
        "processorTruncation",
        "abort");

    if(truncation == "abort")      truncationAbort_ = true;
    else if(truncation == "warn")  truncationAbort_ = false;
    else
    {
        if(Plus::processor::isMaster())
        {
            fatalErrorInFunction
                << "Unknown processorTruncation \"" << truncation << "\".\n"
                << "Available ones are: abort, warn" << endl;
        }
        Plus::processor::abort(0);
    }

    const Foam::scalar pi = Foam::constant::mathematical::pi;

    const Foam::scalar cPhi = Foam::pow
    (
        (pi/6.0) / (Foam::pow(2.0*pi, 1.5)*Foam::max(maxSolidFraction_, SMALL)),
        1.0/3.0
    );

    const Foam::scalar cBeta = minFilterPerDiameter_
                             / (2.0*Foam::sqrt(2.0*Foam::log(2.0)));

    cFloor_ = Foam::max(cPhi, cBeta);

    Foam::Info
        << "    particleGaussian: sigma_i = max( "
        << Green_Text(cFloor_) << "*d_p,  "
        << Green_Text(minSigmaCells_) << "*dx )\n"
        << "      both floors are requirements of the volume-filtered\n"
        << "      formulation (Anderson & Jackson 1967; Ireland & Desjardins,\n"
        << "      JCP 338 (2017) 405), not tuning:\n"
        << "        minFilterPerDiameter " << minFilterPerDiameter_
        << " -> sigma >= " << cBeta << "*d_p   delta_f/d_p (I&D validate 2-16)\n"
        << "        maxSolidFraction  " << maxSolidFraction_
        << "  -> sigma >= " << cPhi  << "*d_p   peak deposited solid fraction\n"
        << "        minSigmaCells     " << minSigmaCells_
        << "  -> sigma >= " << minSigmaCells_
        << "*dx     I&D run delta_f/dx = 4\n"
        << "      binding term: the particle floor above d/dx = "
        << minSigmaCells_/cFloor_ << ", the mesh floor below it\n"
        << (macroLength_ > 0
            ? "        macroscopicLength given -> support is checked against it\n"
            : "        macroscopicLength NOT set -> the Anderson & Jackson upper\n"
              "          bound (filter small vs the flow scale) is UNCHECKED\n")
        << "      delta_f/d_p = " << 2.0*Foam::sqrt(2.0*Foam::log(2.0))*cFloor_
        << " at the particle floor (I&D validate 2 to 16)\n";
}

void pFlow::coupling::particleGaussian::prepareScratch()const
{
    #ifdef _OPENMP
        const size_t nThreads = static_cast<size_t>(omp_get_max_threads());
    #else
        const size_t nThreads = 1;
    #endif

    const size_t nCells = static_cast<size_t>(mesh().nCells());

    if(stamp_.size() != nThreads)
    {
        stamp_.assign(nThreads, std::vector<Foam::label>());
        stampCount_.assign(nThreads, 0);
        frontier_.assign(nThreads, std::vector<Foam::label>());
        nextFrontier_.assign(nThreads, std::vector<Foam::label>());
    }

    for(size_t t = 0; t < nThreads; ++t)
    {
        if(stamp_[t].size() != nCells || stampCount_[t] >= Foam::labelMax - 1)
        {
            stamp_[t].assign(nCells, -1);
            stampCount_[t] = 0;
        }
    }
}

const Foam::scalarField& pFlow::coupling::particleGaussian::cellExtent()const
{
    if(cellExtent_.size() == mesh().nCells()) return cellExtent_;

    const auto& cells  = mesh().cells();
    const auto& faces  = mesh().faces();
    const auto& points = mesh().points();

    cellExtent_.setSize(mesh().nCells(), 0.0);

    forAll(cells, celli)
    {
        Foam::vector lo( Foam::GREAT,  Foam::GREAT,  Foam::GREAT);
        Foam::vector hi(-Foam::GREAT, -Foam::GREAT, -Foam::GREAT);

        for(const Foam::label facei : cells[celli])
        {
            for(const Foam::label pointi : faces[facei])
            {
                lo = Foam::min(lo, points[pointi]);
                hi = Foam::max(hi, points[pointi]);
            }
        }

        const Foam::vector d = hi - lo;

        cellExtent_[celli] =
            Foam::max(d.x(), Foam::max(d.y(), d.z()));
    }

    return cellExtent_;
}

Foam::scalar pFlow::coupling::particleGaussian::maxCellSize()const
{
    if(maxCellSize_ > 0) return maxCellSize_;

    const Foam::scalarField& ext = cellExtent();

    Foam::scalar eMax = 0;
    forAll(ext, celli) eMax = Foam::max(eMax, ext[celli]);

    maxCellSize_ = Foam::max(eMax, Foam::VSMALL);

    return maxCellSize_;
}

void pFlow::coupling::particleGaussian::markProcessorBoundaryCells()const
{
    const Foam::label nCells = mesh().nCells();

    if(static_cast<Foam::label>(procBoundaryCell_.size()) == nCells) return;

    procBoundaryCell_.assign(static_cast<size_t>(nCells), 0);

    const auto& bMesh = mesh().boundaryMesh();

    forAll(bMesh, patchi)
    {
        if(!Foam::isA<Foam::processorPolyPatch>(bMesh[patchi])) continue;

        for(const auto celli : bMesh[patchi].faceCells())
        {
            procBoundaryCell_[static_cast<size_t>(celli)] = 1;
        }
    }
}

Foam::scalar pFlow::coupling::particleGaussian::gaussianTailBeyond
(
    const Foam::scalar rOverSigma
)
{
    if(rOverSigma <= 0) return 1.0;

    const Foam::scalar r = rOverSigma;

    return Foam::erfc(r/Foam::sqrt(2.0))
         + Foam::sqrt(2.0/Foam::constant::mathematical::pi)*r*Foam::exp(-0.5*r*r);
}

void pFlow::coupling::particleGaussian::updateWeights
(
    const Plus::procCMField<real>& parDiameter
)
{
    if(!mesh().hasCellCells())
    {
        mesh().cellCells();
    }

    prepareScratch();
    markProcessorBoundaryCells();

    const auto& cellCells    = mesh().cellCells();
    const auto& parCellIndex = cMesh().parCellIndex();
    auto&       weights      = this->weights();
    const auto& centerMass   = weights.centerMass();
    const size_t numPar      = centerMass.size();

    const Foam::scalarField& cellV = mesh().cellVolumes();
    const Foam::scalarField& cellExt = cellExtent();
    const Foam::vectorField& cellC = mesh().cellCentres();

    Foam::label  layersUsed = 0;
    Foam::scalar procClearance = Foam::GREAT;
    Foam::scalar sigMax = 0;
    Foam::scalar dMax   = 0;

    #pragma omp parallel for schedule (dynamic) \
        reduction(max:layersUsed) reduction(min:procClearance) \
        reduction(max:sigMax) reduction(max:dMax)
    for(size_t i=0; i<numPar; i++)
    {
        auto& parWeights = weights[i];
        parWeights.clear();
        sigma_[i] = 0.0;

        Foam::label targetCellId = parCellIndex[i];

        if( targetCellId < 0 )
        {
            const realx3& cpSeed = centerMass[i];
            targetCellId = cMesh().findNearestCellTree
            (
                cpSeed,
                nSigma_*maxSigmaBound(parDiameter[i])
            );

            if( targetCellId < 0 ) continue;
        }

        #ifdef _OPENMP
            const size_t tid = static_cast<size_t>(omp_get_thread_num());
        #else
            const size_t tid = 0;
        #endif
        auto& stamp    = stamp_[tid];
        auto& frontier = frontier_[tid];
        auto& nextFr   = nextFrontier_[tid];

        const Foam::scalar dCell = cellExt[targetCellId];
        const Foam::scalar sig   = sigma(parDiameter[i], dCell);
        sigma_[i] = sig;

        const Foam::scalar twoS2   = 2.0*sig*sig;
        const Foam::scalar reach   = nSigma_*sig;
        const Foam::scalar reach2  = reach*reach;

        {
            const realx3& cpChk = centerMass[i];
            const Foam::vector dSeed
            (
                cpChk.x() - cellC[targetCellId].x(),
                cpChk.y() - cellC[targetCellId].y(),
                cpChk.z() - cellC[targetCellId].z()
            );

            if(parCellIndex[i] < 0 && Foam::magSqr(dSeed) > reach2)
            {
                sigma_[i] = 0.0;
                continue;
            }
        }
        const Foam::scalar expand  = reach + dCell;
        const Foam::scalar expand2 = expand*expand;

        const realx3& cp_i = centerMass[i];
        const Foam::vector cp{cp_i.x(), cp_i.y(), cp_i.z()};

        Foam::scalar pSubTotal   = 0.0;
        Foam::scalar minProcDist2 = Foam::GREAT;

        const Foam::label mark = ++stampCount_[tid];
        frontier.clear();
        frontier.push_back(targetCellId);
        stamp[targetCellId] = mark;

        Foam::label layer = 0;
        while(!frontier.empty() && layer < maxLayers_)
        {
            nextFr.clear();

            for(const auto celli : frontier)
            {
                const Foam::vector r = cp - cellC[celli];
                const Foam::scalar d2 = Foam::magSqr(r);

                if(d2 <= reach2)
                {
                    const Foam::scalar f = Foam::exp(-d2/twoS2)*cellV[celli];

                    parWeights.push_back({celli, f});
                    pSubTotal += f;

                    if(procBoundaryCell_[static_cast<size_t>(celli)])
                    {
                        minProcDist2 = Foam::min(minProcDist2, d2);
                    }
                }

                if(d2 <= expand2)
                {
                    for(const auto nbr : cellCells[celli])
                    {
                        if(stamp[nbr] != mark)
                        {
                            stamp[nbr] = mark;
                            nextFr.push_back(nbr);
                        }
                    }
                }
            }

            frontier.swap(nextFr);
            ++layer;
        }

        layersUsed = Foam::max(layersUsed, layer);
        sigMax     = Foam::max(sigMax, sig);
        dMax       = Foam::max(dMax, static_cast<Foam::scalar>(parDiameter[i]));

        if(minProcDist2 < Foam::GREAT && sig > 0)
        {
            procClearance =
                Foam::min(procClearance, Foam::sqrt(minProcDist2)/sig);
        }

        (void)pSubTotal;
    }

    if(layersUsed >= maxLayers_ && layersUsed > layersUsed_)
    {
        WARNING
            << "particleGaussian: the outward walk hit the maxLayers cap of "
            << maxLayers_ << ". The kernel is truncated and the deposited "
            << "source is renormalised over an incomplete support. Raise "
            << "particleGaussianInfo/maxLayers, or coarsen the coupling mesh."
            << END_WARNING;
    }
    layersUsed_ = Foam::max(layersUsed_, layersUsed);

    Foam::reduce(sigMax, Foam::maxOp<Foam::scalar>());
    Foam::reduce(dMax,   Foam::maxOp<Foam::scalar>());
    maxSigmaSeen_    = sigMax;
    maxDiameterSeen_ = dMax;

    checkProcessorTruncation(procClearance);

    if(macroLength_ > 0 && !warnedMacro_)
    {
        const Foam::scalar support = nSigma_*maxSigmaSeen_;

        if(support > 0.5*macroLength_ && Plus::processor::isMaster())
        {
            warnedMacro_ = true;

            WARNING
                << "Kernel support nSigma*sigma = " << support*1e3 << " mm is "
                << 100.0*support/macroLength_ << " % of macroscopicLength = "
                << macroLength_*1e3 << " mm.\n"
                << "      Anderson & Jackson (1967) require the averaging "
                   "radius to be small compared\n"
                << "      with the scale of macroscopic variation. The filter "
                   "is averaging over a\n"
                << "      sizeable fraction of the flow feature being resolved, "
                   "so the volume-averaged\n"
                << "      equations are being used at the edge of their "
                   "derivation."
                << END_WARNING;
        }
    }
}

void pFlow::coupling::particleGaussian::checkProcessorTruncation
(
    Foam::scalar procClearance
)const
{
    Foam::reduce(procClearance, Foam::minOp<Foam::scalar>());

    if(!Plus::processor::isParallel()) return;

    const Foam::scalar reach  = nSigma_*maxSigmaSeen_;
    const Foam::scalar shipped = shippingRadius();

    if(maxSigmaSeen_ <= 0 || shipped >= reach) return;

    if(reportedClearance_ < Foam::GREAT) return;
    reportedClearance_ = 0;

    std::ostringstream msg;
    msg << "particleGaussian: the kernel reaches farther than the particle is shipped.\n"
        << "      kernel support  nSigma*sigma      = " << reach*1e3 << " mm\n"
        << "      particle shipped to ranks within  = " << shipped*1e3 << " mm\n"
        << "        (particleMapping/domainExpansionRatio x d_max, with d_max = "
        << maxDiameterSeen_*1e3 << " mm)\n"
        << "      Ranks between those two radii are never sent the particle, so\n"
        << "      they contribute no weights and their share of the support is\n"
        << "      lost. The deposit and the sampled velocity then depend on the\n"
        << "      decomposition again.\n"
        << "      Raise particleMapping/domainExpansionRatio to at least "
        << std::ceil(reach/maxDiameterSeen_) << ".";

    if(truncationAbort_)
    {
        if(Plus::processor::isMaster())
        {
            fatalErrorInFunction << msg.str() << endl;
        }
        Plus::processor::abort(0);
    }
    else if(Plus::processor::isMaster())
    {
        WARNING << msg.str() << END_WARNING;
    }
}


void pFlow::coupling::particleGaussian::localWeightSums
(
    Plus::realProcCMField& sums
)const
{
    const auto& weights = this->weights();
    const size_t numPar = weights.size();

    for(size_t i=0; i<numPar; ++i)
    {
        Foam::scalar t = 0;
        for(const auto& [cellid, w] : weights[i]) t += w;
        sums[i] = static_cast<real>(t);
    }
}

void pFlow::coupling::particleGaussian::scaleWeights
(
    const Plus::realProcCMField& totals
)
{
    auto& weights = this->weights();
    const size_t numPar = weights.size();

    for(size_t i=0; i<numPar; ++i)
    {
        if(weights[i].empty()) continue;

        const Foam::scalar t = Foam::max
        (
            static_cast<Foam::scalar>(totals[i]),
            static_cast<Foam::scalar>(1.0e-300)
        );

        for(auto& [cellid, w] : weights[i]) w /= t;
    }
}
