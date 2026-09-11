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

// Maxey & Riley (1983) Phys. Fluids 26, 883; Lebedev (1976) USSR Comput. Math. Math. Phys. 16, 10

#include <cmath>
#include <cstddef>
#include <vector>

#include "sphereQuadrature.hpp"
#include "processorPolyPatch.H"
#include "emptyPolyPatch.H"

namespace pFlow::coupling
{

pFlow::coupling::sphereQuadrature::sphereQuadrature(const couplingMesh& cm)
:
    cMesh_(cm)
{
    const Foam::scalar c = 1.0/Foam::sqrt(3.0);
    int k = 0;

    for(int axis=0; axis<3; ++axis)
    {
        for(int sgn=0; sgn<2; ++sgn)
        {
            Foam::vector v(0,0,0);
            v[axis] = (sgn==0 ? 1.0 : -1.0);
            dir_[k] = v;
            wAng_[k] = 1.0/15.0;
            ++k;
        }
    }

    for(int sx=0; sx<2; ++sx)
    for(int sy=0; sy<2; ++sy)
    for(int sz=0; sz<2; ++sz)
    {
        dir_[k] = Foam::vector
        (
            (sx==0? c : -c), (sy==0? c : -c), (sz==0? c : -c)
        );
        wAng_[k] = 3.0/40.0;
        ++k;
    }

    rNode_[0] = 0.906179845938664;  wRad_[0] = 0.583666002653407;
    rNode_[1] = 0.538469310105684;  wRad_[1] = 0.416333997346593;
}

void pFlow::coupling::sphereQuadrature::resize(size_t nPar)
{
    pt_.assign(nPar*nPts, Foam::vector::zero);
    cell_.assign(nPar*nPts, -1);
    radius_.assign(nPar, 0.0);

    remoteFound_.assign(nPar*nPts, 0);
    remoteVal_.assign(nPar*nPts, Foam::vector::zero);
    located_ = false;
}

void pFlow::coupling::sphereQuadrature::prepare
(
    size_t              i,
    const Foam::vector& centre,
    Foam::scalar        a,
    Foam::label         hint
)
{
    radius_[i] = a;

    const size_t base = i*nPts;
    Foam::label  seed = hint;

    for(int q=0; q<nAng; ++q)
    {
        const Foam::vector x = centre + a*dir_[q];
        pt_[base+q]   = x;
        cell_[base+q] = cMesh_.findPointInCellTree(x, seed);
        if(cell_[base+q] >= 0) seed = cell_[base+q];
    }

    for(int r=0; r<nRad; ++r)
    {
        for(int q=0; q<nAng; ++q)
        {
            const size_t s = base + nAng*(1+r) + q;
            const Foam::vector x = centre + (a*rNode_[r])*dir_[q];
            pt_[s]   = x;
            cell_[s] = cMesh_.findPointInCellTree(x, seed);
            if(cell_[s] >= 0) seed = cell_[s];
        }
    }
}

void pFlow::coupling::sphereQuadrature::bind(const Foam::volVectorField& F)const
{
    const auto& mesh = cMesh_.mesh();

    mesh.cellCentres();
    mesh.faceCentres();
    mesh.tetBasePtIs();

    interp_.reset(new Foam::interpolationCellPoint<Foam::vector>(F));
    bound_ = &F;

    if(!located_) locateRemote();
    evaluateRemote();
}

// A node outside the owning rank's cells is looked for on every other rank;
// the lowest rank that holds it (a node on a shared face is found on both
// sides) interpolates the bound field there, and the value is summed back to
// the owner. A node that no rank holds lies outside the fluid, as for a
// particle touching a wall: it takes the boundary value of the bound field
// on the nearest boundary face (on a wall, its velocity), from the rank that
// owns that face. A node whose nearest face is on a cyclic patch stays out of
// the average.
void pFlow::coupling::sphereQuadrature::buildBoundaryTree()const
{
    const auto& mesh = cMesh_.mesh();
    const Foam::polyBoundaryMesh& patches = mesh.boundaryMesh();

    Foam::label nFaces = 0;

    forAll(patches, patchi)
    {
        const Foam::polyPatch& pp = patches[patchi];
        if(Foam::isA<Foam::processorPolyPatch>(pp) || Foam::isA<Foam::emptyPolyPatch>(pp)) continue;
        nFaces += pp.size();
    }

    if(nFaces == 0) return;

    Foam::labelList faces(nFaces);
    nFaces = 0;

    forAll(patches, patchi)
    {
        const Foam::polyPatch& pp = patches[patchi];
        if(Foam::isA<Foam::processorPolyPatch>(pp) || Foam::isA<Foam::emptyPolyPatch>(pp)) continue;

        forAll(pp, i)
        {
            faces[nFaces++] = pp.start() + i;
        }
    }

    boundaryTree_.reset
    (
        new Foam::indexedOctree<Foam::treeDataFace>
        (
            Foam::treeDataFace(true, mesh, faces),
            treeBoundBoxExtend(Foam::treeBoundBox(mesh.points()), 1.0e-3),
            8,
            10,
            3.0
        )
    );
}

void pFlow::coupling::sphereQuadrature::locateRemote()const
{
    located_ = true;

    missing_.clear();

    for(size_t i=0; i<radius_.size(); ++i)
    {
        if(radius_[i] <= 0) continue;           // not prepared here: owned elsewhere

        for(int q=0; q<nPts; ++q)
        {
            const size_t s = i*nPts + q;
            if(cell_[s] < 0) missing_.push_back(s);
        }
    }

    const Foam::label nProcs = Foam::Pstream::nProcs();
    const Foam::label me     = Foam::Pstream::myProcNo();

    Foam::List<Foam::pointField> all(nProcs);
    all[me].setSize(Foam::label(missing_.size()));

    forAll(all[me], k)
    {
        all[me][k] = pt_[missing_[k]];
    }

    Foam::Pstream::allGatherList(all);

    offsets_.setSize(nProcs + 1);
    offsets_[0] = 0;

    for(Foam::label p=0; p<nProcs; ++p)
    {
        offsets_[p+1] = offsets_[p] + all[p].size();
    }

    const Foam::label nTot = offsets_[nProcs];

    remotePts_.setSize(nTot);
    remoteCell_.setSize(nTot);
    remoteCell_ = -1;
    wallFace_.setSize(nTot);
    wallFace_ = -1;

    if(nTot == 0) return;

    Foam::labelList holder(nTot, nProcs);
    const Foam::boundBox& bb = cMesh_.mesh().bounds();

    for(Foam::label p=0; p<nProcs; ++p)
    {
        forAll(all[p], k)
        {
            const Foam::label g = offsets_[p] + k;
            remotePts_[g] = all[p][k];

            if(p == me || !bb.contains(all[p][k])) continue;

            const Foam::label c = cMesh_.findPointInCellTree(all[p][k], -1);

            if(c >= 0)
            {
                remoteCell_[g] = c;
                holder[g] = me;
            }
        }
    }

    Foam::Pstream::listCombineReduce(holder, Foam::minEqOp<Foam::label>());

    Foam::label nRemote = 0;
    Foam::label nOut    = 0;

    forAll(holder, g)
    {
        if(holder[g] != me) remoteCell_[g] = -1;
        if(holder[g] < nProcs) ++nRemote;
        else ++nOut;
    }

    // nodes outside the fluid: nearest boundary face over all ranks
    Foam::labelList onWall(nTot, 0);
    Foam::label nWall = 0;

    if(nOut > 0)
    {
        if(!boundaryTree_) buildBoundaryTree();

        const Foam::polyBoundaryMesh& patches = cMesh_.mesh().boundaryMesh();

        Foam::scalarList d2(nTot, Foam::GREAT);
        Foam::labelList  face(nTot, -1);
        Foam::pointField near(nTot, Foam::Zero);

        if(boundaryTree_)
        {
            forAll(holder, g)
            {
                if(holder[g] < nProcs) continue;

                const Foam::pointIndexHit hit =
                    boundaryTree_->findNearest(remotePts_[g], Foam::GREAT);

                if(hit.hit())
                {
                    d2[g]   = Foam::magSqr(hit.point() - remotePts_[g]);
                    face[g] = boundaryTree_->shapes().objectIndex(hit.index());
                    near[g] = hit.point();
                }
            }
        }

        Foam::scalarList d2Min(d2);
        Foam::Pstream::listCombineReduce(d2Min, Foam::minEqOp<Foam::scalar>());

        Foam::labelList faceHolder(nTot, nProcs);

        forAll(face, g)
        {
            if(face[g] >= 0 && d2[g] == d2Min[g]) faceHolder[g] = me;
        }

        Foam::Pstream::listCombineReduce(faceHolder, Foam::minEqOp<Foam::label>());

        const Foam::vectorField& Sf = cMesh_.mesh().faceAreas();

        forAll(faceHolder, g)
        {
            if(faceHolder[g] != me) continue;

            // outside the fluid: beyond the face along its outward normal
            const bool outside = (((remotePts_[g] - near[g]) & Sf[face[g]]) >= 0);

            if(outside && !patches[patches.whichPatch(face[g])].coupled())
            {
                wallFace_[g] = face[g];
                onWall[g] = 1;
            }
        }

        Foam::Pstream::listCombineReduce(onWall, Foam::maxEqOp<Foam::label>());

        forAll(onWall, g)
        {
            if(onWall[g]) ++nWall;
        }
    }

    for(size_t k=0; k<missing_.size(); ++k)
    {
        const Foam::label g = offsets_[me] + Foam::label(k);
        remoteFound_[missing_[k]] = (holder[g] < nProcs || onWall[g]);
    }

    if(nWall > 0 && !reportedWall_)
    {
        reportedWall_ = true;

        Foam::Info
            << "Faxen quadrature: " << nWall << " node(s) of particles touching "
            << "a wall lie outside the fluid and take the boundary value of the "
            << "field on the nearest face." << Foam::endl;
    }

    if(nRemote > 0 && !reportedRemote_)
    {
        reportedRemote_ = true;

        Foam::Info
            << "Faxen quadrature: " << nRemote << " node(s) of particles within "
            << "a radius of a processor boundary lie on the neighbouring rank "
            << "and are evaluated there." << Foam::endl;
    }
}

void pFlow::coupling::sphereQuadrature::evaluateRemote()const
{
    const Foam::label nTot = remotePts_.size();

    if(nTot == 0) return;

    Foam::List<Foam::vector> val(nTot, Foam::vector::zero);

    forAll(remoteCell_, g)
    {
        if(remoteCell_[g] >= 0)
        {
            val[g] = interp_->interpolate(remotePts_[g], remoteCell_[g]);
        }
    }

    const Foam::polyBoundaryMesh& patches = cMesh_.mesh().boundaryMesh();

    forAll(wallFace_, g)
    {
        if(wallFace_[g] >= 0)
        {
            const Foam::labelPair pf = patches.whichPatchFace(wallFace_[g]);
            val[g] = bound_->boundaryField()[pf.first()][pf.second()];
        }
    }

    Foam::Pstream::listCombineReduce(val, Foam::plusEqOp<Foam::vector>());

    const Foam::label me = Foam::Pstream::myProcNo();

    for(size_t k=0; k<missing_.size(); ++k)
    {
        remoteVal_[missing_[k]] = val[offsets_[me] + Foam::label(k)];
    }
}

void pFlow::coupling::sphereQuadrature::unbind()const
{
    interp_.reset();
    bound_ = nullptr;
}

Foam::vector pFlow::coupling::sphereQuadrature::interpolateAt
(
    const Foam::volVectorField& F,
    const Foam::vector&         x,
    Foam::label                 celli
)const
{
    if(bound_ == &F && interp_)
    {
        return interp_->interpolate(x, celli);
    }

    const auto& mesh = cMesh_.mesh();
    const auto& cc   = mesh.cellCentres();

    Foam::scalar wsum = 1.0/Foam::max(Foam::mag(x - cc[celli]), SMALL);
    Foam::vector vsum = F[celli]*wsum;

    const Foam::labelList& nbr = mesh.cellCells(celli);

    forAll(nbr, j)
    {
        const Foam::label cj = nbr[j];
        const Foam::scalar w = 1.0/Foam::max(Foam::mag(x - cc[cj]), SMALL);
        wsum += w;
        vsum += F[cj]*w;
    }

    return vsum/wsum;
}

Foam::vector pFlow::coupling::sphereQuadrature::surfaceAverage
(
    size_t                      i,
    const Foam::volVectorField& F
)const
{
    const size_t base = i*nPts;
    const bool   remote = (bound_ == &F);

    Foam::vector  acc(0,0,0);
    Foam::scalar  wsum = 0.0;

    for(int q=0; q<nAng; ++q)
    {
        const size_t s = base + q;
        const Foam::label c = cell_[s];

        if(c >= 0)
        {
            acc += wAng_[q]*interpolateAt(F, pt_[s], c);
        }
        else if(remote && remoteFound_[s])
        {
            acc += wAng_[q]*remoteVal_[s];
        }
        else
        {
            continue;
        }

        wsum += wAng_[q];
    }

    if(wsum <= 0) return Foam::vector::zero;

    if(wsum < 0.999) ++nIncomplete_;

    return acc/wsum;
}

Foam::vector pFlow::coupling::sphereQuadrature::volumeAverage
(
    size_t                      i,
    const Foam::volVectorField& F
)const
{
    const size_t base = i*nPts;
    const bool   remote = (bound_ == &F);

    Foam::vector  acc(0,0,0);
    Foam::scalar  wsum = 0.0;

    for(int r=0; r<nRad; ++r)
    {
        for(int q=0; q<nAng; ++q)
        {
            const size_t s = base + nAng*(1+r) + q;
            const Foam::label c = cell_[s];
            const Foam::scalar w = wRad_[r]*wAng_[q];

            if(c >= 0)
            {
                acc += w*interpolateAt(F, pt_[s], c);
            }
            else if(remote && remoteFound_[s])
            {
                acc += w*remoteVal_[s];
            }
            else
            {
                continue;
            }

            wsum += w;
        }
    }

    if(wsum <= 0) return Foam::vector::zero;

    if(wsum < 0.999) ++nIncomplete_;

    return acc/wsum;
}

}
