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

#ifndef __sphereQuadrature_hpp__
#define __sphereQuadrature_hpp__

#include <memory>
#include <vector>

#include "OFCompatibleHeader.hpp"
#include "couplingMesh.hpp"
#include "interpolationCellPoint.H"
#include "indexedOctree.H"
#include "treeDataFace.H"

namespace pFlow::coupling
{

class sphereQuadrature
{
public:

    static constexpr int nAng = 14;
    static constexpr int nRad = 2;

    static constexpr int nPts = nAng*(1 + nRad);

private:

    const couplingMesh&         cMesh_;

    Foam::vector                dir_[nAng];
    Foam::scalar                wAng_[nAng];
    Foam::scalar                rNode_[nRad];
    Foam::scalar                wRad_[nRad];

    std::vector<Foam::vector>   pt_;
    std::vector<Foam::label>    cell_;
    std::vector<Foam::scalar>   radius_;

    mutable Foam::label         nIncomplete_ = 0;

    mutable const Foam::volVectorField* bound_ = nullptr;
    mutable std::unique_ptr<Foam::interpolationCellPoint<Foam::vector>> interp_;

    // Nodes outside this rank's cells, evaluated on the rank that holds them;
    // nodes outside the fluid take the boundary value on the nearest face.
    // Global list: the missing nodes of all ranks, in rank order.
    mutable bool                located_ = false;
    mutable bool                reportedRemote_ = false;
    mutable bool                reportedWall_ = false;
    mutable std::vector<size_t> missing_;       // this rank's node slots not found here
    mutable Foam::labelList     offsets_;       // start of each rank's nodes in the global list
    mutable Foam::pointField    remotePts_;     // global list of nodes
    mutable Foam::labelList     remoteCell_;    // cell here of each global node this rank evaluates, -1 otherwise
    mutable Foam::labelList     wallFace_;      // boundary face here whose value a global node takes, -1 otherwise
    mutable std::vector<char>   remoteFound_;   // per node slot: held by another rank, or on a wall
    mutable std::vector<Foam::vector> remoteVal_; // per node slot: value there of the bound field

    // boundary faces of this rank, processor and empty patches excluded
    mutable std::unique_ptr<Foam::indexedOctree<Foam::treeDataFace>> boundaryTree_;

    void buildBoundaryTree()const;

    void locateRemote()const;

    void evaluateRemote()const;

    Foam::vector interpolateAt(
        const Foam::volVectorField& F,
        const Foam::vector&         x,
        Foam::label                 celli)const;

public:

    explicit sphereQuadrature(const couplingMesh& cm);

    void resize(size_t nPar);

    void prepare(size_t i, const Foam::vector& centre, Foam::scalar a, Foam::label hint);

    Foam::vector surfaceAverage(size_t i, const Foam::volVectorField& F)const;

    Foam::vector volumeAverage(size_t i, const Foam::volVectorField& F)const;

    Foam::label countIncomplete()const { return nIncomplete_; }

    // collective: every rank calls it, in the same order
    void bind(const Foam::volVectorField& F)const;
    void unbind()const;

    void resetCount()const { nIncomplete_ = 0; }
};

}

#endif
