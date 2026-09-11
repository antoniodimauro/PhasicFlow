/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow       
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Modified in 2026 by Antonio Di Mauro (T-junction CFD-DEM study);
  changes: git diff cf6d3fb main, github.com/antoniodimauro/PhasicFlow
------------------------------------------------------------------------------  
Licence:
  This file is part of phasicFlow code. It is a free software for simulating 
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions. 
 
  phasicFlow is distributed to help others in their research in the field of 
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/

#include "porosity.hpp"
#include "PstreamReduceOps.H"
#include "procVectorPlus.hpp"
#include "procCommunicationPlus.hpp"
#include "streams.hpp"
#include "unresolvedCouplingSystem.hpp"


pFlow::coupling::porosity::porosity
(
	const unresolvedCouplingSystem& CS,
	const couplingMesh& 			cMesh,
	const Plus::realProcCMField& 	parDiam
)
:
	Foam::volScalarField
	(
		Foam::IOobject
	    (
	        "alpha",
	        Foam::timeName(cMesh.mesh().time()),
	        cMesh.mesh(),
	        Foam::IOobject::MUST_READ,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	cMesh.mesh()
	),
	alphaMin_(lookupDict<Foam::scalar>(CS.unresolvedDict().subDict("porosity"), "alphaMin")),
	uCS_(CS),
	particleDiameter_(parDiam),
	distribution_(CS.distribution())
{

}

const Foam::fvMesh &pFlow::coupling::porosity::mesh() const
{
    return uCS_.cMesh().mesh();
}


const pFlow::coupling::couplingMesh &pFlow::coupling::porosity::cMesh() const
{
    return uCS_.cMesh();
}

const pFlow::Plus::procCMField<Foam::label> &pFlow::coupling::porosity::parCellIndex() const
{
	return uCS_.parCellIndex();
}

pFlow::int32 pFlow::coupling::porosity::numInMesh() const
{
    return cMesh().numInMesh();
}

void pFlow::coupling::porosity::calculatePorosity()
{
	// Keep alpha^n as the old time. internalFieldUpdate writes through the plain field index, which stores no old
	// time, so without this the correctBoundaryConditions() below would store the new field and fvc::ddt(alpha) would
	// be zero (QUESTIONS.md G7). The time-index guard of storeOldTimes() makes the first call of a step store it;
	// later passes, the rollback and the write repeat keep it.
	this->storeOldTimes();

	this->internalFieldUpdate();
	this->correctBoundaryConditions();

	// A particle that is inserted or deleted in this step changes alpha with no displacement flow behind it: for that
	// step the old time is set to the new alpha, so it produces no dalpha/dt.
	// The count is of the particles whose centre is in a cell of this rank: each particle has one host cell, while
	// the rank arrays overlap (extended domains, domainDistribute.cpp), so their sizes change at a border crossing.
	const Foam::label tIdx = this->mesh().time().timeIndex();
	Foam::label nHosted = 0;
	{
		const auto& pci = parCellIndex();
		for(size_t i=0; i<pci.size(); ++i)
		{
			if(pci[i] >= 0) nHosted++;
		}
	}
	const Foam::label nNow = Foam::returnReduce(nHosted, Foam::sumOp<Foam::label>());

	if(tIdx != countTimeIndex_)
	{
		countChanged_   = (nPrevCount_ >= 0 && nNow != nPrevCount_);
		nPrevCount_     = nNow;
		countTimeIndex_ = tIdx;
	}

	if(countChanged_)
	{
		this->oldTime() == *this;
	}

	if(ddtReports_ < 3 && this->nOldTimes() > 0)
	{
		const Foam::scalar dMax = Foam::gMax
		(
			Foam::mag(this->primitiveField() - this->oldTime().primitiveField())()
		);
		Foam::Info<< "    porosity: max|alpha - alpha.oldTime()| = " << dMax
			<< (countChanged_ ? " (particle count changed: no dalpha/dt this step)" : "")
			<< Foam::endl;
		ddtReports_++;
	}
}


pFlow::uniquePtr<pFlow::coupling::porosity> 
pFlow::coupling::porosity::create
(
	const unresolvedCouplingSystem& CS,
	const couplingMesh& 			cMesh,
	const Plus::realProcCMField& 	parDiam
)
{

	Foam::word method = 
		CS.unresolvedDict().subDict("porosity").getOrDefault<Foam::word>("method", "distribution");

	if(method == "distribution")
	{
		method = "porosity<distribution>";
	}
		

	if( dictionaryvCtorSelector_.search(method))
	{
		Foam::Info<<"    Creating porosity method "<< Green_Text(method)<<" ...\n\n";
		return dictionaryvCtorSelector_[method] (CS, cMesh, parDiam);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< method << " dose not exist"
				" for porosity method in "<< CS.unresolvedDict().subDict("porosity").name()
				<<"\nAvaiable ones are: \n"
				,
				dictionaryvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;
}
