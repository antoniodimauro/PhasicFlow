/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
------------------------------------------------------------------------------
Licence:
  This file is part of phasicFlow code. It is a free software for simulating
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions.

  phasicFlow is distributed to help others in their research in the field of
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/

#include "insertionRegion.hpp"
#include "dictionary.hpp"
#include "insertion.hpp"
#include "particles.hpp"
#include "twoPartEntry.hpp"
#include "types.hpp"
#include "shape.hpp"

namespace pFlow
{
template<typename T>
bool
setOneEntry(const twoPartEntry& tpEntry, anyList& varList)
{
	if (getTypeName<T>() != tpEntry.firstPart())
		return false;

	T val = tpEntry.secondPartVal<T>();
	varList.emplaceBack(tpEntry.keyword(), val);

	return true;
}

bool
readOneEtrty(const dataEntry& entry, anyList& varList)
{
	twoPartEntry stField(entry);

	if (!(setOneEntry<real>(stField, varList) ||
	      setOneEntry<realx3>(stField, varList) ||
	      setOneEntry<realx4>(stField, varList) ||
	      setOneEntry<int8>(stField, varList) ||
	      setOneEntry<uint8>(stField, varList) ||
	      setOneEntry<uint32>(stField, varList) ||
	      setOneEntry<uint64>(stField, varList) ||
	      setOneEntry<int32>(stField, varList) ||
	      setOneEntry<int64>(stField, varList)))
	{
		fatalErrorInFunction << "un-supported data type " << stField.firstPart()
		                     << endl;
		return false;
	}
	return true;
}
}

bool
pFlow::insertionRegion::readInsertionRegion(const dictionary& dict)
{
	type_ = dict.getVal<word>("regionType");

	const bool hasCount = dict.containsDataEntry("numberPerInsertion");
	const bool hasRate  = dict.containsDataEntry("rate");

	if (hasCount && hasRate)
	{
		fatalErrorInFunction
		  << "both numberPerInsertion and rate are specified in dictionary "
		  << dict.globalName() << ". They are mutually exclusive: use "
		  << "numberPerInsertion for an exact count on every insertion event, "
		  << "or rate for a particle/s feed." << endl;
		return false;
	}

	if (hasCount)
	{
		countPerEvent_ = dict.getVal<uint32>("numberPerInsertion");
		if (countPerEvent_ == 0u)
		{
			fatalErrorInFunction
			  << "numberPerInsertion must be greater than zero in dictionary "
			  << dict.globalName() << ". Deactivate the region instead." << endl;
			return false;
		}
	}
	else if (hasRate)
	{
	rate_ = dict.getVal<real>("rate");
	}
	else
	{
		fatalErrorInFunction
		  << "neither numberPerInsertion nor rate is specified in dictionary "
		  << dict.globalName() << endl;
		return false;
	}

	pRegion_ = peakableRegion::create(type_, dict.subDict(type_ + "Info"));

	mixture_ = makeUnique<shapeMixture>(
	  dict.subDict("mixture"),
	  insertion_.Particles().getShapes().shapeNameList()
	);

	numInserted_ = mixture_().totalInserted();

	if (dict.containsDictionay("setFields"))
	{
		setFieldDict_ =
		  makeUnique<dictionary>("setFields", dict, dict.subDict("setFields"));
	}

	if (setFieldDict_)
	{
		if (!readSetFieldDict())
		{
			fatalErrorInFunction << "Error in reading dictionary "
			                     << setFieldDict_().globalName() << endl;
			return false;
		}
	}

	return true;
}

bool
pFlow::insertionRegion::writeInsertionRegion(dictionary& dict) const
{
	
	if(usesCountPerEvent())
	{
		if(!dict.add("numberPerInsertion", countPerEvent_))
			return false;
	}
	else if(!dict.add("rate", rate_))
		return false;

	if(!tControl_.write(dict))
		return false;

	if (!dict.add("regionType", type_))
		return false;
	if (pRegion_)
	{
		auto& prDict = dict.subDictOrCreate(type_ + "Info");
		if (!pRegion_().write(prDict))
			return false;
	}

	if (mixture_)
	{
		auto& mixDict = dict.subDictOrCreate("mixture");
		if (!mixture_().write(mixDict))
			return false;
	}

	if(setFieldDict_)
	{
	    if(!dict.addDict("setFields", setFieldDict_()))
			return false;	    
	}

	return true;
}

bool
pFlow::insertionRegion::readSetFieldDict()
{
	wordList Keys = setFieldDict_().dataEntryKeywords();

	for (const auto& key : Keys)
	{
		if (!readOneEtrty(setFieldDict_().dataEntryRef(key), setFieldList_))
		{
			return false;
		}
	}

	return true;
}

pFlow::insertionRegion::insertionRegion(
  const word&      name,
  const insertion& instn
)
  : name_(name),
    dict_(instn.subDict(name)),
    insertion_(instn),
    tControl_(dict_, "insertion")
{
	if (!readInsertionRegion(dict_))
	{
		fatalExit;
	}
}

const pFlow::pointStructure&
pFlow::insertionRegion::pStruct() const
{
	return Insertion().pStruct();
}

pFlow::uint32
pFlow::insertionRegion::numberToBeInserted(uint32 iter, real t, real dt)
{
	if (!tControl_.isInRange(iter, t, dt))
		return 0u;

	if (usesCountPerEvent())
	{
		uint64 k;
		if (tControl_.isTimeStep())
		{
			k = static_cast<uint64>(iter - tControl_.startIter()) /
			    static_cast<uint64>(tControl_.iInterval());
		}
		else
		{
			k = static_cast<uint64>(
			  (t - tControl_.startTime()) / tControl_.rInterval() + 0.5
			);
		}

		const uint64 expected = (k + 1u) * static_cast<uint64>(countPerEvent_);

		return expected > numInserted_
		         ? static_cast<uint32>(expected - numInserted_)
		         : 0u;
	}

	if (tControl_.isTimeStep())
	{
		return static_cast<uint32>(
		  (iter - tControl_.startIter() + tControl_.iInterval()) * dt * rate_ -
		  numInserted_
		);
	}
	else
	{
		return static_cast<uint32>(
		  (t - tControl_.startTime() + tControl_.rInterval()) * rate_ -
		  numInserted_
		);
	}
}
