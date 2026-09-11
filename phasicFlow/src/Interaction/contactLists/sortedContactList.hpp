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

#ifndef __sortedContactList_hpp__ 
#define __sortedContactList_hpp__

#include "sortedPairs.hpp"

namespace pFlow
{

template<typename valueType, typename executionSpace, typename idType>
class sortedContactList
:
	public sortedPairs<executionSpace, idType>
{
public:

	using SortedPairs = sortedPairs<executionSpace, idType>;

	using ValueType 		= valueType;

	using IdType 			= typename SortedPairs::IdType;

	using ExecutionSpace 	= typename SortedPairs::ExecutionSpace;

	using memory_space 		= typename SortedPairs::memory_space;

	using PairType 			= typename SortedPairs::PairType;

	using ContainerType 	= typename SortedPairs::ContainerType;

	class TagReFillPairs{};

protected:

	ViewType1D<ValueType,ExecutionSpace> 	values_;

	uint32 									size0_ = 0;

	ViewType1D<PairType,ExecutionSpace> 	sortedPairs0_;

	ViewType1D<ValueType,ExecutionSpace> 	values0_;

	// C8: the pairs and contact histories at the start of a fluid step, restored before every repeated coupling pass
	// so that the tangential springs advance once per step (QUESTIONS.md G9)
	ViewType1D<PairType,ExecutionSpace> 	savedPairs_;

	ViewType1D<ValueType,ExecutionSpace> 	savedValues_;

	uint32 									savedSize_ = 0;

	bool 									hasSaved_ = false;


	void adjustCapacity()
	{
		if(auto s = this->size(); s > values_.size())
		{
			reallocNoInit(values_, s);
		}
	}

	

	using rpReFillPairs = Kokkos::RangePolicy<
		ExecutionSpace,
		Kokkos::Schedule<Kokkos::Static>,
		Kokkos::IndexType<uint32>,
		TagReFillPairs>;

public:

	TypeInfoNV("sortedContactList");

	sortedContactList(uint32 initialSize =1)
	:
		sortedContactList("sortedContactList", initialSize)
	{}
	
	sortedContactList(const word& name, uint32 initialSize =1)
	:
		SortedPairs(name, initialSize),
		values_(groupNames(name, "values"), SortedPairs::capacity()),
		sortedPairs0_(groupNames(name, "sortedPairs0"), SortedPairs::capacity()),
		values0_(groupNames(name, "values0"), SortedPairs::capacity())
	{}

	bool beforeBroadSearch()
	{
		swapViews(values0_, values_);
		swapViews(sortedPairs0_, this->sortedPairs_);
		size0_ = this->size();
		return SortedPairs::beforeBroadSearch();
	}

	bool afterBroadSearch()
	{
		SortedPairs::afterBroadSearch();
		
		adjustCapacity();

		Kokkos::parallel_for(
			"sortedContactList::reFillPairs",
			rpReFillPairs(0, this->size() ),
			*this
			);
		Kokkos::fence();

		return true;
	}

	void saveState()
	{
		savedSize_ = this->size();
		const uint32 n = max(savedSize_, 1u);
		if(savedPairs_.size() < n)  reallocNoInit(savedPairs_, n);
		if(savedValues_.size() < n) reallocNoInit(savedValues_, n);
		if(savedSize_ > 0)
		{
			const auto rng = Kokkos::make_pair(0u, savedSize_);
			Kokkos::deep_copy(Kokkos::subview(savedPairs_, rng),  Kokkos::subview(this->sortedPairs_, rng));
			Kokkos::deep_copy(Kokkos::subview(savedValues_, rng), Kokkos::subview(values_, rng));
		}
		hasSaved_ = true;
	}

	void restoreState()
	{
		if(!hasSaved_) return;
		if(this->sortedPairs_.size() < savedSize_) reallocNoInit(this->sortedPairs_, savedSize_);
		if(values_.size() < savedSize_)            reallocNoInit(values_, savedSize_);
		if(savedSize_ > 0)
		{
			const auto rng = Kokkos::make_pair(0u, savedSize_);
			Kokkos::deep_copy(Kokkos::subview(this->sortedPairs_, rng), Kokkos::subview(savedPairs_, rng));
			Kokkos::deep_copy(Kokkos::subview(values_, rng),            Kokkos::subview(savedValues_, rng));
		}
		this->size_ = savedSize_;
	}

	INLINE_FUNCTION_HD
	ValueType getValue(uint32 idx)const
	{
		return values_[idx];
	}

	INLINE_FUNCTION_HD
	void setValue(uint32 idx, const ValueType& val)const
	{
		values_[idx] = val;
	}

	INLINE_FUNCTION_HD
	void operator()(TagReFillPairs, uint32 idx)const
	{

		uint32 searchLen = max(size0_/1000u,10u);
		uint32 start 	= idx-min(searchLen,idx);
		uint32 end 	= min(size0_,idx+searchLen);
		auto newPair = this->sortedPairs_[idx];
		if(	auto idx0 = binarySearch(
				sortedPairs0_,
				start,
				end,
				newPair);
				idx0!=static_cast<uint32>(-1))
		{
			values_[idx] = values0_[idx0];
		}
		else if(auto idx0 = binarySearch(
			sortedPairs0_,
			start, 
			end,
			newPair);
			idx0!= static_cast<uint32>(-1) )
		{
			values_[idx] = values0_[idx0];
			
		}
		else
		{
			values_[idx] = ValueType();
		}
	}



}; // sortedContactList



} // pFlow



#endif //__sortedContactList_hpp__
