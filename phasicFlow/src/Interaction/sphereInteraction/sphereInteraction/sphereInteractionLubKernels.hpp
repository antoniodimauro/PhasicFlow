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

// Costa et al. (2015) Phys. Rev. E 92, 053012

#ifndef __sphereInteractionLubKernels_hpp__
#define __sphereInteractionLubKernels_hpp__

#include "sphereInteractionKernels.hpp"
#include "lubrication.hpp"

namespace pFlow::sphereInteractionKernels
{

template<
	typename ContactForceModel, 
	typename ContactListType>
struct ppInteractionLubFunctor
{

	using PairType = typename ContactListType::PairType;
	using ValueType = typename ContactListType::ValueType;

	real dt_;

	ContactForceModel forceModel_;
	ContactListType   tobeFilled_;

	deviceViewType1D<real>  		diam_;
	deviceViewType1D<uint32> 		propId_;
	deviceViewType1D<realx3>  		pos_;
	deviceViewType1D<realx3>  		lVel_;
	deviceViewType1D<realx3>  		rVel_;
	deviceViewType1D<realx3>	cForce_;
	deviceViewType1D<realx3>	cTorque_;

	deviceViewType1D<real>		lubScale_;

	lubricationModels::lubrication	lub_;

	ppInteractionLubFunctor(
		real dt,
		ContactForceModel 		forceModel,
		ContactListType 		tobeFilled,
		deviceViewType1D<real>  		diam,
		deviceViewType1D<uint32> 		propId,
		deviceViewType1D<realx3>  		pos,
		deviceViewType1D<realx3>  		lVel,
		deviceViewType1D<realx3>  		rVel,
		deviceViewType1D<realx3>		cForce,
		deviceViewType1D<realx3> 		cTorque,
		deviceViewType1D<real>  		lubScale,
		const lubricationModels::lubrication& lub )
	:
		dt_(dt),
		forceModel_(forceModel),
		tobeFilled_(tobeFilled),
		diam_(diam),
		propId_(propId),
		pos_(pos),
		lVel_(lVel),
		rVel_(rVel),
		cForce_(cForce),
		cTorque_(cTorque),
		lubScale_(lubScale),
		lub_(lub)
	{}

	INLINE_FUNCTION_HD
	void operator()(const uint32 n)const
	{
		
		if(!tobeFilled_.isValid(n))return;

		auto [i,j] = tobeFilled_.getPair(n);
		
		real Ri = 0.5*diam_[i];
		real Rj = 0.5*diam_[j];
		realx3 xi = pos_[i];
		realx3 xj = pos_[j];
		real dist = length(xj-xi);
		real ovrlp = (Ri+Rj) - dist;
		
		if( ovrlp >0.0 )
		{
			
			auto Vi = lVel_[i];
			auto Vj = lVel_[j];
			auto wi = rVel_[i];
			auto wj = rVel_[j];
			auto Nij = (xj-xi)/dist;
			auto Vr = Vi - Vj + cross((Ri*wi+Rj*wj), Nij);
			
			auto history = tobeFilled_.getValue(n);

			int32 propId_i = propId_[i];
			int32 propId_j = propId_[j];

			realx3 FCn, FCt, Mri, Mrj, Mij, Mji;

			// calculates contact force 
			forceModel_.contactForce(
				dt_, i, j,
				propId_i, propId_j,
				Ri, Rj,
				ovrlp,
				Vr, Nij,
				history,
				FCn, FCt
				);

			forceModel_.rollingFriction(
				dt_, i, j,
				propId_i, propId_j,
				Ri, Rj,
				wi, wj,
				Nij,
				FCn,
				Mri, Mrj
				);

			auto M = cross(Nij,FCt);
			Mij = Ri*M+Mri;
			Mji = Rj*M+Mrj;
			
			auto FC = FCn + FCt;
			
			Kokkos::atomic_add(&cForce_[i].x_,FC.x_);
			Kokkos::atomic_add(&cForce_[i].y_,FC.y_);
			Kokkos::atomic_add(&cForce_[i].z_,FC.z_);

			Kokkos::atomic_add(&cForce_[j].x_,-FC.x_);
			Kokkos::atomic_add(&cForce_[j].y_,-FC.y_);
			Kokkos::atomic_add(&cForce_[j].z_,-FC.z_);

			Kokkos::atomic_add(&cTorque_[i].x_, Mij.x_);
			Kokkos::atomic_add(&cTorque_[i].y_, Mij.y_);
			Kokkos::atomic_add(&cTorque_[i].z_, Mij.z_);

			Kokkos::atomic_add(&cTorque_[j].x_, Mji.x_);
			Kokkos::atomic_add(&cTorque_[j].y_, Mji.y_);
			Kokkos::atomic_add(&cTorque_[j].z_, Mji.z_);
			

			tobeFilled_.setValue(n,history);

		}
		else
		{
			tobeFilled_.setValue(n, ValueType());

			if( lub_.active() )
			{
				const auto Nij = (xj-xi)/dist;
				const auto Vr  = lVel_[i] - lVel_[j];

				const auto FL = lub_.force(
					Ri, Rj,
					-ovrlp,
					Vr, Nij,
					min(lubScale_[i], lubScale_[j]));

				Kokkos::atomic_add(&cForce_[i].x_, FL.x_);
				Kokkos::atomic_add(&cForce_[i].y_, FL.y_);
				Kokkos::atomic_add(&cForce_[i].z_, FL.z_);

				Kokkos::atomic_add(&cForce_[j].x_,-FL.x_);
				Kokkos::atomic_add(&cForce_[j].y_,-FL.y_);
				Kokkos::atomic_add(&cForce_[j].z_,-FL.z_);
			}
		}

	}
};


template<typename ContactListType>
struct ppLubricationSumFunctor
{
	ContactListType   		tobeFilled_;
	deviceViewType1D<real>  	diam_;
	deviceViewType1D<realx3>  	pos_;
	deviceViewType1D<real>  	cSum_;
	lubricationModels::lubrication	lub_;

	ppLubricationSumFunctor(
		ContactListType 		tobeFilled,
		deviceViewType1D<real>  	diam,
		deviceViewType1D<realx3>  	pos,
		deviceViewType1D<real>  	cSum,
		const lubricationModels::lubrication& lub )
	:
		tobeFilled_(tobeFilled),
		diam_(diam),
		pos_(pos),
		cSum_(cSum),
		lub_(lub)
	{}

	INLINE_FUNCTION_HD
	void operator()(const uint32 n)const
	{
		if(!tobeFilled_.isValid(n))return;

		auto [i,j] = tobeFilled_.getPair(n);

		const real Ri = 0.5*diam_[i];
		const real Rj = 0.5*diam_[j];

		const real dist = length(pos_[j]-pos_[i]);
		const real gap  = dist - Ri - Rj;

		const real c = lub_.coefficient(Ri, Rj, gap);

		if(c <= 0) return;

		Kokkos::atomic_add(&cSum_[i], c);
		Kokkos::atomic_add(&cSum_[j], c);
	}
};

template<typename MassView>
struct ppLubricationScaleFunctor
{
	real					dt_;
	MassView				mass_;

	deviceViewType1D<real>			added_;
	bool					haveAdded_ = false;

	deviceViewType1D<real>  		cSum_;
	lubricationModels::lubrication		lub_;

	ppLubricationScaleFunctor(
		real dt,
		MassView mass,
		deviceViewType1D<real> added,
		deviceViewType1D<real> cSum,
		const lubricationModels::lubrication& lub )
	:
		dt_(dt), mass_(mass), added_(added),
		haveAdded_(added.size() >= mass.size() && mass.size() > 0),
		cSum_(cSum), lub_(lub)
	{}

	INLINE_FUNCTION_HD
	void operator()(const uint32 i)const
	{
		const real sum = cSum_[i];

		if(sum <= 0)
		{
			cSum_[i] = 1.0;
			return;
		}

		const real m = haveAdded_ ? (mass_[i] + added_[i]) : mass_[i];

		const real cMax = lub_.maxCoeffSum(m, dt_);

		cSum_[i] = min(static_cast<real>(1.0), cMax/sum);
	}
};


}

#endif
