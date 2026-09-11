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

#include <algorithm>

#include "AdamsMoulton4PEC.hpp"
#include "pointStructure.hpp"
#include "Time.hpp"
#include "vocabs.hpp"

namespace pFlow
{

using rpIntegration = Kokkos::RangePolicy<
		DefaultExecutionSpace,
		Kokkos::Schedule<Kokkos::Static>,
		Kokkos::IndexType<uint32>
		>;

bool AM4_predictAllActive(
	real dt,
	realx3Field_D& y,
	realx3PointField_D& dy,
	realx3PointField_D& dy0,
	realx3PointField_D& dy1,
	realx3PointField_D& dy2,
	realx3PointField_D& dy3)
{
	auto d_y   = y.deviceView();
	auto d_dy  = dy.deviceView();
	auto d_dy0 = dy0.deviceView();
	auto d_dy1 = dy1.deviceView();
	auto d_dy2 = dy2.deviceView();
	auto d_dy3 = dy3.deviceView();
	auto activeRng = dy0.activeRange();

	Kokkos::parallel_for(
		"AdamsMoulton4PEC::predictAllActive",
		rpIntegration(activeRng.start(), activeRng.end()),
		LAMBDA_HD(uint32 i){
			d_y[i] += dt * (
				  static_cast<real>(55.0/24.0) * d_dy[i]
				- static_cast<real>(59.0/24.0) * d_dy1[i]
				+ static_cast<real>(37.0/24.0) * d_dy2[i]
				- static_cast<real>( 9.0/24.0) * d_dy3[i]);
			d_dy0[i] = d_dy[i];
		});
	Kokkos::fence();

	return true;	
}

bool AM4_predictScattered(
	real dt, 
	realx3Field_D& y,
	realx3PointField_D& dy,
	realx3PointField_D& dy0,
	realx3PointField_D& dy1,
	realx3PointField_D& dy2,
	realx3PointField_D& dy3)
{
	auto d_y   = y.deviceView();
	auto d_dy  = dy.deviceView();
	auto d_dy0 = dy0.deviceView();
	auto d_dy1 = dy1.deviceView();
	auto d_dy2 = dy2.deviceView();
	auto d_dy3 = dy3.deviceView();
	auto activeRng = dy0.activeRange();
	const auto& activeP = dy0.activePointsMaskDevice();

	Kokkos::parallel_for(
		"AdamsMoulton4PEC::predictScattered",
		rpIntegration(activeRng.start(), activeRng.end()),
		LAMBDA_HD(uint32 i){
			if(activeP(i))
			{
				d_y[i] += dt * (
					  static_cast<real>(55.0/24.0) * d_dy[i]
					- static_cast<real>(59.0/24.0) * d_dy1[i]
					+ static_cast<real>(37.0/24.0) * d_dy2[i]
					- static_cast<real>( 9.0/24.0) * d_dy3[i]);
				d_dy0[i] = d_dy[i];
			}
		});
	Kokkos::fence();
	
	return true;	
}
bool AM4_correctAllActive(
	real dt,
	realx3Field_D& y,
	realx3PointField_D& dy,
	realx3PointField_D& dy0,
	realx3PointField_D& dy1,
	realx3PointField_D& dy2,
	realx3PointField_D& dy3,
	real damping = 1.0)
{
	auto d_y   = y.deviceView();
	auto d_dy  = dy.deviceView();
	auto d_dy0 = dy0.deviceView();
	auto d_dy1 = dy1.deviceView();
	auto d_dy2 = dy2.deviceView();
	auto d_dy3 = dy3.deviceView();
	auto activeRng = dy0.activeRange();

	Kokkos::parallel_for(
		"AdamsMoulton4PEC::correctAllActive",
		rpIntegration(activeRng.start(), activeRng.end()),
		LAMBDA_HD(uint32 i){
			d_y[i] += damping * dt * (
				  static_cast<real>( 9.0/24.0) * d_dy[i]
				- static_cast<real>(36.0/24.0) * d_dy0[i]
				+ static_cast<real>(54.0/24.0) * d_dy1[i]
				- static_cast<real>(36.0/24.0) * d_dy2[i]
				+ static_cast<real>( 9.0/24.0) * d_dy3[i]);

			d_dy3[i] = d_dy2[i];
			d_dy2[i] = d_dy1[i];
			d_dy1[i] = d_dy0[i];
			d_dy0[i] = d_dy[i];
		});
	Kokkos::fence();

	return true;
}

bool AM4_correctScattered(
	real dt,
	realx3Field_D& y,
	realx3PointField_D& dy,
	realx3PointField_D& dy0,
	realx3PointField_D& dy1,
	realx3PointField_D& dy2,
	realx3PointField_D& dy3,
	real damping = 1.0)
{
	auto d_y   = y.deviceView();
	auto d_dy  = dy.deviceView();
	auto d_dy0 = dy0.deviceView();
	auto d_dy1 = dy1.deviceView();
	auto d_dy2 = dy2.deviceView();
	auto d_dy3 = dy3.deviceView();
	auto activeRng = dy0.activeRange();
	const auto& activeP = dy0.activePointsMaskDevice();

	Kokkos::parallel_for(
		"AdamsMoulton4PEC::correctScattered",
		rpIntegration(activeRng.start(), activeRng.end()),
		LAMBDA_HD(uint32 i){
			if(activeP(i))
			{
				d_y[i] += damping * dt * (
					  static_cast<real>( 9.0/24.0) * d_dy[i]
					- static_cast<real>(36.0/24.0) * d_dy0[i]
					+ static_cast<real>(54.0/24.0) * d_dy1[i]
					- static_cast<real>(36.0/24.0) * d_dy2[i]
					+ static_cast<real>( 9.0/24.0) * d_dy3[i]);

				d_dy3[i] = d_dy2[i];
				d_dy2[i] = d_dy1[i];
				d_dy1[i] = d_dy0[i];
				d_dy0[i] = d_dy[i];
			}
		});
	Kokkos::fence();

	return true;
}

}

pFlow::AdamsMoulton4PEC::AdamsMoulton4PEC
(
	const word& baseName,
	pointStructure& pStruct,
	const word& method,
	const realx3Field_D& initialValField,
	bool  keepHistory
)
:
	integration(baseName, pStruct, method, initialValField, keepHistory),
	dy0_
	(
		objectFile
		(
			groupNames(baseName, "dy0"),
			pStruct.time().integrationFolder(),
			objectFile::READ_IF_PRESENT,
			keepHistory ? objectFile::WRITE_ALWAYS : objectFile::WRITE_NEVER
		),
		pStruct,
		zero3,
		zero3
	),
	dy1_
	(
		objectFile
		(
			groupNames(baseName, "dy1"),
			pStruct.time().integrationFolder(),
			objectFile::READ_IF_PRESENT,
			keepHistory ? objectFile::WRITE_ALWAYS : objectFile::WRITE_NEVER
		),
		pStruct,
		zero3,
		zero3
	),
	dy2_
	(
		objectFile
		(
			groupNames(baseName, "dy2"),
			pStruct.time().integrationFolder(),
			objectFile::READ_IF_PRESENT,
			keepHistory ? objectFile::WRITE_ALWAYS : objectFile::WRITE_NEVER
		),
		pStruct,
		zero3,
		zero3
	),
	dy3_
	(
		objectFile
		(
			groupNames(baseName, "dy3"),
			pStruct.time().integrationFolder(),
			objectFile::READ_IF_PRESENT,
			keepHistory ? objectFile::WRITE_ALWAYS : objectFile::WRITE_NEVER
		),
		pStruct,
		zero3,
		zero3
	),
	initialValField_(initialValField),
	boundaryList_(pStruct, method, *this)
{
}

void pFlow::AdamsMoulton4PEC::updateBoundariesSlaveToMasterIfRequested()
{
	dy0_.updateBoundariesSlaveToMasterIfRequested();
	dy1_.updateBoundariesSlaveToMasterIfRequested();
	dy2_.updateBoundariesSlaveToMasterIfRequested();
	dy3_.updateBoundariesSlaveToMasterIfRequested();
}

bool pFlow::AdamsMoulton4PEC::predict
(
	real dt,
	realx3PointField_D& y,
	realx3PointField_D& dy
)
{
	bool success = false;
	if(dy0_.isAllActive())
	{
		success = AM4_predictAllActive(dt, y.field(), dy, dy0_, dy1_, dy2_, dy3_);
	}
	else
	{
		success = AM4_predictScattered(dt, y.field(), dy, dy0_, dy1_, dy2_, dy3_);
	}
	return success;
}

bool pFlow::AdamsMoulton4PEC::predict
(
	real dt,
	realx3Field_D& y,
	realx3PointField_D& dy
)
{
	bool success = false;
	if(dy0_.isAllActive())
	{
		success = AM4_predictAllActive(dt, y, dy, dy0_, dy1_, dy2_, dy3_);
	}
	else
	{
		success = AM4_predictScattered(dt, y, dy, dy0_, dy1_, dy2_, dy3_);
	}
	return success;
}

bool pFlow::AdamsMoulton4PEC::correct
(
	real dt,
	realx3PointField_D& y,
	realx3PointField_D& dy,
	real damping
)
{
	bool success = false;
	if(dy0_.isAllActive())
	{
		success = AM4_correctAllActive(dt, y.field(), dy, dy0_, dy1_, dy2_, dy3_, damping);
	}
	else
	{
		success = AM4_correctScattered(dt, y.field(), dy, dy0_, dy1_, dy2_, dy3_, damping);
	}

	success = success && boundaryList_.correct(dt, y, dy);

	return success;
}

bool pFlow::AdamsMoulton4PEC::correctPStruct
(
	real dt,
	pointStructure& pStruct,
	realx3PointField_D& vel
)
{
	bool success = false;
	if(dy0_.isAllActive())
	{
		success = AM4_correctAllActive(dt, pStruct.pointPosition(), vel, dy0_, dy1_, dy2_, dy3_);
	}
	else
	{
		success = AM4_correctScattered(dt, pStruct.pointPosition(), vel, dy0_, dy1_, dy2_, dy3_);
	}

	success = success && boundaryList_.correctPStruct(dt, pStruct, vel);

	return success;
}

void pFlow::AdamsMoulton4PEC::resetHistory()
{
	dy0_.fill(zero3);
	dy1_.fill(zero3);
	dy2_.fill(zero3);
	dy3_.fill(zero3);
}

void pFlow::AdamsMoulton4PEC::saveHistory()
{
	auto v0 = dy0_.deviceViewAll();
	auto v1 = dy1_.deviceViewAll();
	auto v2 = dy2_.deviceViewAll();
	auto v3 = dy3_.deviceViewAll();

	svSize_ = static_cast<uint32>(v0.size());

	if( sv0_.size() < svSize_ )
	{
		Kokkos::realloc(sv0_, svSize_);
		Kokkos::realloc(sv1_, svSize_);
		Kokkos::realloc(sv2_, svSize_);
		Kokkos::realloc(sv3_, svSize_);
	}

	using pair_t = Kokkos::pair<size_t,size_t>;
	const pair_t rng(0, static_cast<size_t>(svSize_));

	Kokkos::deep_copy(Kokkos::subview(sv0_, rng), v0);
	Kokkos::deep_copy(Kokkos::subview(sv1_, rng), v1);
	Kokkos::deep_copy(Kokkos::subview(sv2_, rng), v2);
	Kokkos::deep_copy(Kokkos::subview(sv3_, rng), v3);
}

void pFlow::AdamsMoulton4PEC::restoreHistory()
{
	if( svSize_ == 0 )
	{
		resetHistory();
		return;
	}

	auto v0 = dy0_.deviceViewAll();
	auto v1 = dy1_.deviceViewAll();
	auto v2 = dy2_.deviceViewAll();
	auto v3 = dy3_.deviceViewAll();

	const size_t n = std::min<size_t>(svSize_, v0.size());

	using pair_t = Kokkos::pair<size_t,size_t>;
	const pair_t rng(0, n);

	Kokkos::deep_copy(Kokkos::subview(v0, rng), Kokkos::subview(sv0_, rng));
	Kokkos::deep_copy(Kokkos::subview(v1, rng), Kokkos::subview(sv1_, rng));
	Kokkos::deep_copy(Kokkos::subview(v2, rng), Kokkos::subview(sv2_, rng));
	Kokkos::deep_copy(Kokkos::subview(v3, rng), Kokkos::subview(sv3_, rng));
}