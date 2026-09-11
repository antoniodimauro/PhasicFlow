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

#ifndef __AdamsMoulton4PEC_hpp__
#define __AdamsMoulton4PEC_hpp__


#include "integration.hpp"
#include "pointFields.hpp"
#include "boundaryIntegrationList.hpp"

namespace pFlow
{

class AdamsMoulton4PEC
:
	public integration
{
private:

	realx3PointField_D dy0_;

	realx3PointField_D dy1_;

	realx3PointField_D dy2_;

	realx3PointField_D dy3_;

	const realx3Field_D& initialValField_;

	deviceViewType1D<realx3> sv0_, sv1_, sv2_, sv3_;

	uint32                   svSize_ = 0;

	boundaryIntegrationList boundaryList_;

public:

	/// Type info
	TypeInfo("AdamsMoulton4PEC");

	// - Constructors
		
		/// Construct from components
		AdamsMoulton4PEC(
			const word& baseName,
			pointStructure& pStruct,
			const word& method,
			const realx3Field_D& initialValField,
			bool  keepHistory);

		/// Destructor 
		~AdamsMoulton4PEC() override = default;

		add_vCtor(
			integration,
			AdamsMoulton4PEC,
			word);

	// - Methods

		void updateBoundariesSlaveToMasterIfRequested() override;

		word method() const override
		{
			return "AdamsMoulton4PEC";
		}

		bool predict(
			real dt, 
			realx3PointField_D& y,
			realx3PointField_D& dy) override;

		bool predict(
			real dt,
			realx3Field_D& y,
			realx3PointField_D& dy) override;

		bool correct(
			real dt, 
			realx3PointField_D& y,
			realx3PointField_D& dy,
			real damping = 1.0) override;

		bool correctPStruct(
			real dt,
			pointStructure& pStruct,
			realx3PointField_D& vel) override;

		struct deviceHistory
		{
			deviceViewType1D<realx3> dy0, dy1, dy2, dy3;
		};

		deviceHistory getDeviceHistory()
		{
			return { dy0_.deviceView(), dy1_.deviceView(),
			         dy2_.deviceView(), dy3_.deviceView() };
		}

		bool correctBoundaryOnly(
			real dt, realx3PointField_D& y, realx3PointField_D& dy)
		{
			return boundaryList_.correct(dt, y, dy);
		}

		bool correctPStructBoundaryOnly(
			real dt, pointStructure& pStruct, realx3PointField_D& vel)
		{
			return boundaryList_.correctPStruct(dt, pStruct, vel);
		}

		void resetHistory();

		void saveHistory();

		void restoreHistory();

};

} // pFlow

#endif
