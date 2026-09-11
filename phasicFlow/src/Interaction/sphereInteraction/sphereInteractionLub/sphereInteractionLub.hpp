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

#ifndef __sphereInteractionLub_hpp__
#define __sphereInteractionLub_hpp__

#include "sphereInteraction.hpp"
#include "lubrication.hpp"
#include "sphereInteractionLubKernels.hpp"

namespace pFlow
{

template<
	typename contactForceModel,
	typename geometryMotionModel,
	template <class, class, class> class contactListType>
class sphereInteractionLub
:
	public sphereInteraction<contactForceModel, geometryMotionModel, contactListType>
{
public:

	using Base = sphereInteraction<contactForceModel, geometryMotionModel, contactListType>;

	using ContactForceModel = typename Base::ContactForceModel;

	using MotionModel       = typename Base::MotionModel;

	using ContactListType   = typename Base::ContactListType;

	using rpPPInteraction   = typename Base::rpPPInteraction;

private:

	lubricationModels::lubrication      lubrication_;

	deviceViewType1D<real>              lubScale_;

	void reportLubrication()const;

	void updateLubricationScale(uint32 numPoints, uint32 lastItem);

	deviceViewType1D<real> addedMassView()const;

protected:

	bool sphereSphereInteraction() override;

public:

	TypeInfoTemplate13("sphereInteractionLub", ContactForceModel, MotionModel, ContactListType);

	sphereInteractionLub(
		systemControl& control,
		const particles& prtcl,
		const geometry& geom);

	add_vCtor
	(
		interaction,
		sphereInteractionLub,
		systemControl
	);
};

}

#include "sphereInteractionLub.cpp"

#endif
