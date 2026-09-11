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


template<typename cFM,typename gMM,template <class, class, class> class cLT>
pFlow::sphereInteractionLub<cFM,gMM, cLT>::sphereInteractionLub
(
	systemControl& control,
	const particles& prtcl,
	const geometry& geom
)
:
	Base(control, prtcl, geom)
{
	if( !lubrication_.read(*this) )
	{
		fatalErrorInFunction<<
		"could not read the lubrication settings.\n";
		fatalExit;
	}

	reportLubrication();

	if( lubrication_.active() )
	{
		bool anyBoundaryPairs = false;

		for(uint32 i=0; i<6; ++i)
		{
			if(this->activeBoundaries_[i]) anyBoundaryPairs = true;
		}

		if( anyBoundaryPairs )
		{
			fatalErrorInFunction<<
			"the lubrication model is active and at least one boundary "
			"allocates particle-particle pairs (a processor or periodic "
			"boundary).\n"
			"The boundary sphere-sphere kernel applies the contact model only: "
			"a pair split across that boundary would feel no squeeze film, and "
			"its coefficient would be missing from the per-particle stability "
			"bound. Run the DEM undecomposed and without periodic boundaries, "
			"or switch the lubrication model off.\n";
			fatalExit;
		}
	}
}

template<typename cFM,typename gMM,template <class, class, class> class cLT>
pFlow::deviceViewType1D<pFlow::real>
pFlow::sphereInteractionLub<cFM,gMM, cLT>::addedMassView()const
{
	if( this->time().lookupObjectName("addedMass") )
	{
		return this->time().template lookupObject<realPointField_D>("addedMass").deviceViewAll();
	}
	return deviceViewType1D<real>();
}

template<typename cFM,typename gMM,template <class, class, class> class cLT>
void pFlow::sphereInteractionLub<cFM,gMM, cLT>::reportLubrication()const
{
	if( !lubrication_.active() ) return;

	const auto  diams  = this->sphParticles_.spheres().boundingDiameter();
	const auto  masses = this->sphParticles_.spheres().mass();
	const real  dt     = this->dt();
	const real  sRatio = this->contactSearch_().sizeRatio();

	REPORT(1)<<"Lubrication (Costa et al. 2015 eq (34)) is "<<
	Yellow_Text("active")<<" for particle-particle pairs."<<END_REPORT;

	real needRatio = static_cast<real>(1.0);
	bool capped    = false;

	for(uint32 a=0; a<diams.size(); ++a)
	{
		for(uint32 b=a; b<diams.size(); ++b)
		{
			const real Ra = 0.5*diams[a];
			const real Rb = 0.5*diams[b];
			const real Rl = lubricationModels::lubrication::reducedRadius(Ra,Rb);

			const real mMin = min(masses[a], masses[b]);

			const real epsRes = lubrication_.minResolvedEps(Rl, mMin, dt);

			needRatio = max(needRatio,
				static_cast<real>(1.0) + lubrication_.range(Ra,Rb)/(Ra+Rb));

			REPORT(2)<<"pair ("<<a<<","<<b<<"): R_l = "<<Rl*1.0e6<<
			" um, smallest eps this dt can carry = "<<epsRes<<
			" (gap "<<epsRes*Rl*1.0e9<<" nm)"<<END_REPORT;

			if(epsRes > static_cast<real>(1.0e-2)) capped = true;
		}
	}

	if( sRatio < needRatio )
	{
		WARNING<<"contactSearch sizeRatio is "<<sRatio<<
		", but the lubrication range epsOut = "<<lubrication_.epsOut()<<
		" needs at least "<<needRatio<<
		" for the pair to be listed while the film is still acting.\n"<<
		"The film outside sizeRatio is simply missing: raise sizeRatio, or "
		"lower epsOut to match what the search reaches."<<END_WARNING;
	}

	if( lubrication_.cflLub() > static_cast<real>(1.0) )
	{
		WARNING<<"cflLub is "<<lubrication_.cflLub()<<
		", which asks the assembled lubrication damper for more than the "
		"0.158 negative-real-axis stability limit of the AdamsMoulton4 PEC "
		"scheme this solver integrates with (not the 1.285 of the same scheme "
		"in PECE mode, and not forward Euler's 2).\n"<<
		"The only thing left holding the run together is the added-mass margin "
		"m_bare/m_eff, which is a property of the density ratio and not of the "
		"model. Use cflLub <= 1."<<END_WARNING;
	}

	REPORT(2)<<"cap: sum_j c_ij <= "<<lubrication_.cflLub()<<
	" * 0.158 * m_eff_i/(2 dt), so the realised a*dt <= "<<
	lubrication_.cflLub()*static_cast<real>(0.158)<<
	" for any number of simultaneous films (AdamsMoulton4 PEC limit 0.158)."
	<<END_REPORT;

	REPORT(2)<<"the eps above are the BARE-mass, isolated-pair figures: "
	"m_added is not known until the first coupling call, and at run time the "
	"cap is formed on m_bare + m_added, so a coupled run resolves a finer "
	"film than this by m_eff/m_bare."<<END_REPORT;

	if( capped )
	{
		WARNING<<"the lubrication coefficient is capped by the timestep above "
		"eps = 1e-2 for at least one pair, so the timestep and not the surface "
		"roughness is setting the thinnest film that is resolved.\n"<<
		"The arrest itself survives this: integrating an approach with the cap "
		"in place still brings the pair to a dead stop without the surfaces "
		"meeting. What is lost is the shape of the last part of the approach, "
		"which is compressed into fewer steps. Reduce dt if the film thickness "
		"itself is the quantity of interest."
		<<END_WARNING;
	}
}

template<typename cFM,typename gMM,template <class, class, class> class cLT>
void pFlow::sphereInteractionLub<cFM,gMM, cLT>::updateLubricationScale
(
	uint32 numPoints,
	uint32 lastItem
)
{
	if( !lubrication_.active() ) return;

	if( lubScale_.size() < numPoints )
	{
		Kokkos::realloc(lubScale_, numPoints);
	}

	Kokkos::deep_copy(lubScale_, static_cast<real>(0));

	pFlow::sphereInteractionKernels::ppLubricationSumFunctor
		lubSum(
			this->ppContactList_(),
			this->sphParticles_.diameter().deviceViewAll(),
			this->sphParticles_.pointPosition().deviceViewAll(),
			lubScale_,
			lubrication_
			);

	Kokkos::parallel_for(
		"ppLubricationSum",
		rpPPInteraction(0,lastItem),
		lubSum
		);

	Kokkos::fence();

	auto massView = this->sphParticles_.mass().deviceViewAll();

	auto addedView = addedMassView();

	pFlow::sphereInteractionKernels::ppLubricationScaleFunctor
		<decltype(massView)>
		lubScale(this->dt(), massView, addedView, lubScale_, lubrication_);

	Kokkos::parallel_for(
		"ppLubricationScale",
		rpPPInteraction(0,numPoints),
		lubScale
		);

	Kokkos::fence();
}

template<typename cFM,typename gMM,template <class, class, class> class cLT>
bool pFlow::sphereInteractionLub<cFM,gMM, cLT>::sphereSphereInteraction()
{
	auto lastItem = this->ppContactList_().loopCount();

	updateLubricationScale(this->sphParticles_.diameter().deviceViewAll().size(), lastItem);

	// create the kernel functor 
	pFlow::sphereInteractionKernels::ppInteractionLubFunctor 
		ppInteraction(
			this->dt(),
			this->forceModel_(),
			this->ppContactList_(),
			this->sphParticles_.diameter().deviceViewAll(),
			this->sphParticles_.propertyId().deviceViewAll(),
			this->sphParticles_.pointPosition().deviceViewAll(),
			this->sphParticles_.velocity().deviceViewAll(),
			this->sphParticles_.rVelocity().deviceViewAll(),
			this->sphParticles_.contactForce().deviceViewAll(),
			this->sphParticles_.contactTorque().deviceViewAll(),
			lubScale_,
			lubrication_
			);
	
	Kokkos::parallel_for(
		"ppInteraction",
		rpPPInteraction(0,lastItem),
		ppInteraction
		);

	Kokkos::fence();	

	return true;
}
