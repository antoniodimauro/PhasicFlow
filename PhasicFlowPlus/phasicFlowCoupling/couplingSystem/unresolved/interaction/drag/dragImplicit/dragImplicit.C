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

#include "dragImplicit.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "fluidAcceleration.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"

pFlow::coupling::dragImplicit::dragImplicit
(
	const unresolvedCouplingSystem& uCS,
	const porosity& 				prsty
)
:
    drag(uCS, prsty),
    undisturbedSampler_(uCS, "undisturbedFlow"),
    quad_(prsty.cMesh()),
    spScratch_("dragSpScratch", real(0), uCS.centerMass()),
    upScratch_("dragUpScratch", realx3(0), uCS.centerMass()),
    fbScratch_("dragFbScratch", realx3(0), uCS.centerMass()),
    ueScratch_("dragUeScratch", realx3(0), uCS.centerMass()),
    wxScratch_("dragWxScratch", realx3(0), uCS.centerMass()),
    ubScratch_("dragUbScratch", realx3(0), uCS.centerMass()),
    basset_(uCS, drag::getDict(uCS)),
    pair_(uCS, drag::getDict(uCS))
{
	Foam::Info
		<< "    Undisturbed-flow sampling: "
		<< Green_Text(undisturbedSampler_.method())
		<< " (use distribution when d_p is not small vs. the cell size)\n";

	const Foam::word selfInduced = dict().getOrDefault<Foam::word>
	(
		"selfInducedVelocity",
		"analytic"
	);

	if(selfInduced == "analytic")
	{
		selfInduced_ = true;
		maxBeta_ = dict().getOrDefault<Foam::scalar>("maxSelfInducedFraction", 0.95);

		const Foam::word fvMethod =
			pFlow::coupling::fluidVelocityMethod(
				Porosity().uCS().unresolvedDict());

		selfInducedSampleFactor_ =
			pFlow::coupling::selfInducedSampleFactor(fvMethod);

		Foam::Info
			<< "    Self-induced velocity: "
			<< Green_Text("analytic")
			<< " (Ireland & Desjardins 2017 eq 36/39: the drag is\n"
			   "      evaluated at the undisturbed slip, and eq 33 restores the\n"
			   "      undisturbed void fraction)\n"
			<< "      fluidVelocity " << fvMethod << " -> the sample is "
			<< (selfInducedSampleFactor_ < 1.0 ?
			      "a kernel average, deficit 1/sqrt(2) of the blob-centre value" :
			      "the value at the particle, deficit = the blob-centre value")
			<< "\n";

		if(selfInducedSampleFactor_ < 1.0)
		{
			if(Plus::processor::isMaster())
			{
    			WARNING
    				<< "fluidVelocity is `distribution`, so the deficit present in "
    				<< "the sample is 1/sqrt(2) of the blob-centre value. That "
    				<< "factor is NOT stated by Ireland & Desjardins (2017) or by "
    				<< "Balachandar et al. (2019) -- both read the velocity at the "
    				<< "particle -- it is derived here from their regularised "
    				<< "Stokeslet, and is exact only in the Stokes limit. Use "
    				<< "fluidVelocity cell or interpolate if you would rather not "
    				<< "rely on it." << END_WARNING;
			}
		}
	}
	else if(selfInduced == "none")
	{
		selfInduced_ = false;
		Foam::Info
			<< "    Self-induced velocity: "
			<< Yellow_Text("none")
			<< " (the particle is left sitting in its own wake)\n";
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			fatalErrorInFunction
				<< "Unknown selfInducedVelocity \"" << selfInduced << "\" in "
				<< dict().name()
				<< "\nAvailable ones are: analytic, none"
				<< endl;
		}
		Plus::processor::abort(0);
	}

	faxen_ = dict().getOrDefault<Foam::Switch>("faxenCorrection", false);

	if(faxen_)
	{
		const Foam::word vMethod = fluidVelocityMethod(uCS.unresolvedDict());

		const Foam::word sMethod = fluidFieldSampler::readMethod(uCS);

		if(vMethod == "distribution" || sMethod == "distribution")
		{
			if(Plus::processor::isMaster())
			{
				fatalErrorInFunction
					<< "faxenCorrection yes requires fluidVelocity and "
					<< "fluidSampling to be cell or interpolate, not "
					<< "distribution (fluidVelocity = " << vMethod
					<< ", fluidSampling = " << sMethod << ").\n"
					<< "A kernel average already applies a finite-size "
					<< "correction of (sigma^2/2)*grad^2 u. The Faxen average "
					<< "is the finite-size correction, evaluated exactly over "
					<< "the particle, so a kernel-averaged sample would count "
					<< "it twice."
					<< endl;
			}
			Plus::processor::abort(0);
		}

		Foam::Info
			<< "    Faxen correction: "
			<< Green_Text("on")
			<< " (exact Faxen averages over the particle:"
			   " surface for drag, volume for added mass)\n";
	}

	conservativeReaction_ = dict().getOrDefault<Foam::Switch>
	(
		"conservativeReaction", false
	);

	selfInducedLift_ = dict().getOrDefault<Foam::Switch>
	(
		"selfInducedLift", false
	);

	pointForceReaction_ = dict().getOrDefault<Foam::Switch>
	(
		"pointForceReaction", false
	);

	if(pointForceReaction_)
	{
		Foam::Info
			<< "    Drag reaction shape: " << Green_Text("point force")
			<< " (Maxey & Patel 2001). The implicit Sp*(U - u_p) is\n"
			   "      kept and Sp*(U - Ubar) is removed by deferred"
			   " correction, so the fluid receives K(x)*F\n"
			   "      without the force dipole Sp*grad(U) that damps the"
			   " local shear and rotation.\n";
	}

	{
		const Foam::word rotDrag = dict().getOrDefault<Foam::word>
		(
			"rotationDrag", "none"
		);

		if(rotDrag == "Gotoh1990")
		{
			rotationDrag_ = true;
			rotVorticitySampler_ =
				makeUnique<fluidFieldSampler>(uCS, "rotationDragVorticity");
			rotStrainSampler_ =
				makeUnique<fluidFieldSampler>(uCS, "rotationDragStrain");

			Foam::Info
				<< "    Rotation drag: " << Green_Text("Gotoh1990")
				<< " (Gotoh 1990; Candelier, Mehlig & Magnaudet 2019 eq. 4.10).\n"
				   "      Rotation part of the local flow Omega = (|omega| - sigma)/2,"
				   " sigma = sqrt(2 S:S) (the split\n"
				   "      of the rotation lift); Ta = Omega a^2/nu. Drag factor"
				   " max(f(Re), 1 + 0.524 Ta^1/2):\n"
				   "      the larger of the Oseen and the rotational screening"
				   " (exact for Re_p << Ta^1/2 << 1).\n"
				   "      The self-induced velocity is screened the same way at"
				   " a_s = sqrt(pi)*sigma_kernel:\n"
				   "      min(psiOseen, 1/(1 + 0.524 Ta_s^1/2)). Pure shear:"
				   " no change.\n";

			const Foam::word comb = dict().getOrDefault<Foam::word>
			(
				"rotationDragCombination", "max"
			);

			if(comb == "quadrature")
			{
				rotationDragQuadrature_ = true;

				Foam::Info
					<< "      Combination: " << Green_Text("quadrature")
					<< ", 1 + sqrt((f - 1)^2 + (g - 1)^2) in place of max(f, g):"
					   " inertia and rotation both\n"
					   "      thin the boundary layer (Davis 1992, Table 2).\n";
			}
			else if(comb != "max")
			{
				FatalErrorInFunction
					<< "rotationDragCombination " << comb << " unknown;"
					<< " valid: max, quadrature" << Foam::exit(Foam::FatalError);
			}

			const Foam::word geo = dict().getOrDefault<Foam::word>
			(
				"rotationDragGeostrophic", "none"
			);

			if(geo == "Stewartson1953")
			{
				rotationDragGeostrophic_ = true;

				Foam::Info
					<< "      Geostrophic drag: " << Green_Text("Stewartson1953")
					<< ", g += Ta/(1 + Ta) (0.49*4/9) Ta (Stewartson 1953;"
					   " Mason 1975),\n"
					   "      weight Omega/(Omega + nu/a^2); the screening of the"
					   " self-induced velocity is unchanged.\n";
			}
			else if(geo != "none")
			{
				FatalErrorInFunction
					<< "rotationDragGeostrophic " << geo << " unknown;"
					<< " valid: none, Stewartson1953" << Foam::exit(Foam::FatalError);
			}

			rotationSelfInducedLateral_ = dict().getOrDefault<Foam::Switch>
			(
				"rotationSelfInducedLateral", Foam::Switch(false)
			);

			if(rotationSelfInducedLateral_)
			{
				Foam::Info
					<< "      Self-induced velocity, lateral part: " << Green_Text("on")
					<< ", beta_perp = beta 0.0517 T/(1 + 0.524 T), T = Ta_s^1/2, along\n"
					   "      omega_hat x slip, removed from the sampled slip (Gotoh 1990;"
					   " Candelier, Mehlig &\n      Magnaudet 2019, eq. 4.10).\n";
			}
		}
		else if(rotDrag != "none")
		{
			FatalErrorInFunction
				<< "rotationDrag " << rotDrag << " unknown; valid: none, Gotoh1990"
				<< Foam::exit(Foam::FatalError);
		}
	}

	if(selfInducedLift_)
	{
		Foam::Info
			<< "    Self-induced velocity: " << Green_Text("drag and lift")
			<< " fed back (Esmaily & Horwitz 2018). The drag sees\n"
			   "      u_eff + beta/sp0*F_lift; the lift uses the undisturbed"
			   " slip amp*(u_eff - u_p).\n";
	}

	if(conservativeReaction_)
	{
		Foam::Info
			<< "    Drag reaction on the fluid: " << Green_Text("conservative")
			<< ". The implicit source Sp*(U - u_p) integrates to\n"
			   "      sp*(Ubar - u_p), Ubar the kernel average of U; the particle"
			   " gets\n      sp*(u_eff - u_p) plus the wall terms. The difference"
			   " sp*(Ubar - u_eff) - wallExtra\n      is added explicitly over"
			   " the same kernel, so the fluid receives exactly\n      minus"
			   " the force on the particle.\n";
	}

	wallCorr_ = dict().getOrDefault<Foam::Switch>("wallCorrection", false);

	if(wallCorr_)
	{
		minWallRadii_ = dict().getOrDefault<Foam::scalar>("minWallRadii", 2.0);

		minBlobWallRadii_ = dict().getOrDefault<Foam::scalar>
		(
			"minBlobWallRadii", minWallRadii_
		);

		const Foam::word wModel = dict().getOrDefault<Foam::word>
		(
			"wallDragModel", "Stokes"
		);

		if(wModel == "finiteRe")
		{
			wallFiniteRe_ = true;
		}
		else if(wModel == "gapReynolds")
		{
			wallFiniteRe_ = true;
			wallGapRe_    = true;
			wallGapReScale_ = dict().getOrDefault<Foam::scalar>("wallGapReScale", 1.25);
			Foam::Info
				<< "    Near-wall normal factor: " << Yellow_Text("gapReynolds")
				<< ", 1 + (lambda_Brenner - 1)/(1 + Re_h/c), c = " << wallGapReScale_ << ",\n"
				   "      Re_h = |u_rel| h/nu on the gap h (after Vasseur & Cox 1977);\n"
				   "      parallel factor as finiteRe.\n";
		}
		else if(wModel != "Stokes")
		{
			fatalErrorInFunction
				<< "wallDragModel must be Stokes, finiteRe or gapReynolds, got "
				<< wModel << " in " << dict().name() << endl;
			Plus::processor::abort(0);
		}

		wallNormalImplicit_ = dict().getOrDefault<Foam::Switch>
		(
			"wallNormalImplicit", false
		);

		if(wallNormalImplicit_)
		{
			Foam::Info
				<< "    Near-wall normal resistance: " << Green_Text("implicit")
				<< ", in the drag coefficient weighted by cos^2 of the\n"
				   "      slip to the wall normal (exact for normal and parallel"
				   " slip). Wall distance taken\n      at the particle centre."
				   " minWallRadii can go down to ~1.01.\n";
		}

		Foam::Info
			<< "    Near-wall drag: " << Green_Text("on")
			<< " (Brenner 1961 normal; Goldman, Cox & Brenner 1967\n"
			   "      I and II for translation, rotation and shear parallel to"
			   " the wall.\n      Frozen below h/a = "
			<< minWallRadii_ << ", where none of them covers the lubrication"
			   " limit.\n      The self-induction wall factor is frozen"
			   " separately, below h/a_eff = "
			<< minBlobWallRadii_ << ",\n      because a_eff = sqrt(pi/2)*sigma"
			   " is a different length from a.\n      The wall-induced lift is"
			   " not part of the drag: the lift class adds the\n"
			   "      translation-induced lift of Zeng et al. (2009) when"
			   " wallCorrection is on.)\n";

		if(wallFiniteRe_)
		{
			Foam::Info
				<< "    Near-wall drag model: " << Green_Text("finiteRe") << ". Normal:"
				   " Ryu & Matsudaira, Chem. Eng. Sci. 65 (2010)\n"
				   "      4913 eq (6), Brenner's factor reduced by"
				   " (1 + sqrt(Re))^(0.39 (H/R)^-0.85),\n"
				   "      fitted for Re <= 40, 1.4 <= H/R <= 9. Parallel:"
				   " Zeng et al., Phys. Fluids\n"
				   "      21 (2009) 033302 eqs (16)-(18), on the RELATIVE"
				   " parallel velocity (their shear\n      and translation"
				   " fits do not superpose at finite Re).\n"
				   "      Both relative to the unbounded drag 1 + 0.15 Re^0.687."
				   " The rotation-\n      induced force and all wall torques"
				   " keep their Stokes factors.\n";
		}
	}

	const Foam::word form = dict().getOrDefault<Foam::word>
	(
		"undisturbedFlowForce",
		"pressureGradient"
	);

	if(form == "materialDerivative")
	{
		undisturbedIsMaterialDerivative_ = true;
		Foam::Info
			<< "    Undisturbed-flow force: "
			<< Green_Text("materialDerivative")
			<< " (F = rho_f*V_p*Du/Dt; Maxey-Riley form)\n";
	}
	else if(form == "pressureGradient")
	{
		undisturbedIsMaterialDerivative_ = false;
		Foam::Info
			<< "    Undisturbed-flow force: "
			<< Green_Text("pressureGradient")
			<< " (F = -V_p*grad(p))\n";
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			fatalErrorInFunction
				<< "Unknown undisturbedFlowForce \"" << form << "\" in "
				<< dict().name()
				<< "\nAvailable ones are: pressureGradient, materialDerivative"
				<< endl;
		}
		Plus::processor::abort(0);
	}
}

void pFlow::coupling::dragImplicit::calculateDragForce
(
    const fluidAveraging&,
    const solidAveraging&,
    const Plus::realProcCMField&,
    const distributionBase&,
    Plus::realx3ProcCMField&
)
{
    fatalErrorInFunction
        << "an implicit drag model was selected outside the momentumSemiImplicit "
           "coupling system; use the explicit sphereDrag models with the "
           "momentum coupling system."
        << endl;
    Plus::processor::abort(0);
}

pFlow::uniquePtr<pFlow::coupling::dragImplicit> pFlow::coupling::dragImplicit::create
(
	const unresolvedCouplingSystem& uCS,
	const porosity& 				prsty
)
{
	const auto& dragDict = drag::getDict(uCS);
	auto shapeName = uCS.shapeTypeName();
	auto dType = lookupDict<Foam::word>(dragDict, "model");
	Foam::word dragType = angleBracketsNames(shapeName+"DragImplicit", dType);

	if( drag::couplingSystemvCtorSelector_.search(dragType))
	{
		Foam::Info<<"    Creating drag force "<<Green_Text(dragType)<<" ...\n\n";
		auto basePtr = drag::couplingSystemvCtorSelector_[dragType] (uCS, prsty);
		auto* implicitPtr = dynamic_cast<dragImplicit*>(basePtr.release());
		return uniquePtr<dragImplicit>(implicitPtr);
	}

	if(Plus::processor::isMaster())
	{
		printKeys
		(
			fatalError << "Ctor Selector "<< dragType << " does not exist. \n"
			<<"Available ones are: \n\n"
			,
			drag::couplingSystemvCtorSelector_
		);
	}
	Plus::processor::abort(0);
	return nullptr;
}

void pFlow::coupling::dragImplicit::commitPairForces
(
	const Plus::realx3ProcCMField& lift
)const
{
	if(!pair_.active()) return;

	const auto& parCellInd = this->parCellIndex();

	for(size_t i=0; i<parCellInd.size() && i<lift.size(); ++i)
	{
		if(parCellInd[i] < 0) continue;

		pairModel().addForce(i, Foam::vector(lift[i].x(), lift[i].y(), lift[i].z()));
	}

	pairModel().commit();
}

void pFlow::coupling::dragImplicit::addBassetAddedMass
(
	Plus::realProcCMField& addedMass
)const
{
	if(!basset_.active()) return;

	const auto& parCellInd = this->parCellIndex();
	const size_t nPar = parCellInd.size();

	#pragma omp parallel for schedule(static)
	for(size_t i=0; i<nPar; ++i)
	{
		if(parCellInd[i] < 0) continue;

		addedMass[i] += static_cast<real>(basset_.addedMass(i));
	}
}

void pFlow::coupling::dragImplicit::reportDisturbanceRegime
(
	Foam::scalar dp,
	Foam::scalar sigma
)const
{
	if(reportedRegime_ || sigma <= 0 || dp <= 0) return;
	reportedRegime_ = true;

	const Foam::scalar dfil = 2.0*Foam::sqrt(2.0*Foam::log(2.0))*sigma;
	const Foam::scalar fac  = (dfil > 0) ? (dp*dp)/(18.0*dfil*dfil) : 0.0;

	Foam::Info
		<< Blue_Text("Disturbance regime: ")
		<< "delta_f/d_p = " << Yellow_Text(dfil/dp)
		<< " (Ireland & Desjardins 2017 validate 2 to 16),\n"
		<< "      St_delta = (rho_p/rho_f) * " << Yellow_Text(fac)
		<< "  (their eq 42; validated 1.1 to 17).\n"
		<< "      Below that range their correction over-corrects -- see"
		   " particleGaussian.hpp.\n"
		<< Foam::endl;
}

void pFlow::coupling::dragImplicit::reportWall
(
	Foam::scalar factorObserved,
	Foam::scalar minHoverA
)const
{
	if(!wallCorr_) return;

	Foam::reduce(factorObserved, Foam::maxOp<Foam::scalar>());
	Foam::reduce(minHoverA,      Foam::minOp<Foam::scalar>());

	if(factorObserved > reportedWall_ + 0.05)
	{
		reportedWall_ = factorObserved;

		Foam::Info
			<< Blue_Text("Near-wall drag: ")
			<< "largest factor " << Yellow_Text(factorObserved)
			<< " at h/a = " << Yellow_Text(minHoverA) << Foam::endl;

		if(factorObserved > 1.2 && !warnedUnsteadyWall_)
		{
			warnedUnsteadyWall_ = true;

			if(Plus::processor::isMaster())
			{
    			WARNING
    				<< "The near-wall factors have reached " << factorObserved
    				<< ", so the wall is a leading-order influence on the steady "
    				<< "force. The history force, the rotational history couple, "
    				<< "the added mass and the added inertia are still using their "
    				<< "UNBOUNDED kernels: the wall corrections implemented here "
    				<< "are steady-drag factors (" << (wallGapRe_ ? "gapReynolds" : wallFiniteRe_ ? "finiteRe" : "Stokes")
    				<< ") and say nothing about the unsteady problem."
    				<< END_WARNING;
			}
		}

		if(minHoverA < minWallRadii_ && !warnedLubrication_)
		{
			warnedLubrication_ = true;

			if(Plus::processor::isMaster())
			{
    			WARNING
    				<< "A particle is within h/a = " << minHoverA
    				<< " of a wall. Brenner (1961) eq (2.19) is exact but its "
    				<< "series converges ever more slowly there, and the parallel "
    				<< "resistance of Goldman, Cox & Brenner (1967) diverges "
    				<< "logarithmically at contact, which no table can carry. The "
    				<< "factors are held at their h/a = " << minWallRadii_
    				<< " values; the lubrication regime is NOT modelled."
    				<< END_WARNING;
			}
		}
	}
}

void pFlow::coupling::dragImplicit::reportWallRe
(
	Foam::scalar ReMaxNearWall
)const
{
	if(!wallCorr_ || !wallFiniteRe_ || warnedWallRe_) return;

	Foam::reduce(ReMaxNearWall, Foam::maxOp<Foam::scalar>());

	if(ReMaxNearWall > 40.0)
	{
		warnedWallRe_ = true;

		if(Plus::processor::isMaster())
		{
			WARNING
				<< "A particle within H/R = 9 of a wall has Re = " << ReMaxNearWall
				<< ". The normal wall factor of Ryu & Matsudaira (2010) is fitted"
				<< " for Re <= 40 and is extrapolated beyond it." << END_WARNING;
		}
	}
}

void pFlow::coupling::dragImplicit::reportConservative
(
	Foam::scalar corrOverDrag
)const
{
	if(!conservativeReaction_ || reportedConservative_) return;

	Foam::reduce(corrOverDrag, Foam::maxOp<Foam::scalar>());

	reportedConservative_ = true;

	Foam::Info
		<< Blue_Text("Drag reaction: ")
		<< "first conservative correction, largest |correction|/|drag| = "
		<< Yellow_Text(corrOverDrag) << Foam::endl;
}

void pFlow::coupling::dragImplicit::reportWallSelfInduction
(
	Foam::scalar gObserved,
	Foam::scalar minHOverAEff
)const
{
	if(!wallCorr_ || !selfInduced_) return;

	Foam::reduce(gObserved,    Foam::minOp<Foam::scalar>());
	Foam::reduce(minHOverAEff, Foam::minOp<Foam::scalar>());

	if(gObserved < reportedWallBeta_ - 0.05)
	{
		reportedWallBeta_ = gObserved;

		Foam::Info
			<< Blue_Text("Wall self-induction: ")
			<< "the wall has cut the self-induced velocity to "
			<< Yellow_Text(gObserved)
			<< " of its unbounded value, at h/a_eff = "
			<< Yellow_Text(minHOverAEff) << Foam::endl;

		if(minHOverAEff < minBlobWallRadii_ && !warnedBlobWall_)
		{
			warnedBlobWall_ = true;

			if(Plus::processor::isMaster())
			{
    			WARNING
    				<< "The deposition kernel is within h/a_eff = "
    				<< minHOverAEff << " of a wall, where a_eff = sqrt(pi/2)*sigma "
    				<< "is the blob's Stokes radius. Gamma_k of Pakseresht, "
    				<< "Esmaily & Apte (2020) is the wall drag factor of that "
    				<< "blob, and those factors are not defined once it would "
    				<< "intersect the wall; the value is frozen at h/a_eff = "
    				<< minBlobWallRadii_ << ". The volume-filtered picture is in "
    				<< "trouble of its own there too -- the source is being "
    				<< "spread over a support the wall cuts into. Narrow the "
    				<< "kernel, or read the near-wall drag as an estimate."
    				<< END_WARNING;
			}
		}
	}
}

void pFlow::coupling::dragImplicit::reportFaxen(Foam::scalar ratioObserved)const
{
	if(!faxen_) return;

	Foam::reduce(ratioObserved, Foam::maxOp<Foam::scalar>());

	if(ratioObserved > reportedFaxen_ + 0.05)
	{
		reportedFaxen_ = ratioObserved;

		Foam::Info
			<< Blue_Text("Faxen correction: ")
			<< "surface mean departs from the point value by "
			<< Yellow_Text(100.0*ratioObserved)
			<< Yellow_Text(" % at a particle")
			<< Foam::endl;

		if(ratioObserved > 1.0)
		{
			if(Plus::processor::isMaster())
			{
    			WARNING
    				<< "The Faxen surface mean now differs from the point value by "
    				<< "(" << 100.0*ratioObserved << " %). The leading-order "
    				<< "finite-size expansion is being used outside its range: the "
    				<< "particle spans more of the velocity profile than a "
    				<< "quadratic fit to it describes. Reduce d_p, or accept that "
    				<< "the point-particle closure is at its limit here."
    				<< END_WARNING;
			}
		}
	}
}

void pFlow::coupling::dragImplicit::reportRotationDrag
(
	Foam::scalar factorObserved,
	Foam::scalar taObserved,
	Foam::scalar psiRatioObserved,
	Foam::scalar lateralObserved
)const
{
	if(!rotationDrag_) return;

	Foam::reduce(factorObserved, Foam::maxOp<Foam::scalar>());
	Foam::reduce(taObserved, Foam::maxOp<Foam::scalar>());
	Foam::reduce(psiRatioObserved, Foam::minOp<Foam::scalar>());
	Foam::reduce(lateralObserved, Foam::maxOp<Foam::scalar>());

	if(factorObserved > reportedRotationFactor_ + 0.02 || reportedRotationFactor_ == 0)
	{
		reportedRotationFactor_ = Foam::max(factorObserved, Foam::scalar(1e-12));

		Foam::Info
			<< "Rotation drag: largest factor " << factorObserved
			<< " on the standard drag, largest Ta " << taObserved
			<< "; self-induced velocity screened to at least "
			<< psiRatioObserved << " of its non-rotating value";
		if(rotationSelfInducedLateral_)
		{
			Foam::Info << "; lateral part up to " << lateralObserved << " of the slip";
		}
		Foam::Info << Foam::endl;
	}
}

void pFlow::coupling::dragImplicit::reportSelfInduction
(
	Foam::scalar betaObserved,
	Foam::label  anyWidth
)const
{
	if(!selfInduced_) return;

	Foam::reduce(anyWidth, Foam::maxOp<Foam::label>());
	Foam::reduce(betaObserved, Foam::maxOp<Foam::scalar>());

	if(anyWidth == 0)
	{
		if(!warnedNoWidth_)
		{
			warnedNoWidth_ = true;
			if(Plus::processor::isMaster())
			{
    			WARNING
    				<< "The momentum source is being deposited by a method with no "
    				<< "defined kernel width, so the self-induced velocity cannot "
    				<< "be divided out and the drag is under-predicted. Set "
    				<< "momentumExchange to distribution with a Gaussian-family "
    				<< "distributionMethod, or accept it with "
    				<< "selfInducedVelocity none."
    				<< END_WARNING;
			}
		}
		return;
	}

	if(betaObserved > reportedBeta_ + 0.05)
	{
		reportedBeta_ = betaObserved;

		if(betaObserved >= maxBeta_)
		{
			if(Plus::processor::isMaster())
			{
    			WARNING
    				<< "Self-induced velocity has reached the cap: beta = "
    				<< betaObserved << " >= " << maxBeta_ << ". The particle is "
    				<< "essentially moving in its own wake and the correction is "
    				<< "being clamped -- the drag is no longer trustworthy. "
    				<< "particleGaussian's minFilterPerDiameter keeps the kernel "
    				<< "several diameters wide, where beta is well below the cap, "
    				<< "so reaching it means either that floor has been lowered or "
    				<< "that the deposit is being made by a method whose reported "
    				<< "width is narrower than the one it actually uses."
    				<< END_WARNING;
			}
		}
		else
		{
			Foam::Info
				<< "    Self-induced velocity: beta_max = " << betaObserved
				<< ", slip corrected by " << 1.0/(1.0-betaObserved) << "x\n";
		}
	}
}

Foam::tmp<Foam::volVectorField>
pFlow::coupling::dragImplicit::pressureGradient(const Foam::volScalarField& rho)const
{
	if(undisturbedIsMaterialDerivative_)
	{
		const auto& U =
			Porosity().mesh().lookupObject<Foam::volVectorField>("U");

		return -rho*materialDerivativeU(U);
	}

	return drag::pressureGradient(rho);
}

const Foam::volVectorField& pFlow::coupling::dragImplicit::materialDerivativeU
(
	const Foam::volVectorField& U
)const
{
	if(!tDDtU_.valid())
	{
		tDDtU_ = materialDerivative(U);
	}

	return tDDtU_();
}

void pFlow::coupling::dragImplicit::setLiftFeedback
(
	const Plus::realx3ProcCMField& lift,
	const Plus::centerMassField&   centres
)const
{
	if(!selfInducedLift_) return;

	const size_t n = lift.size();

	liftFeedback_.resize(n);
	liftFeedbackPos_.resize(n);

	for(size_t i=0; i<n; ++i)
	{
		liftFeedback_[i]    = Foam::vector(lift[i].x(), lift[i].y(), lift[i].z());
		liftFeedbackPos_[i] =
			Foam::vector(centres[i].x(), centres[i].y(), centres[i].z());
	}
}

bool pFlow::coupling::dragImplicit::liftFeedbackFor
(
	size_t              i,
	const Foam::vector& centre,
	Foam::scalar        radius,
	Foam::vector&       fL
)const
{
	if(!selfInducedLift_ || i >= liftFeedback_.size()) return false;

	if(Foam::mag(centre - liftFeedbackPos_[i]) > radius) return false;

	fL = liftFeedback_[i];
	return true;
}

bool pFlow::coupling::dragImplicit::undisturbedVelocity(size_t i, Foam::vector& u)const
{
	if(!selfInducedLift_ || i >= undisturbedSet_.size() || !undisturbedSet_[i])
	{
		return false;
	}

	u = undisturbedVel_[i];
	return true;
}
