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

#ifndef __costaCF_hpp__
#define __costaCF_hpp__

#include "types.hpp"
#include "symArrays.hpp"

namespace pFlow::cfModels
{

template<bool limited=true>
class costa
{
public:

	struct contactForceStorage
	{
		realx3 overlap_t_ = 0.0;
	};

	struct costaProperties
	{
		real 	en_ 	= 1000.0;

		real 	et_		= 1000.0;

		real    mu_ 	= 0.00001;

		INLINE_FUNCTION_HD
		costaProperties(){}

		INLINE_FUNCTION_HD
		costaProperties(real en, real et, real mu):
			en_(en), et_(et), mu_(mu)
		{}

		INLINE_FUNCTION_HD
		costaProperties(const costaProperties&)=default;

		INLINE_FUNCTION_HD
		costaProperties& operator=(const costaProperties&)=default;

		INLINE_FUNCTION_HD
		~costaProperties() = default;
	};

protected:

	using CostaArrayType = symArray<costaProperties>;

	int32 					numMaterial_ = 0;

	ViewType1D<real>  		rho_;

	CostaArrayType  		costaProperties_;

	real 					nSteps_ = 8;

	static constexpr real 	tangentialMassFactor_ = static_cast<real>(2.0/7.0);

	bool readCostaDictionary(const dictionary& dict)
	{
		auto en = dict.getVal<realVector>("en");
		auto et = dict.getVal<realVector>("et");
		auto mu = dict.getVal<realVector>("mu");

		nSteps_ = dict.getValOrSet<real>("collisionTimeSteps", static_cast<real>(8));

		if( nSteps_ < static_cast<real>(8) )
		{
			fatalErrorInFunction<<
			"collisionTimeSteps is "<<nSteps_<<
			", and 8 is the minimum. Costa et al. (2015) prescribe N = 8, and "
			"the predictor-corrector this code integrates the contact with "
			"does not merely lose accuracy below it, it ADDS energy. Driving a "
			"dry ballistic impact through the AdamsMoulton4 PEC scheme with a "
			"prescribed en of 0.97 returns, WITHOUT added mass:\n"
			"    N     3      4      5      6      7\n"
			"    e   2.583  0.805  2.538  1.462  1.317\n"
			"and with the added mass of a coupled run at rho_p/rho_f = 0.2:\n"
			"    e   1.752  1.164  1.010  0.999  0.991\n"
			"The behaviour is not even monotonic in N -- N = 4 happens to be "
			"dissipative uncoupled while N = 3 and N = 5 return more than twice "
			"the incoming energy -- so no value below 8 can be trusted on the "
			"strength of its neighbours. A contact that returns more than it "
			"receives will pump a dense cluster until the run fails, so this is "
			"refused rather than warned about.\n"
			"To take a larger timestep, raise dt and leave N at 8: the contact "
			"stays resolved and the overlap grows as dt^2.\n";
			return false;
		}

		if( nSteps_ < static_cast<real>(12) )
		{
			WARNING<<"collisionTimeSteps is "<<nSteps_<<
			". Below N = 12 the PEC predictor-corrector returns a restitution "
			"ABOVE the prescribed en for a particle with NO added mass (1.068 "
			"at N = 8, 1.462 at N = 6 for en = 0.97), so a free impact in an "
			"uncoupled run gains energy.\n"
			"In a coupled run the added mass stretches the contact and the "
			"realised restitution stays below one (0.988 at N = 8 for "
			"m_added/m_bare = 2.5), but it is then the added mass and not this "
			"dictionary that is setting it. Use N >= 12 unless the density "
			"ratio is known and fixed."<<END_WARNING;
		}

		auto nElem = en.size();

		if(nElem != et.size())
		{
			fatalErrorInFunction<<
			"sizes of en("<<nElem<<") and et("<<et.size()<<") do not match.\n";
			return false;
		}

		if(nElem != mu.size())
		{
			fatalErrorInFunction<<
			"sizes of en("<<nElem<<") and mu("<<mu.size()<<") do not match.\n";
			return false;
		}

		ForAll(i, en)
		{
			if( en[i] <= 0 || en[i] > 1 || et[i] <= 0 || et[i] > 1 )
			{
				fatalErrorInFunction<<
				"en and et must lie in (0,1]: the coefficients of eq (18) are "
				"built from their logarithms.\n";
				return false;
			}
		}

		uint32 nMat;
		if( !CostaArrayType::getN(nElem, nMat) )
		{
			fatalErrorInFunction<<
			"sizes of properties do not match a symetric array with size ("<<
			numMaterial_<<"x"<<numMaterial_<<").\n";
			return false;
		}
		else if( numMaterial_ != nMat)
		{
			fatalErrorInFunction<<
			"size mismatch for porperties. \n"<<
			"you supplied "<< numMaterial_<<" items in materials list and "<<
			nMat << " for other properties.\n";
			return false;
		}

		Vector<costaProperties> prop("prop", nElem);
		ForAll(i,en)
		{
			prop[i] = {en[i], et[i], mu[i]};
		}

		costaProperties_.assign(prop);

		return true;
	}

	static const char* modelName()
	{
		if constexpr (limited)
		{
			return "costaLimited";
		}
		else
		{
			return "costaNonLimited";
		}
		return "";
	}

public:

	TypeInfoNV(modelName());

	INLINE_FUNCTION_HD
	costa(){}

	costa(int32 nMaterial, const ViewType1D<real>& rho, const dictionary& dict)
	:
		numMaterial_(nMaterial),
		rho_("rho",nMaterial),
		costaProperties_("costaProperties",nMaterial)
	{
		Kokkos::deep_copy(rho_,rho);
		if(!readCostaDictionary(dict))
		{
			fatalExit;
		}
	}

	INLINE_FUNCTION_HD
	costa(const costa&) = default;

	INLINE_FUNCTION_HD
	costa(costa&&) = default;

	INLINE_FUNCTION_HD
	costa& operator=(const costa&) = default;

	INLINE_FUNCTION_HD
	costa& operator=(costa&&) = default;

	INLINE_FUNCTION_HD
	~costa()=default;

	INLINE_FUNCTION_HD
	int32 numMaterial()const
	{
		return numMaterial_;
	}

	INLINE_FUNCTION_HD
	real collisionTimeSteps()const
	{
		return nSteps_;
	}

	INLINE_FUNCTION_HD
	void contactForce
	(
		const real dt,
		const uint32 i,
		const uint32 j,
		const uint32 propId_i,
		const uint32 propId_j,
		const real Ri,
		const real Rj,
		const real ovrlp_n,
		const realx3& Vr,
		const realx3& Nij,
		contactForceStorage& history,
		realx3& FCn,
		realx3& FCt
	)const
	{
		auto prop = costaProperties_(propId_i,propId_j);

		real vrn = dot(Vr, Nij);
		realx3 Vt = Vr - vrn*Nij;

		history.overlap_t_ += Vt*dt;

		const real c43pi = static_cast<real>(4.0/3.0)*Pi;

		real mi = c43pi*pow(Ri,static_cast<real>(3.0))*rho_[propId_i];
		real mj = c43pi*pow(Rj,static_cast<real>(3.0))*rho_[propId_j];

		real me = (mi*mj)/(mi+mj);

		real Tn = nSteps_*dt;

		real lnEn = log(prop.en_);
		real lnEt = log(prop.et_);

		real kn   = me*(Pi*Pi + lnEn*lnEn)/(Tn*Tn);
		real etan = static_cast<real>(-2.0)*me*lnEn/Tn;

		real met  = tangentialMassFactor_*me;
		real kt   = met*(Pi*Pi + lnEt*lnEt)/(Tn*Tn);
		real etat = static_cast<real>(-2.0)*met*lnEt/Tn;

		FCn = (-kn * ovrlp_n - etan * vrn)*Nij;
		FCt = -kt * history.overlap_t_ - etat*Vt;

		real ft = length(FCt);
		real ft_fric = prop.mu_ * length(FCn);

		if(ft > ft_fric)
		{
			if( length(history.overlap_t_) > static_cast<real>(0.0))
			{
				if constexpr (limited)
				{
					FCt *= (ft_fric/ft);
					history.overlap_t_ = - (FCt/kt);
				}
				else
				{
					FCt = (FCt/ft)*ft_fric;
				}
			}
			else
			{
				FCt = 0.0;
			}
		}
	}
};

}

#endif
