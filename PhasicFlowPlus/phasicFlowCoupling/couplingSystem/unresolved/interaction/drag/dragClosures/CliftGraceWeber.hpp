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

// Clift, Grace & Weber (1978) Table 5.2

#ifndef __CliftGraceWeber_hpp__
#define __CliftGraceWeber_hpp__

#include "OFCompatibleHeader.hpp"
#include "typeInfo.hpp"

namespace pFlow::coupling
{

class CliftGraceWeber
{
	Foam::scalar 		residualRe_;

public:

	TypeInfoNV("CliftGraceWeber");

	CliftGraceWeber(const Foam::dictionary& dict);

	~CliftGraceWeber() = default;

	static inline
	Foam::scalar stokesFactor(Foam::scalar Re)
	{
		const Foam::scalar w = Foam::log10(Re);

		if(Re <= 0.01)
		{
			return 1.0 + (3.0/16.0)*Re/24.0;
		}
		if(Re <= 20.0)
		{
			return 1.0 + 0.1315*Foam::pow(Re, 0.82 - 0.05*w);
		}
		if(Re <= 260.0)
		{
			return 1.0 + 0.1935*Foam::pow(Re, 0.6305);
		}
		// 260 < Re <= 1500
		const Foam::scalar CD = Foam::pow(10.0, 1.6435 - 1.1242*w + 0.1558*w*w);
		return CD*Re/24.0;
	}

	inline
	Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep)const
	{
		const Foam::scalar Rec = Foam::max(Re, residualRe_);
		// Wen & Yu (1966)
		return stokesFactor(Rec)*Foam::pow(ep, -3.65);
	}

	inline
	Foam::scalar operator()(Foam::scalar Re, Foam::scalar ep)const
	{
		return dimlessDrag(Re, ep);
	}
};

}

#endif
