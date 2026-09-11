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

// Schiller & Naumann (1933) Z. Ver. Dtsch. Ing. 77, 318

#ifndef __SchillerNaumann_hpp__
#define __SchillerNaumann_hpp__

#include "OFCompatibleHeader.hpp"

#include "typeInfo.hpp"

namespace pFlow::coupling
{

class SchillerNaumann
{
	Foam::scalar 		residualRe_;

public:

	TypeInfoNV("SchillerNaumann");

	SchillerNaumann(const Foam::dictionary& dict);

	~SchillerNaumann() = default;

	inline
	Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep)const
	{
		const Foam::scalar Rec = Foam::max(Re, residualRe_);

		const Foam::scalar f = 1.0 + 0.15*Foam::pow(Rec, 0.687);

		return f*Foam::pow(ep, -3.65);
	}

	inline
	Foam::scalar operator()(Foam::scalar Re, Foam::scalar ep)const
	{
		return dimlessDrag(Re, ep);
	}
};

}

#endif
