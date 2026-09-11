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

#include "CliftGraceWeber.hpp"

pFlow::coupling::CliftGraceWeber::CliftGraceWeber(const Foam::dictionary& dict)
:
	residualRe_(lookupDict<Foam::scalar>(dict, "residualRe"))
{
}
