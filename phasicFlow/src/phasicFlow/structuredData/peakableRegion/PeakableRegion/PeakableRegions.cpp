/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow       
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
  Modified in 2026 by Antonio Di Mauro (T-junction CFD-DEM study);
  changes: git diff cf6d3fb main, github.com/antoniodimauro/PhasicFlow
------------------------------------------------------------------------------  
Licence:
  This file is part of phasicFlow code. It is a free software for simulating 
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions. 
 
  phasicFlow is distributed to help others in their research in the field of 
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/


#include "PeakableRegion.hpp"
#include "geometricRegions.hpp"

template class pFlow::PeakableRegion<pFlow::box>;

template class pFlow::PeakableRegion<pFlow::sphere>;

template class pFlow::PeakableRegion<pFlow::cylinder>;
template class pFlow::PeakableRegion<pFlow::pointGeom>;
