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

#include "sphereDragImplicit.hpp"
#include "sphereDrag.hpp"
#include "SchillerNaumann.hpp"
#include "DiFelice.hpp"
#include "Rong.hpp"
#include "ErgunWenYu.hpp"
#include "Beetstra.hpp"
#include "noneDrag.hpp"
#include "CliftGraceWeber.hpp"

template class pFlow::coupling::sphereDrag<pFlow::coupling::SchillerNaumann>;

template class pFlow::coupling::sphereDragImplicit<pFlow::coupling::SchillerNaumann>;
template class pFlow::coupling::sphereDragImplicit<pFlow::coupling::DiFelice>;
template class pFlow::coupling::sphereDragImplicit<pFlow::coupling::ErgunWenYu>;
template class pFlow::coupling::sphereDragImplicit<pFlow::coupling::Rong>;
template class pFlow::coupling::sphereDragImplicit<pFlow::coupling::Beetstra>;
template class pFlow::coupling::sphereDragImplicit<pFlow::coupling::noneDrag>;
template class pFlow::coupling::sphereDragImplicit<pFlow::coupling::CliftGraceWeber>;
