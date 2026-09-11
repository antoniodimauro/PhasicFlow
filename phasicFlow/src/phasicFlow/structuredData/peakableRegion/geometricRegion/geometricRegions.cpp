// Modified in 2026 by Antonio Di Mauro (T-junction CFD-DEM study);
// changes: git diff cf6d3fb main, github.com/antoniodimauro/PhasicFlow



#include "geometricRegions.hpp"


template class pFlow::geometricRegion<pFlow::box>;

template class pFlow::geometricRegion<pFlow::sphere>;
template class pFlow::geometricRegion<pFlow::pointGeom>;
