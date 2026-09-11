# Modified in 2026 by Antonio Di Mauro (T-junction CFD-DEM study);
# changes: git diff cf6d3fb main, github.com/antoniodimauro/PhasicFlow

set(SourceFiles 
domainDistribute/domainDistribute.cpp
DEMSystem/DEMSystem.cpp
sphereDEMSystem/sphereFluidParticles.cpp
sphereDEMSystem/sphereDEMSystem.cpp
sphereDEMSystem/sphereFluidParticlesImplicit.cpp
sphereDEMSystem/sphereDEMSystemSemiImplicit.cpp
grainDEMSystem/grainFluidParticles.cpp
grainDEMSystem/grainDEMSystem.cpp

sphereDEMSystem/thermalSphereDEMSystem.cpp
sphereDEMSystem/thermalSphereFluidParticles.cpp

)

set(link_libs Kokkos::kokkos phasicFlow Particles Geometry Property Interaction Interaction Utilities)

pFlow_add_library_install(DEMSystems SourceFiles link_libs)



