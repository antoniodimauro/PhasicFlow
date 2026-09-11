# PhasicFlow

CFD-DEM stack for a T-junction study: the DEM library `phasicFlow` and the
OpenFOAM coupling library `PhasicFlowPlus` in one tree, with the unresolved,
four-way coupled solver `unresolvedSphereSemiImplicitPFPlus` for spheres.

## Origin

`phasicFlow/` and `PhasicFlowPlus/` are the projects of www.cemf.ir (GPL-3),
taken at phasicFlow `454ba4d2` and PhasicFlowPlus `bbf2110`. The first commit
is the two trees as distributed, the second the modifications
(`git diff cf6d3fb main`). Added files: Copyright (C) Antonio Di Mauro, GPL-3.
Kokkos 4.4.01 is a git submodule (`phasicFlow/thirdParty/kokkos`).

## Build

OpenFOAM v2406 or later, CMake 3.16 or later, GCC 10.1 or later.

```bash
source <openfoam>/etc/bashrc
./Allwmake -j 8
```

Everything installs into `$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS`.

## Run

```bash
./init
mpirun --bind-to none -np 12 unresolvedSphereSemiImplicitPFPlus -parallel
```

`settings/settingsDict` needs `integrationMethod AdamsMoulton4PEC;` and
`libs ("libpFCouplingUtilities.so");`. `--bind-to none` lets the DEM use its
OpenMP threads.
