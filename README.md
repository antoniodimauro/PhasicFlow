# PhasicFlow

A modified version of **phasicFlow** and **PhasicFlowPlus**, the CFD-DEM codes of CEMF (www.cemf.ir).

Developed by Antonio Di Mauro, MSc in Aerospace Engineering, for the master's thesis carried out at the
Università degli Studi di Napoli Federico II under the supervision of Prof. Matteo Chiatto.

The original codes:

- phasicFlow, the DEM code, by Hamidreza Norouzi: https://github.com/PhasicFlow/phasicFlow
- PhasicFlowPlus, the CFD-DEM coupling with OpenFOAM, by its core developers Hamidreza Norouzi and Bahram Haddadi,
  with contributions from Alireza Hosseini and Nima Joghataei: https://github.com/PhasicFlow/PhasicFlowPlus

Both are free software under GPL-3, and so is this repository. If you use it, please cite the PhasicFlow article:
H. R. Norouzi, PhasicFlow: A parallel, multi-architecture open-source code for DEM simulations,
Computer Physics Communications 291 (2023) 108821, https://doi.org/10.1016/j.cpc.2023.108821

CFD-DEM stack for a T-junction study: the DEM library `phasicFlow` and the
OpenFOAM coupling library `PhasicFlowPlus` in one tree, with the unresolved,
four-way coupled solver `unresolvedSphereSemiImplicitPFPlus` for spheres.

## Origin

`phasicFlow/` and `PhasicFlowPlus/` are the projects of www.cemf.ir (GPL-3),
taken at phasicFlow `454ba4d2` and PhasicFlowPlus `bbf2110`. The first commit
is the two trees as distributed, the second the modifications
(`git diff cf6d3fb main`); every modified file says so in its header.
Added files: Copyright (C) Antonio Di Mauro, GPL-3. An added file derived from a distributed one
also keeps its copyright (www.cemf.ir) and its authors, and names the original file in its header.
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
