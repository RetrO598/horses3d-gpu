# VKI turbine cascade with total-condition inflow

This case prescribes inlet total pressure 90400 Pa, total temperature 404.1 K
and zero inlet angle, with outlet static pressure 52300 Pa. The blade is an
adiabatic no-slip wall, `periodic_1` and `periodic_2` form the periodic pair,
and `front` and `back` are slip walls.

The case and output names omit the turbulence-model suffix. The existing
physical configuration remains `Flow equations = NSSA`, Reynolds number
6.0e6; it requires the Spalart-Allmaras executable and matching problem library.

## Mesh and layout

```text
VKITurbine.control
MESH/vki_inviscid_p3_mesh.h5
SETUP/ProblemFile.f90
```

The active input is a HOPR-format HDF5 mesh with geometry order 3.
The filename's `inviscid` describes the supplied mesh name, not the equations
selected by this control file. The solution polynomial order remains 2.
The solver must be built with `WITH_HDF5=YES` and the corresponding HDF5
dependency settings. All seven mesh boundary names are defined in the control
file. The previous unused Gmsh mesh and its generated geometry/Tecplot files
were moved to the ignored local `RESULTS/legacy_mesh/` directory during cleanup;
they are not required to run this case.

## Build and run

Run `./configure` from `Solver/` to generate the shared test Makefile,
solver links and `RESULTS/`. From this case directory, for a matching macOS
CPU build:

```sh
make -C SETUP libproblemfile_nssa.dylib COMPILER=gfortran PLATFORM=MACOSX COMM=SEQUENTIAL
./horses3d.nssa VKITurbine.control
```

On Linux use the target `libproblemfile_nssa.so` and `PLATFORM=LINUX`.
Compiler and communication settings must match the solver; clean case objects
before switching builds. Output and restart paths are `RESULTS/VKITurbine.hsol`;
restart is disabled by default. Existing local results retain their historical
filenames. This cleanup does not validate convergence or mesh suitability for
wall-resolved turbulence calculations.
