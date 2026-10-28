# Inviscid bump with total-condition inflow

This case uses the supplied FEST-3D smooth-bump mesh, converted to ASCII
Gmsh 2.2 with one spanwise layer. The mesh contains 4608 linear hexahedra.
The inlet specifies total pressure and temperature; the outlet specifies
static pressure. All four remaining surfaces are slip walls.

The case runs through the compressible NS executable with
`Reynolds number = 0.0` for zero physical viscosity. Its reference Mach
number is approximately 0.5. The control file specifies degree 2,
Gauss-Lobatto nodes, Roe interface flux and RK3.

## Build and run

The layout follows the other test cases:

```text
InviscidBump.control
MESH/smooth_bump_97x49x2.msh
SETUP/ProblemFile.f90
```

Run `./configure` from `Solver/` to generate `SETUP/Makefile`, solver links
and `RESULTS/`. From this case directory, for a matching macOS CPU build:

```sh
make -C SETUP ns COMPILER=gfortran PLATFORM=MACOSX COMM=SEQUENTIAL
./horses3d.ns InviscidBump.control
```

Use `PLATFORM=LINUX` on Linux. Match the solver's compiler and communication
settings, and clean case objects before switching builds. No conversion script
is needed to run the supplied mesh. Output is `RESULTS/InviscidBump.hsol`;
restart is disabled by default. Generated Makefiles, libraries, mesh caches
and results follow the repository-wide ignore rules.
