# T106A: single-passage inviscid case

`MESH/T106A_1layer.msh` is an ASCII Gmsh 2.2 mesh built from the user-supplied
`T106A_cascade_v01.iges`. It contains **only linear hexahedra** in the fluid
volume and quadrilateral boundary faces, matching HORSES3D's mesh reader.
The mesh is ready to use; Gmsh is only needed to regenerate it.

## Operating conditions and case layout

The operating reference is [HOW5 CS2: DNS/LES of T106 cascades](https://how5.cenaero.be/sites/how5.cenaero.be/files/2025-01/CS2_DNSLEST106_0_0.pdf).
For T106A it specifies **isentropic exit Mach 0.405**, **inlet angle 45.5 degrees**,
**gamma=1.4**, **Pr=0.71**, and an isentropic exit/chord Reynolds number of 60000.
This case uses its Mach number, inlet direction and pressure-ratio specification
as an **inviscid analogue**: `Reynolds number = 0.0` remains the solver's
zero-viscosity switch. It is not a reproduction of the required viscous DNS/LES.

HOW5 allows absolute inlet total pressure and temperature to be chosen.
The existing reference scales **p0=101325 Pa** and **T0=300 K** are retained.
With the solver's air gas constant R=287.15 J/(kg K), the outlet static pressure
and reference state are obtained from:

```text
factor = 1 + (gamma - 1) * M2s**2 / 2
T2s = T0 / factor
p2 = p0 / factor**(gamma / (gamma - 1))
V2s = M2s * sqrt(gamma * R * T2s)
```

This gives **T2s=290.47109570538487 K** and **p2=90500.80687392135 Pa**;
`p0/p2=1.1196032775835698`. The control file uses **Mref=0.405** and this
isentropic exit state as its reference and uniform initial state. Reference
Mach is not a prescribed inlet Mach number, and the actual computed outlet
Mach need not equal the isentropic reference exactly.

### Inlet angle and geometry

The [official HOW5 blade coordinates](https://how5.cenaero.be/sites/how5.cenaero.be/files/2025-01/T106_blade_coordinates.dat)
were compared with the original IGES section used for this mesh. The direction
from the minimum-x leading-edge point to the maximum-x trailing-edge point is
approximately **-30.256 degrees** for the HOW5 coordinates and **-30.262 degrees**
for the IGES sampling. A planar similarity fit also finds only about **0.022
degrees** of rotation. Thus there is no evidence of a 7.8-degree coordinate
rotation requiring conversion of HOW5's 45.5 degrees to 37.7 degrees.

Both the global `AOA theta` and inlet `aoatheta` are **+45.5 degrees**, measured
from the mesh +x direction towards +y; `aoaphi=0`. The velocity direction is
`[cos(45.5 deg), sin(45.5 deg), 0]`. The earlier **37.7-degree** value belonged
to the Cambridge low-speed operating condition; retaining it would mix two
operating conditions. The upstream geometric extension remains usable with
this inlet incidence; its tangent does not prescribe the inlet velocity.

The current Cambridge-derived mesh, physical dimensions and 10 mm single
spanwise layer are retained. HOW5's unit-chord geometry, spanwise periodic
DNS/LES domain and 10%-axial-chord span are not introduced. Front/back remain
slip walls for the requested spanwise-uniform Euler calculation. No cylinders
or inlet turbulence are added. This mesh is not certified as geometrically
identical to the HOW5 mesh merely by the orientation comparison.

## Build and run

The case follows the standard test layout: `T106A.control`,
`MESH/T106A_1layer.msh`, and `SETUP/ProblemFile.f90`.
Run `./configure` from `Solver/` to generate the shared test Makefile,
solver links and `RESULTS/` directory. From this case directory:

```sh
make -C SETUP ns COMPILER=gfortran PLATFORM=MACOSX COMM=SEQUENTIAL
./horses3d.ns T106A.control
```

For Linux use `PLATFORM=LINUX`; compiler and communication settings must match
those used for the solver. Clean the case objects before switching builds.
The mesh is supplied ready to run; no mesh-generation script or CAD dependency
is needed. Output is written to `RESULTS/T106A.hsol`; restart is disabled.

The numerical configuration is degree-2 split-form DGSEM with Chandrasekar
volume flux, Gauss-Lobatto nodes, Rusanov interface flux and SSPRK33.
CFL=DCFL=0.25, 10000 steps, output interval 100 and tolerance 1e-13.
`SETUP/ProblemFile.f90` initializes a uniform reference state with zero source.

## Geometry and scope

The IGES has surface geometry rather than a closed fluid solid. A planar
blade contour at z=0 was extracted from its CAD curves and used to construct
the passage between that blade's suction side and the adjacent blade's
pressure side. The adjacent blade is translated by **158 mm** in y. Upstream
and downstream extensions terminate at the IGES inlet and outlet positions:
**x=-105 mm** and **x=275 mm**.

As requested, this is a **single passage, single spanwise layer**, without
the upstream small cylinders. The original 375 mm span is replaced by a
**10 mm** layer for spanwise-uniform inviscid calculation. All mesh coordinates
are exported in **metres**; reference length 1 m preserves these dimensions.

The mesh is a **direct structured H-grid**, with three logical blocks:

```text
inlet block          blade passage block          outlet block
 24 x 32 x 1             128 x 32 x 1               32 x 32 x 1
```

The blocks share complete conforming node rows and form a single regular
184 x 32 x 1 logical grid. Every interior point of the planar grid has exactly
four incident quadrilaterals. Native eight-node hexahedra are created directly
by extruding these i-j cells across one spanwise layer. There is no triangular
or polyhedral precursor and no centroid/midpoint subdivision.

CAD curves are sampled at approximately 0.15 mm, then each blade branch is
resampled into 128 wall edges with mild leading/trailing-edge clustering.
The rounded caps are split using support directions along the upstream and
downstream passage extensions (37.7 and -63.2 degrees). Those construction
angles locate geometric block junctions; the inlet velocity remains 45.5 degrees.
Both blade branches are retained, including the rounded leading/trailing edges.
The linear wall edges approximate the original curved CAD contour.

A discrete elliptic equation supplies the interior seed. Fixed-boundary
optimization then improves corner angles and limits excessive size variation
and in-plane edge ratios. The default is 4000 geometry optimization iterations;
this is a fixed smoothing budget, not a flow convergence criterion.
The delivered grid passes the generator's minimum scaled-Jacobian threshold
of 0.35. This mesh has no viscous boundary-layer refinement.

| Quantity | Previous mesh | Structured mesh |
| --- | ---: | ---: |
| Nodes | 33534 | 12210 |
| Hexahedra | 16509 | 5888 |
| Boundary quadrilaterals | 33532 | 12208 |
| Minimum corner scaled Jacobian | 0.281137 | 0.397603 |

The hexahedron count is reduced by approximately 64.3%.

## Boundary names

| Physical surface | Faces | Intended condition |
| --- | ---: | --- |
| inlet | 32 | TotalInflow, with the desired total pressure/temperature and inlet direction |
| outlet | 32 | Outflow, with the desired static pressure |
| blade | 256 | FreeSlipWall |
| periodic_lower | 56 | Periodic, coupled to periodic_upper |
| periodic_upper | 56 | Periodic, coupled to periodic_lower |
| front | 5888 | FreeSlipWall |
| back | 5888 | FreeSlipWall |

The periodic surfaces each comprise two disconnected extensions, before and
after the blade. They match through the translation **(0, 0.158, 0) m**, including
matching face connectivity. Periodic nodes use the master's shared coordinates
exactly to satisfy the solver's strict comparison tolerance.

Use the `.ns` solver with `Flow equations = NS`, `Reynolds number = 0.0`,
`LES model = none`, and disabled gravity for zero physical viscosity and heat
conduction. Retain `compute gradients = .true.`, `gradient variables = state`,
and `viscous discretization = BR1` for the existing GPU implementation.
Use `mesh file name = MESH/T106A_1layer.msh` when running from this directory.
The supplied operating values and their assumptions are described above;
an IGES file alone does not specify them. The periodic definitions are:

```text
#define boundary periodic_lower
   type = Periodic
   coupled boundary = periodic_upper
#end
#define boundary periodic_upper
   type = Periodic
   coupled boundary = periodic_lower
#end
```

## Mesh provenance and validation scope

The supplied mesh was constructed from `T106A_cascade_v01.iges`, using CAD
sampling and a structured H-grid. The original CAD SHA-256 was recorded as
`221d6dc9a4664ae236d6c50f52d95321cae6c13e4d88a8b8645f1d02c9a08ba0`.
The case retains the ready-to-run mesh; geometry-generation tools and local
backup paths are not part of the test distribution.

This cleanup checks case layout, mesh boundary names, element connectivity
and configuration references. It does not establish steady convergence,
GPU performance, spatial convergence or agreement with the HOW5 viscous case.
