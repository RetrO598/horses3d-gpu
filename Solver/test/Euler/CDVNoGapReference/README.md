# Continuous CDV nozzle reference

`CDVNoGapReference.control` uses `MESH/cdv_planar_half_p2.msh`, the continuous
500-element nozzle mesh, as a reference for the narrow-gap calculation.

- Inlet: fixed pt=6894.757293168 Pa, Tt=55.55555555555556 K and direction +x.
- Outlet: fixed static p=6756.86214730464 Pa.
- Other surfaces: slip walls. Mass-flow control is disabled.

It uses the NS executable at Reynolds number zero, solution orders (2,2,1),
Roe flux, RK3 and CFL 0.75. It starts a fresh 100000-step calculation with
snapshots every 10000 steps.

Run `./configure` from `Solver/` first. From this case directory, with a matching
macOS sequential CPU solver build:

```sh
make -C SETUP ns COMPILER=gfortran PLATFORM=MACOSX COMM=SEQUENTIAL
../../../bin/horses3d.ns CDVNoGapReference.control > RESULTS/CDVNoGapReference.log 2>&1
```

Use PLATFORM=LINUX on Linux and compiler/communication options matching the solver.
Outputs use `RESULTS/CDVNoGapReference`. Compare against the gap case at common
physical sections outside the gap, using converged states or consistent physical
times; identical iteration counts need not represent identical physical times.
