# CDV narrow gap with mass-flow and total-temperature coupling

`CDVGapNoWork.control` uses `MESH/cdv_fan_turbine_1cell_gap_p2.msh`.
The gap extends from x=0.1243120176591453 to 0.1296879823408547 m.

- `fan_inlet`: extraction outflow, initially 6676.328911229583 Pa, 100 Pa above
  the continuous-nozzle reference pressure. The controller adjusts this pressure
  to balance the two gap mass flows, with relative tolerance 0.05%.
- `fan_outlet`: reinjection total inflow, fixed pt=6894.750052334262 Pa and angle
  0.05988380066238474 degrees. Its total temperature is updated every full step
  from the extracted mass-weighted total temperature.
- External inlet: pt=6894.757293168 Pa and Tt=55.55555555555556 K.
- External outlet: fixed static p=6756.86214730464 Pa; other surfaces are slip walls.

Temperature feedback approximates a no-work mean state, with one-step explicit
coupling lag; it does not enforce instantaneous discrete energy-flux equality.
The input uses the NS executable at Reynolds number zero, orders (2,2,1), Roe
flux, RK3 and CFL 0.75. It starts a fresh 100000-step calculation with snapshots
every 10000 steps. Check `.massflow-control` and `.total-temperature-control`
for actual balance and temperature matching.

Run `./configure` from `Solver/` first. From this case directory, with a matching
macOS sequential CPU solver build:

```sh
make -C SETUP ns COMPILER=gfortran PLATFORM=MACOSX COMM=SEQUENTIAL
../../../bin/horses3d.ns CDVGapNoWork.control > RESULTS/CDVGapNoWorkCoupled.log 2>&1
```

Use PLATFORM=LINUX on Linux and compiler/communication options matching the solver.
Outputs use `RESULTS/CDVGapNoWorkCoupled`; older calculation results are preserved.
