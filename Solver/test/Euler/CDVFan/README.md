# CDV nozzle with a prescribed fan

`CDVFan.control` uses `MESH/cdv_fan.msh`, a split nozzle with a gap from
x=0.1143 to 0.1397 m. The two fluid regions are disconnected.

- `fan_inlet`: extraction outflow; its static pressure is adjusted to balance
  extraction and reinjection mass flows.
- `fan_outlet`: total-condition inflow with prescribed ideal fan pressure ratio
  1.05 and the corresponding isentropic total-temperature rise.
- External inlet and outlet conditions are fixed; the other surfaces are slip walls.

This case uses the NS executable with Reynolds number zero, solution orders
(2,2,1), Roe flux and RK3. It starts a fresh 100000-step calculation. Mass-flow
balance must be checked in `.massflow-control`; finishing the iteration budget
alone does not establish convergence.

Run `./configure` from `Solver/` first. From this case directory, with a matching
macOS sequential CPU solver build:

```sh
make -C SETUP ns COMPILER=gfortran PLATFORM=MACOSX COMM=SEQUENTIAL
../../../bin/horses3d.ns CDVFan.control > RESULTS/CDVFan.log 2>&1
```

Use PLATFORM=LINUX on Linux and compiler/communication options matching the solver.
Outputs are stored under `RESULTS/CDVFan`.
