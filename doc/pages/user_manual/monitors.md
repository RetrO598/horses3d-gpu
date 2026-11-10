title: Monitors
---

[TOC]

The monitors are specified individually as blocks in the control file.
The only general keyword that can be specified is explained in the table below:

| Keyword                 | Description                                                                                   | Default value |
|-------------------------|-----------------------------------------------------------------------------------------------|---------------|
| monitors flush interval | *INTEGER*: Iteration interval to flush the monitor information to the monitor files.          | 100           |


## Residual Monitors

## Statistics Monitor 

```markdown
#define statistics
   initial time      = 1.d0
   initial iteration = 10
   sampling interval = 10
   dump interval     = 20
   @start
#end
```

By default, the statistic monitor will average following variables:

<div class="multicol">
    <div style="float: left; width: 33.33%;">
        <ul>
            <li>u</li>
            <li>v</li>
            <li>w</li>
        </ul>
    </div>
    <div style="float: left; width: 33.33%;">
        <ul>
            <li>uu</li>
            <li>vv</li>
            <li>ww</li>
        </ul>
    </div>
    <div style="float: left; width: 33.33%;">
        <ul>
            <li>uv</li>
            <li>uw</li>
            <li>vw</li>
        </ul>
    </div>
    <div style="clear: both;"></div>
</div>



A keyword preceded by @ is used in real-time to signalize the solver what it must do with the statistics computation:

<div class="multicol">
    <ul>
        <li>@start</li>
        <li>@pause</li>
        <li>@stop</li>
        <li>@reset</li>
        <li>@dump</li>
    </ul>
</div>

After reading the keyword, the solver performs the desired action and marks it with a star, e.g. @start*.

@note 
Real-time keywords may not work in parallel MPI computations. It depends on how the system is configured.
@endnote

## Probes

```markdown
#define probe 1
   name     = SomeName
   variable = SomeVariable
   position = [0.d0, 0.d0, 0.d0]
#end
```

| Keyword   | Description                                                 | Default value       |
|-----------|-------------------------------------------------------------|---------------------|
| name      | *CHARACTER*: Name of the monitor.                           | **Mandatory Keyword** |
| variable  | *CHARACTER*: Variable to be monitored. Implemented options are: pressure, velocity, u, v, w, mach, k. | **Mandatory Keyword** |
| position  | *REAL(3)*: Coordinates of the point to be monitored.        | **Mandatory Keyword** |


## Surface Monitors

```markdown
#define surface monitor 1
   name              = SomeName
   marker            = NameOfBoundary
   variable          = SomeVariable
   reference surface = 1.d0
   direction         = [1.d0, 0.d0, 0.d0]
#end
```
| Keyword           | Description                                                                                      | Default value       |
|-------------------|--------------------------------------------------------------------------------------------------|---------------------|
| name              | *CHARACTER*: Name of the monitor.                                                                | **Mandatory Keyword** |
| marker            | *CHARACTER*: Name of the boundary where a variable will be monitored.                             | **Mandatory Keyword** |
| variable          | *CHARACTER*: Variable to be monitored. Implemented options are: mass-flow, flow, pressure-force, viscous-force, force, lift, drag, pressure-average. | **Mandatory Keyword** |
| reference surface | *REAL*: Reference surface [area] for the monitor. Needed for "lift" and "drag" computations.    | --                  |
| direction         | *REAL(3)*: Direction in which the force is going to be measured. Needed for "pressure-force", "viscous-force" and "force". Can be specified for "lift" (default [0.d0,1.d0,0.d0]) and "drag" (default [1.d0,0.d0,0.d0]). | --                  |

## Volume monitors

Volume monitors compute the average of a quantity in the whole domain. They can be scalars(s) or vectors(v).

```markdown
#define volume monitor 1
   name     = SomeName
   variable = SomeVariable
#end
```

| Keyword  | Description                                                                                                    | Default value       |
|----------|----------------------------------------------------------------------------------------------------------------|---------------------|
| name     | *CHARACTER*: Name of the monitor.                                                                              | **Mandatory Keyword** |
| variable | *CHARACTER*: Variable to be monitored. The variable can be scalar (s) or vectorial (v). Implemented options are: kinetic energy (s), kinetic energy rate (s), enstrophy (s), entropy (s), entropy rate (s), mean velocity (s), velocity (v), momentum (v), source (v). | **Mandatory Keyword** |


## Load Balancing Monitors

Load balancing monitors compute the DOF in each partition of the mesh to check the unbalance. 

```markdown
#define load balancing monitor 1
   name     = SomeName
   variable = SomeVariable
#end
```

| Keyword  | Description                                                                                                     | Default value       |
|----------|-----------------------------------------------------------------------------------------------------------------|---------------------|
| name     | *CHARACTER*: Name of the monitor.                                                                               | **Mandatory Keyword** |
| variable | *CHARACTER*: Variable to be monitored. Implemented options are: max dof per partition, min dof per partition, avg dof per partition, absolute dof unbalancing, relative dof unbalancing. | **Mandatory Keyword** |


The variable *absolute dof unbalancing* is computed as follows:

\[
\text{absolute dof unbalancing} = \text{max dof per partition} - \text{min dof per partition}
\]

The variable *relative dof unbalancing* is computed as follows:

\[
\text{relative dof unbalancing} = \frac{100 \times \text{absolute dof unbalancing}}{\text{avg dof per partition}}
\]




## Boundary mass-flow transient detection

Enable `mass flow statistics = .true.` in an inlet or outlet boundary that has
an active `mass-flow` surface monitor. The detector stores each actual complete-step
flow observation, in order, without a fixed time grid, temporal interpolation or
subsampling. CFL-dependent observation times are retained as metadata; the ADF
regression uses the discrete observation sequence.

Defaults are `mass flow statistics adf lags = 200`,
`mass flow statistics minimum samples = 1000` and
`mass flow statistics hop samples = 1000`. The minimum cannot be less than 1000
and must provide residual degrees of freedom for the selected lag. No check is
performed before that minimum.

History uses contiguous arrays that grow automatically when their initial capacity
is exceeded; old samples are retained until a reset.

Checkpoints store the current detector history and averaging state directly,
with configuration, mesh, iteration and time consistency checks during restore.
A pressure update still resets both observed boundaries;
both must pass before their control averaging begins. Time-weighted averaging is
unchanged; the removal of interpolation concerns the ADF input history.

For controlled boundaries, the manager tracks waiting for detection, detection
passed, and block averaging as explicit states. A completed averaging block
contains its mean, standard deviation, time interval and sample count; the latest
two blocks are exposed to the controller. Pressure-stage numbering belongs to the
controller, which requires successful checks after the current stage began.

Detection passes when the p-value is below `mass flow statistics significance`
and the relative standard deviation is below
`mass flow statistics relative std tolerance`. These per-boundary settings apply
to both standalone monitoring and controller-managed monitoring. The controller
uses the detector result without imposing a second set of detection thresholds.

Observed zones are bound to enabled monitor entries once during controller
construction. Subsequent reads, averaging and resets use those fixed entries.
The manager advances the group from the same complete-step observations and
starts all averaging blocks together, so the controller does not revalidate
detector epochs, observation steps or block alignment at runtime.

Mass-flow balance controls the outlet pressure using mass-flow statistics only;
CFD residuals are not a controller convergence criterion. After every pressure
change, all observed boundaries must pass detection after the stage began and
`minimum settle steps` must elapse before block averaging starts. During block
averaging, the detector is not rerun: block-to-block drift and standard deviation
gate response quality. Poor-quality blocks are skipped while averaging continues.

A flow response is accepted when its absolute error is at most
`absolute mass tolerance (kg/s)` plus `relative mass tolerance` times the reference
flow scale (bounded below by `flow floor (kg/s)`). The controller then collects two
fresh blocks at the same pressure. Only a second acceptable response within the
flow tolerance marks the balance complete. A failed confirmation is reused for
the next pressure proposal.

There is no controller timeout, outer evaluation limit, or statistics retry limit.
The obsolete `maximum settle steps`, `maximum outer evaluations`,
`maximum statistics retries`, and `final cfd residual tolerance` options are no
longer used. Pressure limits and `maximum pressure step (pa)` still apply. When
pressure exploration is restricted or there is no throughflow, the current
pressure is held and observation continues. A restriction is logged once until
its reason changes or pressure is updated. Nonfinite responses and invalid input
remain fatal errors.

Once balanced, the controller holds the last pressure while the solver continues
according to its original convergence, time and iteration limits. The controller
does not extend the flow solve. If the solver finishes first, the control log
reports `calculation-ended-unbalanced` and the final pressure. Accepted flow errors
are recorded in the response rows. Pressure-response history grows automatically; repeated
observations at the current history point replace its response.

The `.massflow-control` log records one `control-start` event per run, one
`mean-start` event for each group averaging start (including confirmation),
accepted constraint responses, pressure updates/restrictions and the final result.
Control events include the iteration and pressure-stage number. Per-boundary ADF
diagnostics and averaging data are written only to the boundary `.massflow` logs;
they are not copied into the control log. The ADF `check` rows also supply detection
pass/fail outcomes without separate detection-transition messages. Standalone
monitors retain their own `average-start` event.

Boundary statistics and pressure control share one extension inside each `.hsol`
snapshot. The time integrator writes the indexed boundary records, then an optional
`MASSFLOW-CONTROL` record, and seals a single footer locating the extension. No
separate controller checkpoint is created, and the field-data offsets are unchanged.

The control record stores the current pressure, flow scale, pressure-stage number,
and only the valid pressure/flow-error response points. A compact configuration
summary contains the constraint, target, boundary names/order and flow coefficients.
Those must match on resume. Tolerances, relaxation, averaging duration and pressure
step limits may change; current and historical pressures must fit the new bounds.
Record offsets, lengths, counts and finite numerical values are checked without
fingerprinting the solution or duplicating its iteration/time metadata.

`restart policy = resume` restores completed responses and the current pressure;
stability and partial averages are collected again, including at a previously
balanced pressure. `restart policy = reset` discards controller history and uses
the pressure configured in the control file. The current record format is required
for resume; previous controller checkpoint formats are not supported.
Pressure changes do not overwrite solver residuals. Enabling balance does not
disable p-adaptation, FMG initialization or IBM.
