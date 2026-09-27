# High JONSWAP HOS failure diagnosis, 2026-09-27

Scope: read frozen inputs, deployed-source control flow and existing logs;
run a bounded remote MATLAB initial-field audit. No HOS rerun, source/binary
change, threshold relaxation, new filter or changed physical input occurred.

## Confirmed termination

Run: `/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-r4gl-80tp-20260926-v2`.
Low kpHs/2=.02 completed all four phases at 1100.8 s, with 5505 records
per phase. High kpHs/2=.12 failed at phi270, MPI rank 5, exit 1. Its last
saved probe time is 5.6 s. The family controller stopped its other phases
and the pipeline terminated before either GL postprocessing stage.

The last solver message is:

```text
Slope along x is too large...  3.001E+00   1.008E+01
STOP 1
```

The exact deployed source, at
`hos-mpi-check-20260926-v1/source/sources/HOS/HOS-ocean.f90`, prints
`time_cur, slopex_max`. Thus 3.001 is internal nondimensional time, not a
slope or a 3 s physical time. **10.08 is the cached maximum x slope**, and
the hard stop is `MAXVAL(ABS(etax)) > 10.0`. For this dimensional finite-depth
initial-condition route, T_out=sqrt(h/g); the internal time corresponds to
approximately 5.7 physical seconds. The guard is an extreme-slope stop,
not an independently certified physical-breaking criterion.

The measured high-family RSS peak was about 7.73 GiB against its 16 GiB
guard. The recorded cause is the slope STOP, not the memory guard, GL
postprocessing or mail. Failure mail was accepted by the local mail system
at 2026-09-27 15:42:52 UTC, without inbox-delivery confirmation.

## Important control-flow limitation

The deployed `HOS-ocean.f90` calls `check_slope` before each adaptive RK
attempt (around line 639). `RK_adapt_2var_3D_in_mo_lin` repeatedly calls
`solveHOS_lin` on trial stages; `solveHOS_lin` writes the shared eta
derivative arrays through `phisxy_etaxy` (resol_HOS.f90 around line 1299).
The RK routine combines its stages into returned spectra without refreshing
those derivative arrays from that final combination.

If a step is rejected, the main loop keeps the old accepted `a_eta/a_phis`
and reduces the step, but the shared derivative arrays still contain a
trial-stage evaluation. With breaking models disabled, the main loop also
does not explicitly refresh them from accepted spectra before the next
slope check. Therefore **the slope guard can inspect a trial-stage cache,
including a rejected trial**, rather than the accepted solution.

This is verified source-control-flow behavior, not proof that the final
failed attempt was rejected: the production log did not preserve that flag,
the corresponding error value or accepted/trial spectra at the stop. It
prevents treating the message as proof of physical wave breaking. RK error
itself is reduced across MPI ranks with MPI_ALLREDUCE/MAX; no rank-local
error-decision discrepancy was found in the inspected routine.

## Initial fields checked remotely

`audit_failed_high_initial.m` read the frozen low/high E/P arrays and all
32 MPI input slabs per family. Eta and true-surface-psi text exports match
the saved arrays with maximum absolute discrepancy **0**. This verifies
file export, not the exactness of the second-order initial model or every
aspect of solver-internal normalization.

The following x slopes were computed by Fourier differentiation and
twofold Fourier interpolation of the initial field:

| Family / phase | Maximum initial absolute eta (m) | Maximum initial absolute x slope | Second-order eta / first-order eta, L2 |
| --- | ---: | ---: | ---: |
| Low, all phases | 1.607--1.895 | 0.0565--0.0650 | 2.67--2.69% |
| High 0 | 11.542 | 0.4768 | 16.00% |
| High 90 | 11.232 | 0.6030 | 16.12% |
| High 180 | 13.536 | 0.5686 | 16.00% |
| High 270 | 11.898 | 0.5074 | 16.12% |

No initial field is already near the stop threshold 10. For high phi270,
modes with radial k>3kp account for about 2.76% of eta variance but 18.91%
of x-slope variance. These are variance-spectrum diagnostics, not physical
total-energy fractions. The retained parents extend to k/kp=4.604 and
f/fp=2.459. Initial short waves are consequently relevant to slope behavior,
but initial spectra alone do not establish a later high-wavenumber cascade.

The high case is kpHs/2=.12 with Hs=8.602 m. It is much stronger than the
earlier phase-only randomization with potential focusing label Akp=.12 and
Hs about .408 m; success of that weak realization does not validate this one.

## What remains unresolved

The run uses M=5, qx=qy=3 partial dealiasing, no dissipation/breaking model,
and Ta=0 with only 11+20+22 initial fields. High-order startup adjustment,
spatial/aliasing effects and strong nonlinear steepening are plausible
contributors; none has been isolated by a controlled replay. Official
HOS guidance also identifies numerical-parameter sensitivity without
dissipation, but recommends q=3 as a common robust starting point. Full
dealiasing must not be advertised as a proven fix:
<https://lheea.gitlab.io/HOS-Ocean/choice-numerical-parameters.html>.

The two-second preflight demonstrated startup and initial-probe agreement,
not stability through the roughly 5.7 s failure or through 80Tp. The
production run retained probes and full initial fields, not evolving full
eta/psi spectra or a restart at the failure; phi270's energy file contains
only buffered header content. The accepted failed-time field cannot be
recovered from those records.

The most discriminating next step is an isolated short replay of phi270,
instrumented to retain the last accepted spectra, trial spectra, previous
RK error/rejection flag and both slope values. Because the derivatives
involve MPI FFTs, any recomputation in diagnostics must be coordinated
across all ranks, not called only on the rank that exceeds the guard.
Keep the threshold and evolution equations unchanged for that replay.
Only after that should separate short comparisons examine dealiasing,
resolution or initialization. Do not restart the whole 80Tp campaign or
raise the slope threshold to hide the failure.

Audit output stays remote under
`/home/lxy/green-laplace-unidirectional-time-series/results/hos_failure_audit/20260927-v1/`.
`initial-fields/report.json` and `initial_diagnostics.csv` contain the exact
numbers; `source/SHA256SUMS` and `source/BASE_COMMIT` identify the snapshot.
The MATLAB audit completed in 56.47 s, at about 1.35 GiB child peak RSS,
with no Code Analyzer messages. Existing run files and unrelated GPU work
were not modified.
