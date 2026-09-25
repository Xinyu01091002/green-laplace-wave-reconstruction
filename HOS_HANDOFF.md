# HOS-Ocean handoff

Current local project:
`C:/Users/spet5947/Documents/green-laplace-unidirectional-time-series`

Branch: `codex/unidirectional-time-series`.
Repository: <https://github.com/Xinyu01091002/green-laplace-wave-reconstruction>.
Remote alias: `60.188.112.99:60093`, user `lxy`.
Remote project: `/home/lxy/green-laplace-unidirectional-time-series`.

## OW3D background campaign

The user explicitly authorized sixteen independent OW3D jobs and completion
email, while HOS work continues in this project. Sixteen means **one focused
wavegroup plus three random realizations, each at four global phases**, not
sixteen duplicated jobs or sixteen four-phase families.

Run ID: `ow3d16-kpd1-akp012-20260925T212822Z`.
Run root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/ow3d16-kpd1-akp012-20260925T212822Z`.
Controller screen: `gl-ow3d16-20260925T212822Z`.
Email screen: `gl-ow3d16-mail-20260925T212822Z`.

The frozen `manifest.json` contains the execution source commit, source and
input archive SHA-256 values, copied solver hash, expected per-case times,
individual input hashes and hashes of sixteen distinct physical payloads.
Its `source/`, `bin/`, `inputs/` and `cases/` are separate from the active
project. Do not edit or overwrite a running case, or retarget it to changing
HOS source. No automatic local fetch or remote deletion is configured.

| Family | Seed | Global phases |
| --- | --- | --- |
| Focused wavegroup | None | 0, 90, 180, 270 degrees |
| Random realization 1 | 20260925 | 0, 90, 180, 270 degrees |
| Random realization 2 | 20260926 | 0, 90, 180, 270 degrees |
| Random realization 3 | 20260927 | 0, 90, 180, 270 degrees |

All use kpd=1, kp=.0279 1/m, h=35.8422939068 m, g=9.81 m/s2,
50 lambda by 20 lambda, 1025x257x17 native nodes, dt=.2 s and end time
220 s. There is no extra depth/steepness/dt sweep. The physical field grid
has 1024x256 unique points. Directional Gaussian **amplitude** width is
25 degrees; this is not a 25-degree energy RMS width.

The randomized fields preserve the focused family's first-order modal
amplitudes exactly apart from roundoff: no taper or Hs normalization.
Akp=.12 labels potential focusing steepness. Linear spatial sigma is
0.1019037931 m, Hs=4*sigma=0.4076151725 m and kp*Hs/2=.005686231656.
Do not label them as random sea-state steepness .12.

All sixteen OW3D runs use **native walls**, including the untapered random
fields. They are finite-time wall-bounded experiments, not automatically
the same boundary-value problem as periodic HOS-Ocean. Check reflected
components and boundary influence before using these fields as GL/HOS
accuracy references. Initial eta/psi pairs are independent MF12 order 2;
there are no GL high-order reference inputs, MF12 third-order calculations,
or historical VWA/Stokes-style higher-order additions.

Native kinematics saves all x and y indices 125--133, with time indices
2--1101. At t=0, read eta/psi from EP_00000.bin: this binary writes an
unsolved volume phi in its t=0 kinematics record. Sparse EP is saved every
55 s. Data remain dimensional and unaligned; no gain, bias, phase or shift
is fitted. First-order phase-sector extraction is not exact perturbation
order separation after nonlinear propagation.

## Status and email

Query the real controller state, not the existence of a screen alone:

```sh
run=/home/lxy/green-laplace-unidirectional-time-series-runs/ow3d16-kpd1-akp012-20260925T212822Z
cat "$run/status.txt"
cat "$run/status.json"
tail -30 "$run/runtime.log"
cat "$run/mail-watcher-status.txt"
pgrep -x ow3d
```

`status.json` contains per-case PID, state, observed RSS high-water mark and
kinematics file bytes. `started.utc` and `finished.utc` are atomic markers.
No forced simulation timeout replaces the requested 220 s evolution.
The earlier 300 s two-step diagnostic timed out before an advance; it is
not a measured per-step duration or a successful propagation test.

On native process exit, the controller requires exit 0, `JOB IS COMPLETE`
and final time 220 s within 1e-8 s. MATLAB then validates all eta/phi output
records and saves `cases/<id>/processed/surface_strip.mat` with t, x, y,
eta, psi and an IO report. MATLAB extraction is serial, so it does not
start sixteen postprocessors on top of active OW3D jobs. Surface-potential
differences against sparse EP are reported without correction; warnings
do not become an accuracy claim.

After all jobs and validators close, the root `summary.json`, `exit-code.txt`
and `finished.utc` are written. The email watcher then sends one completion
or failure notice using the existing `/usr/local/bin/run.sh` recipient
configuration. No recipient or authentication information is copied into
this project. `mail-status.txt=accepted_by_local_mail` means local mail
acceptance, not confirmed inbox delivery. No test email is sent.

## HOS work to do in this project

1. Read `research/directional_wave_data/HOS_OCEAN_ASSESSMENT.md`. HOS has not
   yet been installed, run or benchmarked in this task. Use a fixed official
   release or commit, record its binary hash, and prefer the documented
   Linux prebuilt binary if suitable. No Fortran compiler was found on PATH.
2. Use the identical parent spectra and independent eta/true-surface-psi
   fields from the frozen run's `inputs/<family>/initial_fields.mat`.
   They contain C, kx, ky, om, E/P for four phases and exact conventions.
   Do not reuse the old tapered random prototype. For random fields, HOS
   periodic evolution is the intended comparison; it differs from OW3D walls.
3. Target the same 1024x256 unique horizontal physical grid and physical
   domain; verify the selected HOS version's real/modal indexing. Do not
   import OW3D duplicate endpoints, ghost nodes or Nz. Check dimensional
   `Initial surface quantities` input, true surface-potential meaning,
   Fourier normalization, signs, spatial origin and the t=0 round trip.
4. Initial assessment proposes HOS M=5 and qx=qy=3, with output every .2 s.
   M is an evolution parameter, not an eta33 label; q=3 is partial dealiasing
   at M=5. Record adaptive integration tolerance and nonlinear activation
   explicitly. Do not silently apply a long initial ramp, breaking model
   or dissipation that changes the supplied state. The cancelled OW3D dt
   sweep must not be reinstated.
5. Check eta/phi output and one short HOS case before making performance or
   accuracy claims. Prefer surface physical/modal output and remote MATLAB
   extraction; verify probe columns before assuming they contain phi.
   Keep GL no-Stokes kernels and all required surface Taylor terms unchanged.
6. While OW3D runs, allow for its observed startup memory: 16 times 23.38 GiB
   plus 25% is about 468 GiB. This is a budget, not a complete-run peak or
   throughput measurement. Recheck live RAM/CPU before starting HOS jobs.

The active local branch may evolve for HOS work; the OW3D campaign executes
its frozen source and copied native binary independently. Preserve the run,
all failed pilot evidence and historical SWORD directories. The completion
mail is attached to this OW3D run, not to future HOS work.
