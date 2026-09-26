# HOS-Ocean handoff

## Completed checkpoint: 2026-09-26

Both directional focused-wavegroup HOS families and full-record GL processing
are complete. This section supersedes the historical launch status below.
See `research/directional_wave_data/hos_low_run/COMPLETED.md` and the saved
compact metrics alongside it. Raw records and plots remain in ignored artifacts
or the original remote run directories.


## ACTIVE: Akp=.02 directional run with automatic GL, 2026-09-26

Read `research/directional_wave_data/hos_low_run/README.md` first.
Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1`.
Screen: hos-akp002-auto-gl-20260926. HOS progress: status.json. Whole workflow:
pipeline-status.json. Started about 13:10:33 UTC. All four initial probe
checks passed and live nonzero-time advancement was verified.

User explicitly required step 3 to be short: MPI4 and MPI8 each ran ONLY
2 physical seconds, taking 51.48 / 29.57 wall seconds. Output/probe agreement
passed (max eta/psi relative L2 9.19e-11 / 7.11e-11). Production uses only
chosen MPI8: four concurrent phases, 32 HOS compute ranks, CPU IDs 8--39.
Settings otherwise match the completed .12 case: kph1, M5/q3, abs tol1e-12,
Ta0, 1024x256 over50x20lambda, duration220s, output.2s. Fresh MF12 order2
initial fields were recomputed at Akp=.02, not a linear rescale of total eta.

Only 8 OW3D jobs remain: focused wavegroup plus random seed20260925, each
four phases. Seed20260926/20260927 jobs were intentionally terminated on user
request, all files retained. Original frozen controller reports exit-15 as
failed; use its user-requested-reduction-20260926.json for cancellation intent.
Do not restart canceled jobs. No GPU computation was involved.

After low HOS completion and native-record checks, pipeline.py automatically
runs complete GL time-series processing for low and then high cases. The
reserved outputs are full-gl-comparison-20260926-v1 in each run; avoid creating
those paths independently while the pipeline is active. Both use identical
input-defined frequency rules and symmetric prediction/reference filtering.
Low GL results are not yet available at this launch checkpoint.

Akp=.12 HOS is fully complete: 1101 records through220s in all four phases.
Concurrent batch wait was 2h33m19s. Per-phase wall: 2:23:53 / 2:25:01 /
2:33:16 / 2:32:31 for 0/90/180/270 deg. Historical stale progress fields
showing200s do not override the verified final native records.


## Partial directional GL comparison completed; production continues

A frozen common 0--173.6 s prefix (869 samples) was taken at 2026-09-26
06:13:49 UTC from v2. Read
`research/directional_wave_data/hos_production/PARTIAL_COMPARISON.md`.
GL J8, 7.5-degree directional allocation plus observed first harmonic gives
main-group eta22 L2 percentages: center 4.22797%; y=-/+52.782 m 3.10274/3.14256%;
y=-/+70.376 m 2.58541/2.62474%. All points pass the one-third amplitude gate.
Observed and predicted second harmonics use an identical output projection;
raw comparisons are also retained. Source-derived first-harmonic frequency
support is fixed before scoring, not target-tuned.

A 10 s shorter center record changes the GL main-group prediction by .08291%
and the error metric from 4.22797% to 4.21392%. This is a useful preliminary
result, not the final 220 s result. No HOS run was paused or changed. Extra
postprocessing used one CPU and about 1.51 GiB peak RSS for 113.40 s wall.
Remote preview directory under the active run: partial-comparison-20260926-v1.
Local plot/metrics: artifacts/hos_ocean/directional-partial/.


## ACTIVE directional production: 2026-09-26 04:12:14 UTC

Four phases, four MPI ranks each, are now running. Read
`research/directional_wave_data/hos_production/LAUNCH.md` for full settings.
Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2`.
Screen: hos-directional-4x4-v2-20260926. Status is in status.json; do not infer
completion from screen existence. All four initial-probe checks passed and
actual nonzero-time advancement was observed. Full 220 s evolution remains
in progress; no GL comparison result is claimed yet.

Settings: kph=1, Akp=.12, 25-degree amplitude angular sigma, frozen MF12
order-2 eta/psi (11+20+22), 1024x256 points over 50x20 lambda_p, M5/q3,
Cash-Karp adaptive RK5(4), absolute tolerance 1e-12, Ta=0. Run 220 s =
15.986Tp; nominal focus 110 s. Five native-grid eta probes, output .2 s.
No dense spatial or volumic output; initial full eta/psi and spectrum retained.
Each phase is pinned to a distinct four-CPU set within CPU IDs 16--31.
Combined sampled HOS RSS about 3.9 GiB, 8 GiB guard; existing 16 OW3D jobs
are unchanged. Per-rank time/RSS and automatic output validation are enabled.

IMPORTANT: v1 production was stopped by its initial-probe guard. Upstream
MPI probe IO omitted broadcast of nprobes. v2 broadcasts it and flushes
probe records; no evolution equation changed. Four-rank phase-90 native-field
comparison passed at t=0,.2 s to max 1.59e-15 m. Use the v2 copied binary:
SHA256 93d64252f732449fad9d2171979c91e13b0b9de1d17a73d7587fa00721efcbb6.
The old MPI binary remains valid for the recorded spatial-field comparison,
but must not be used for multi-rank probe output. Failed v1 evidence is kept.


## Latest: directional MPI 1-versus-4 check passed, 2026-09-26

Read `research/unidirectional_time_series/hos_ocean/mpi/README.md`.
Same 1024x256 directional phase, M5/q3, 2 s physical evolution. MPI-1:
135.67 s; MPI-4: 67.25 s (2.017x speedup). Sampled aggregate job RSS including
launcher/descendants: .901 / .973 GiB. Sum rank CPU time: 135.04 / 267.25 s.
Four ranks shorten single-task wait but use about twice the total CPU time.
All 11 eta/psi frames agree between MPI-1 and MPI-4 to relative L2 <4e-14;
both agree with prior serial reference to <1e-10. No fitted correction.

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-mpi-check-20260926-v1`.
Binary: build/sources/HOS-Ocean; launch using its private environment.json.
SHA-256: `b00f8d20a4d6ba82c472bb239e8577f1e6c1b5a1c6320edb3257831f75c2edde`.
MPI/FFTW-MPI were unpacked privately, not system-installed. The existing
IO-only precision patch is retained. Initial files are y-slab partitions;
outputs must be assembled using rank-local coordinates.

Two short MPI jobs completed sequentially, maximum four compute ranks.
No 8-rank run or multi-phase MPI production campaign was launched. Existing
16 OW3D processes were not altered; available RAM stayed about 472 GiB.
Local summary: artifacts/hos_ocean/mpi-results/summary.json.

## Latest result: GL time-series check against HOS and existing OW3D, 2026-09-26

Completed the user-authorized near-matched single-point eta22 test. Read
`research/unidirectional_time_series/hos_ocean/timeseries/RESULTS.md` first.
Settings: kph=1, Alpha=1, Akp=.02/.12, 4096 periodic points, 68 lambda_p,
focus at .67L and elapsed 40Tp; no31 initialization remains in HOS. Four
phases per family, 50Tp duration, sample Tp/40, 641 samples over 32--48Tp.
Original GL6/8/12/MF12/VWA/Walker processing is frozen and reused; temporal
Hilbert-four-phase signs were retained. The old OW3D prediction matrices
were reproduced to about 2e-15 relative norm, confirming processing parity.

Main-group GL6 errors: HOS .26350% / 3.97363%, OW3D .25921% / 3.93437%.
GL12: HOS .10855% / 3.88777%, OW3D .10379% / 3.84863% (low / high).
The same output support filter on prediction/reference changes these
numbers negligibly; no fit or target-selected cutoff was used.
This supports HOS as a practical reference for this eta22 time-series test;
it does not resolve prior eta33 questions or certify every wave condition.

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1`.
Local plots/metrics: `artifacts/hos_ocean/timeseries-results/`.
HOS time per phase: low 292--299 s, high 383--389 s; peak about 18.4 MiB.
This run's 9 HOS processes (8 production + 1 IO smoke) sum to 2726.69 s
process wall time, about 688.76 s active batch waiting with concurrent phases.
All HOS runs so far: 42 processes, cumulative process wall time 4895.00 s.
MATLAB times and RSS are recorded separately. No new OW3D job was launched,
and the frozen 16-case directional OW3D campaign was not modified.

A new IO-only build extends probe text to ES25.16E3. Binary SHA-256:
`13c4bae40855f3d3a29da12d5567559c4a41e415c27bc5f67d5b927db6c2f060`.
The old probe-header parsing failure and subsequent resume are preserved.
Raw data remain remote. No commit/push was performed for this work.

## Four-phase audit update, 2026-09-26

Read research/unidirectional_time_series/hos_ocean/deepwater/SPECTRUM_AUDIT.md.
Original spatial Hilbert-four-phase matrix agrees with the previously used
positive-k operator (eta relative difference <=7.1e-12). No material sign or
first/third swap was found. Earlier generic DFT/extra-phase explanations were
not a demonstrated diagnosis. The original plot also applies harmonic
bandpass filtering; the previous raw HOS comparison did not reproduce that
full workflow. Unfiltered separated spatial spectra are now plotted and saved.
No new HOS simulations were run in this audit. All 33 HOS executions so far
sum to 2168.31 s process wall time; the 28 deep-water runs sum to 1629.02 s.


## Latest experiment: 2026-09-26, unidirectional deep-water four phases

User requested spectral MF12 initialization without 31 and reconstruction
of second/third harmonics from current first harmonic at 3Tp and 20Tp,
low then high steepness, with runtime and memory measurement. This new scope
supersedes earlier order-two-only/no-MF12-third-order limits for this experiment.

Completed Akp=.02/.12, Alpha=1, kp=.0279, kph=20, 4096 unique points over
68 lambda_p, nominal focus 40Tp, M=5/qx=5/Ta=0. Initialize eta/true-surface
psi with 11+20+22+33, zero muStar and linear frequencies. All eight final
phase runs reached 20Tp (240.20014 s); old directional OW3D jobs were not modified.

Read `research/unidirectional_time_series/hos_ocean/deepwater/RESULTS.md`.
Final full-domain L2 percentages at 3Tp/20Tp:
low second .009535/.026933, third 7.96591/4.94667;
high second .358064/1.032746, third 8.66296/6.74302.
Absolute tolerance 1e-10 was inadequate for low-amplitude harmonics; final
runs are tol16/akp002 and tol14/akp012 with checked tighter-tolerance differences.
Do not reuse the original coarse errors as physical results. The remaining
third-sector error is not proven to be startup error: four-phase -1/+3 aliasing
and evolving mixed-sign content remain unresolved; an eight-phase audit is a
possible next diagnostic, not already authorized or executed as a new campaign.

Run root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-deepwater-fourphase-20260926-v1`.
Local summaries/plots: `artifacts/hos_ocean/deepwater-results/`.
Per-phase final HOS wall times: low 214--219 s, high 72 s; RSS about 18.35 MiB.
MATLAB preparation/analysis peaks are about 1.6 GiB; see performance.csv.
All initial and tolerance-audit runs (28 HOS executions total) and logs remain
remote. No physical-accuracy, exact infinite-depth or long-time stability
certification is claimed. The older directional 220 s HOS proposal below
is separate and has not been launched.

## Precision update: 2026-09-25 22:35 UTC

The IO precision issue identified below is now resolved for external
`Initial surface quantities` and ordinary physical free-surface output.
Read `research/unidirectional_time_series/hos_ocean/PRECISION.md` first.

The upstream-v2.1.0 patch changes four IO lines only. Same-compiler baseline,
patched/legacy-input and patched/full-input 2 s pilots all passed. Raw t=0
relative L2 errors with full input: eta 8.65e-16, psi 8.04e-16. All eleven
frames were finite; same-input output differences satisfy legacy rounding
bounds. No fit or alignment was applied. The equations and zero-mode policy
are unchanged. Original binaries, runs and source snapshots are preserved.

Next external-surface solver:
`/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0/precision-v1/build-precise/sources/HOS-Ocean`.
SHA-256: `21a171f03c32a49ae05b62e4896d5b081c26a895372552d6ba563ba8fd44f9b4`.
Use the deployment's private library paths and one thread per initial job.
The pilots are complete; full 220 s HOS families have not been launched.
This is IO and short-run consistency validation, not physical-accuracy or
long-time stability certification. Earlier pending-precision instructions
below are historical and superseded by this update.

## HOS deployment update: 2026-09-25 22:18 UTC

HOS-Ocean v2.1.0 is now deployed remotely. The official startup example and
our 1024x256, M=5, qx=qy=3 wavegroup phase-0 pilot through 2 s completed.
The latter took 126.77 s wall time and 0.870 GiB peak RSS on one thread.
All eleven eta/psi frames were read and validated for finite values and times.
This is deployment/IO evidence, not physical-accuracy certification.

Read `research/unidirectional_time_series/hos_ocean/README.md` for hashes,
paths, scripts and the six-significant-digit input/output limitation.
The unmodified solver also resets the tiny mean modes introduced by input
quantization. No gain, bias or phase fitting was applied to saved data.
The full 220 s HOS campaign has not started. Preserve this pilot, improve and
verify IO precision before production high-order comparisons. The sixteen
OW3D processes were still running at the final check and were not modified.

The original handoff below predates this deployment. Its statements that
HOS has not been installed or run are historical; use the update above.

Current local project:
`C:/Users/spet5947/Documents/green-laplace-unidirectional-time-series`

Branch: `codex/unidirectional-time-series`.
Repository: <https://github.com/Xinyu01091002/green-laplace-wave-reconstruction>.
Remote alias: `60.188.112.99:60093`, user `lxy`.
Remote project: `/home/lxy/green-laplace-unidirectional-time-series`.

## OW3D background campaign

Launch verified: **2026-09-25 21:37:46 UTC**. Sixteen actual `ow3d` processes
were observed running, and the separate mail watcher was observed waiting.
Execution source commit: `1d4edd7d255b4ad7c50d2681db7ec9b8ee2450f5`.
This launch observation is not completion or propagation validation; consult
the live status files below. Later handoff-only commits do not change the run.

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
