# Sequential phase-only random HOS trials

User approved starting on 2026-09-26 after selecting the recommended
phase-only comparison. No Hs renormalization and no spatial taper.
High/low labels .12/.02 refer to potential focusing steepness, not kp Hs/2.
Seed 20260925 is shared across amplitude levels and across all four globally
shifted phases (0/90/180/270 degrees).

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-random-fourphase-20260926-v1`

Screen: `hos-random-sequential-20260926`. `pipeline-status.json` tracks the
whole queue. Per-family `akp012/status.json` and `akp002/status.json` track
HOS. Order: high HOS, high GL, low fresh initialization/HOS, low GL.
GL failure is recorded and does not discard the completed HOS data or block
the next amplitude; HOS failure stops the queue for investigation.

High HOS started 2026-09-26 15:28:08 UTC. All four initial-probe checks passed
(max absolute error below 2e-15 m); all advanced beyond t=0. At 15:28:28 UTC,
the four log snapshots had reached 0.4 s. Sampled peak aggregate RSS was
4448984 KiB, about 4.24 GiB. No completed random propagation/result is claimed
at this launch checkpoint. Low is queued, not yet running.

## Settings and resources

- kp=.0279 /m; kph=1; g=9.81; Tp=13.7619971 s.
- 50 x 20 wavelengths, 11260.1887 x 4504.0755 m; unique grid 1024 x 256.
- Mean direction +x; Gaussian amplitude angular sigma 25 degrees; unchanged
  radial support and modal magnitudes from the focused families.
- M5, qx=qy=3; adaptive Cash-Karp RK5(4), absolute tolerance 1e-12; Ta=0.
- Duration 220 s (15.986 Tp); probe interval .2 s, 1101 samples; five probes.
- MF12 order2 eta and true-surface psi, 11+20+22; no initial31/33.
- Four phases x eight MPI ranks on disjoint CPU sets 8--39. MATLAB on CPU40.
  Retain the eight OW3D jobs. No new MPI benchmark or GPU task.
- Aggregate HOS sampled RSS guard 8 GiB; wall/RSS records saved for all ranks,
  controllers, initialization and GL. Launch check: 48 CPUs, 605 GiB available;
  after HOS startup about 601 GiB remained available.
- High linear Hs=.4076151725 m, kp Hs/2=.005686231656. Low expected Hs=
  .0679358621 m, kp Hs/2=.000947705276. Low is freshly constructed with MF12;
  assert its first-order complex C equals high C/6, including random phases.

## GL comparison

Reuse the existing GL J8 implementation and initial-complex-directional-shape
allocation, with 15/7.5-degree discretizations and observed separated first
harmonic. The target is the second-harmonic phase sector. No exact perturbation
order claim, no fitting, no inversion of directions from a single probe.

Replace the focused envelope-peak window with the predeclared fixed interior
interval [3 Tp, 220-3 Tp] = [41.2860,178.7140] s. Full-record raw metrics remain
available. Reference and predictions get identical frequency projection.
The same fixed interval is used for the 10 s shorter center-record diagnostic.
It measures record-end/frequency-grid sensitivity; it does not establish that
initial transients have vanished. Directional cancellation/conditioning is
reported, with no numerical floor or empirical amplitude correction.
The existing one-third reference-amplitude eligibility rule is retained.

Outputs: each family's `random-gl-comparison-v1`, including metrics.csv,
report.json, full_comparison.png/pdf and the remote comparison.mat.

## Frozen source

`package_source.py` creates the deployment snapshot from the committed focused
workflow plus these new queue/preparation scripts. Local snapshot:
`artifacts/hos_ocean/random-launch-source-v2/`, with SHA256 source-manifest.json
and base commit c27f4c650b0c6a0959afe4dc85cd820a706bacd9. The first local assembly
attempt hit a Windows default-text-encoding assertion before deployment; v2
uses explicit UTF-8 reads. Only v2 was deployed and launched.
Raw data remain remote. No existing source/run directory was overwritten.

## User update: high-only and completion email

The user cancelled the queued low-amplitude run and explicitly requested email.
Only the active Akp=.12 focusing-label random case and its GL comparison remain.
Remote cancel-low.json guards prepare_random_run before any low initialization.
The already-running legacy Python loop is not reloaded by edits: on reaching
low preparation it raises USER_CANCELLED_LOW_BEFORE_INITIALIZATION. A detached
finish_high_only.py watcher reconciles that intentional cancellation, preserves
high results, writes summary/finished markers, and invokes completion_email.sh.
It reports success only if high HOS and GL both completed; other outcomes are
reported as needing attention. No active solver or OW3D process is stopped.
The live guard/notification sources are recorded by SHA256 in cancel-low.json.
Email reuses the remote /usr/local/bin/run.sh recipient without copying it into
source control. One success/failure notice; local mail acceptance is not proof
of inbox delivery. New Hs-normalized random campaigns have not been launched.
