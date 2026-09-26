# Akp=.02 directional HOS run and automatic GL comparison

Completed: see [COMPLETED.md](COMPLETED.md). The launch notes below are historical.

User approved retaining one OW3D wavegroup and one random realization, then
performing a SHORT 4/8-process check before launching the lower-steepness
HOS family. No full 220 s benchmark was run at both process counts.

## OW3D resource change

Eight retained cases are wavegroup_phi000/090/180/270 and
random_s20260925_phi000/090/180/270. The eight seed-20260926/20260927 cases
were sent SIGTERM only after their PIDs, executable, parent and working
directories were checked. All files were preserved. Retained processes were
confirmed live; no controller-wide termination or restart was performed.
Available RAM rose to approximately 604 GiB.

The frozen 16-case controller labels native exit -15 as failed. These eight
are intentional user cancellations, documented in
`ow3d16-kpd1-akp012-20260925T212822Z/user-requested-reduction-20260926.json`
and each stopped case's user-cancelled.json. Do not interpret them as numerical
failures or automatically restart them. The first 3-second process-exit check
was too early during memory release; subsequent verification confirmed all
eight selected processes exited without SIGKILL and all eight retained jobs
remained live.

## Fresh low-steepness initial fields

The frozen geometry/spectrum is reused, but MF12 second-order coefficients
and all four eta/true-surface-psi initial fields are recomputed for Akp=.02.
The total Akp=.12 field was not divided by six. An independent homogeneity
audit against separately scaled linear and quadratic components gives
relative L2 1.57e-15 (eta) and 1.52e-15 (psi).

Unchanged physical/numerical settings: kph=1, kp=.0279, g=9.81,
50x20 lambda_p, 1024x256 unique points, amplitude angular sigma 25 degrees,
nominal focus 110 s, M5/qx=qy=3, Ta=0, absolute tolerance 1e-12,
220 s duration (15.986Tp), .2 s probe output, five original probe locations.
Initial content remains MF12 11+20+22, without initial 31/33.

## Step 3 was strictly a short benchmark

Same phase-0 low-steepness input and executable, each run only 2 physical
seconds (11 output frames at .2 s), sequentially:

| MPI ranks | Wall time | Peak sampled job RSS |
| --- | ---: | ---: |
| 4 | 51.4794 s | 1016748 KiB |
| 8 | 29.5748 s | 1117236 KiB |

The 8-rank case reduces measured wait by about 43% (1.74x speedup relative
to 4). All frames and probes passed the planned 4/8 equivalence checks;
maximum eta/psi relative L2 differences are 9.19e-11 / 7.11e-11. These are
short-run measurements, not a full-campaign speed guarantee. Selection is
based on wall-clock time and resource budget, not on reducing CPU-hours.

## Active production and automatic processing

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1`.
Screen session: `hos-akp002-auto-gl-20260926`.
Production started 2026-09-26 around 13:10:33 UTC.

Only the selected 8-rank configuration was launched for 220 s. Four phases
run concurrently: 32 HOS compute ranks plus 8 retained OW3D jobs, nominally
40 of 48 CPU slots. Rank affinity uses CPU IDs 8--39 in disjoint phase sets;
CPU 40 is reserved for sequential MATLAB processing. HOS aggregate sampled
RSS was about 4.25 GiB shortly after launch; the 8 GiB job-RSS guard remains.
All four initial probe checks passed and nonzero advancement was observed.

`status.json` describes HOS processes and progress. `pipeline-status.json`
describes the whole task. Once all four HOS cases exit zero and their 1101
records are validated, the pipeline automatically freezes complete records
and runs the directional GL time-series comparison on one CPU. It then runs
the same full-record processing for the already completed Akp=.12 dataset,
so both use the same rules. Outputs are reserved under each run's
`full-gl-comparison-20260926-v1`; do not independently create those folders
while the pipeline is active.

The GL path retains original J8 quadrature, 15/7.5-degree directional
allocation, observed first-harmonic input, and the input-defined shared
frequency support. Both raw and identically filtered prediction/reference
metrics are saved, with the one-third off-centerline amplitude eligibility
rule. No gain/phase/offset fitting is introduced. The complete low-steepness
GL result is pending HOS completion, not claimed already finished.

## Completed Akp=.12 timing

Fresh native-data checks confirmed all four .12 cases have 1101 records
ending at physical time 220 s. Per-phase maximum rank wall times:
0 deg 2:23:53; 90 deg 2:25:01; 180 deg 2:33:16; 270 deg 2:32:31.
The concurrent batch ran 04:12:14--06:45:33 UTC: 2:33:19 actual waiting time.
The original status last_logged_time values can lag after completion;
final probe records and completion states are authoritative.

Local compact evidence is in `artifacts/hos_ocean/akp002-run/`.
Inputs, raw outputs, benchmarks and resource files remain remote. No GPU
execution is used: the released resources and MPI allocation are CPU/RAM.