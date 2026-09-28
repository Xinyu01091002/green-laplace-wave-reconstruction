# JONSWAP kp Hs / 2 = 0.06 continuation

## Verified completion, 2026-09-28

The **HOS 20Tp family completed**. All four phases exited0 with1377 finite
six-column probe records from0 to275.2s at.2s spacing. A fresh read of each
raw `probes.dat` confirms those conditions, equality with the corresponding
postprocessing CSV, and agreement of each CSV SHA256 with `snapshot.json`.
No HOS processes from this run remain. This verifies native completion and
record integrity, not independent spatial/time convergence or physical
accuracy of the nonlinear simulation.

| Phase | HOS finished UTC, 2026-09-28 | Wall seconds | Samples | Exit |
| --- | --- | ---: | ---: | ---: |
| 0 | 05:41:28 | 45011.86 | 1377 | 0 |
| 90 | 05:40:57 | 44981.05 | 1377 | 0 |
| 180 | 05:57:19 | 45963.22 | 1377 | 0 |
| 270 | 05:58:40 | 46044.31 | 1377 | 0 |

The four-phase wall time was about12h47m24s; aggregate HOS peak sampled RSS
was8102964KiB (about7.73GiB). Frozen execution source commit:
`802d1f67d931b92aa4a99e6f7fc3a4c7fa774e74`.
The previously queued postprocessor finished at06:00:48 UTC, exit0, and the
whole pipeline reports `completed_HOS_and_GL`. Mail was accepted by the local
MTA at06:00:51 UTC (07:00:51 British Summer Time); inbox delivery is not
independently confirmed. No additional email or numerical job was started
in this completion audit. Raw and derived time records remain remote.

### The automatic GL result is historical, not the new method

The old queue automatically called its frozen
`source/research/directional_wave_data/gl_directional_sum_time.m`, which
explicitly forms ordered-pair arrays. This is the legacy implementation now
prohibited by the latest user instructions in
`HANDOFF_TIME_SERIES_GL_20260928.md`. The completion audit only read existing
reports/source and checked files; it did not rerun that executor or load its
prediction fields. Do not use the successful pipeline label as evidence
that the requested new low-cost direct-time GL derivation is complete, and
do not rerun this historical queue/comparison as a new benchmark.

For provenance, the saved legacy report at
`medium/gl-time-comparison-v1/report.json` has SHA256
`86feb1cd79043d4d75050c17876e64fe9c5c932088ba6a8c3cdcb714759882d9`.
Its fixed3--17Tp-window relative L2 values are recorded below, without new
evaluation or tuning:

| y offset (m) | Legacy raw L2 (%) | Legacy common sum-band L2 (%) |
| --- | ---: | ---: |
| 0 | 26.7147 | 26.7304 |
| -52.7821 | 11.4855 | 11.5235 |
| +52.7821 | 11.5059 | 11.4948 |
| -70.3762 | 9.9023 | 9.9621 |
| +70.3762 | 11.2850 | 11.2388 |

All five probes pass the former one-third reference-amplitude eligibility
rule (ratios.761--1.129 off center). Removing the final10s changes the
center legacy prediction by55.6916% relative L2 in the same interior window,
while its first input changes4.3843% and filtered reference2.3608%.
That is sensitivity to record length, not a55.7% GL--HOS error measurement.
The center uses38 retained frequencies; positive-frequency input energy
retention across probes is99.4512--99.7327%, below the user's newer99.9%
full-record requirement. These limitations preclude treating this historical
comparison as an accepted reconstruction validation. The HOS records remain
available as observations for the new derivation and subsequent validation.

## Run design: first 20Tp, requested 2026-09-27

The user shortened this trial from approximately80Tp to approximately20Tp.
The completed run root is
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-kphs006-20tp-20260927-v1`.
Screen: `hos-jonswap-kphs006-20tp`. HOS started **17:11:16 UTC** on2026-09-27.
Read `pipeline-status.json` and `medium/status.json` there for terminal state.

The new native solver endpoint is **275.2s**, 19.9971Tp, with the same.2s
output interval and **1377 samples** per probe. Fixed comparison window:
**3--17Tp = [41.2859913,233.9539505]s**. Full-record raw metrics and the
existing10s-shorter-record sensitivity check remain. This is a shorter
initial trial; the window does not certify that startup transients vanish.
Its temporal frequency-bin spacing is about four times that of the80Tp
design. No kernel, filter, amplitude, phase, direction or grid was changed.

The deployed HOS reads the duration at startup; its inspected evolution loop
has no stop-time reload. The former80Tp parent/controller and only their four
MPI jobs were therefore stopped at17:10:51 UTC, before a fresh20Tp launch.
Their last saved/logged times were8.2/8.2/8.0/8.0s. There was no HOS failure.
Old inputs, logs and partial outputs remain in the80Tp directory, with
`duration-change-before.json` preserving the prior status/process snapshot.
Its terminal status is `stopped_for_user_duration_change`; it sends no
failure email for this intentional cancellation. OW3D and other work are
untouched. The old record has no restart state; the new record starts at0s.

`duration-input-check.json` verifies that all32 new eta/psi input slabs and
all four probe-position files are byte-identical to the former80Tp inputs.
Only duration, sample count, requested-duration metadata and scoring-window
settings differ. All four new YAML inputs specify275.2s. Code Analyzer,
declared-first-order checks, MPI export equality and live initial probes pass.
The unchanged HOS controller expects1377 finite records and native exit0;
it will not mistake an interrupted partial record for completed20Tp.
GL postprocessing and the existing one-shot completion/failure email follow.

At17:13:58 UTC all four new phases had finite records through.6s and eight
live ranks each; peak aggregate RSS was7.7235GiB. No old80Tp HOS ranks
remained, all eight OW3D processes were still present, and the intentional
80Tp cancellation had sent no failure email.

Frozen local snapshot:
`artifacts/hos_ocean/jonswap-kphs006-20tp-20260927-v1/`.
`launch-source-manifest.json` records base commit
`fe8d3272979f255be9309ca61e189702642b45f8`, both source hashes and the20Tp
request. `previous-retirement.json` records the intentional stop, and
`retire_previous_v2.py` is retained with its hash in the old before-snapshot.
Actual launch in the new root:

```sh
python3 retire_previous_v2.py
screen -dmS hos-jonswap-kphs006-20tp bash -c \
  'exec python3 -u run_jonswap_rescaled.py > runtime.log 2>&1'
```

## Historical 80Tp launch, superseded by the run above

User requested this amplitude on 2026-09-27 after the .12 family failed.
This is a new four-phase family, preserving the previous random realization,
spectral support, depth, grid and numerical parameters. The old .12 failure
and the completed .02 family are retained. Neither family is restarted.

Remote run:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-kphs006-80tp-20260927-v2`.

Screen: `hos-jonswap-kphs006-80tp-v2`.
Whole workflow: `pipeline-status.json`; HOS: `medium/status.json`.
HOS started **2026-09-27 16:47:28 UTC**, four phases with eight MPI ranks each,
on CPUs 8--39. MATLAB preparation and subsequent GL use CPU 40. The eight
existing OW3D processes are untouched; all new raw data and processing stay
remote. No new HOS evolution code, filter or slope threshold was introduced.

At 16:50:22 UTC all four phases had saved finite probe records through .8 s,
five samples each, with eight live solver ranks per phase and no reported
failure. Aggregate peak RSS was 8094204 KiB (approximately 7.72 GiB).
This is verified nonzero advancement, not passage beyond the old 5.7 s
failure time or completion of the requested record.

## Initial conditions and consistency

- kp=.0279 /m, h=35.8422939068 m, kph=1, Tp=13.7619970876 s.
- **kp Hs / 2=.06**, with linear spatial-variance Hs=4.30107526882 m.
  This is not the older potential-focusing Akp label.
- Same 11128 retained parents, JONSWAP gamma=3.3, seed=20260925, complex phases
  and directional distribution as the .12 input. Four global phase shifts
  are 0, 90, 180 and 270 degrees.
- Retain initialization 11 + R4 20 + GL16 22, for eta and true surface psi;
  no initial 31/33. First order scales by 1/2, second order by 1/4.
  The total nonlinear field is not simply divided by two.
- MATLAB reconstructs first-order eta/psi directly from declared C/kx/ky/om,
  checks them against the source four-phase odd sector, then scales the
  frozen order-two remainder. The reused fields are generated initial
  conditions, not evolved HOS reference fields.
- Relative first-order recovery errors are at most 1.01e-16; recovery of
  the separately saved .02 inputs by the same homogeneity is within 1.71e-16.
  Hs, phase identities, finite fields and zero spatial mean pass.
- All 32 exported MPI slabs match the new saved eta/psi arrays with maximum
  absolute error exactly zero. Four live initial-probe checks pass below
  4e-14 m. These checks establish implementation consistency, not long-time
  stability or independent physical accuracy.

| Phase | Initial maximum spatial wave height, m | Maximum crest on native grid, m | Maximum absolute x slope |
| --- | ---: | ---: | ---: |
| 0 | 8.9433 | 5.2734 | .19259 |
| 90 | 8.5369 | 5.1389 | .23084 |
| 180 | 8.7242 | 6.1183 | .21912 |
| 270 | 9.7729 | 5.3977 | .20639 |

Wave height is crest minus trough between adjacent zero-upcrossings along
periodic +x transects, with fourfold x Fourier interpolation. It describes
the initial spatial field, not the maximum temporal wave over the run.
The largest native-grid wave is 9.7572 m. Slopes use Fourier differentiation
and twofold interpolation in x/y. Eta second/first-order L2 ratios are
8.00--8.06%. Full metrics are in `initialization-audit.json`.

## Computation and completion criteria

Same 50x20 peak-wavelength periodic box, unique grid 1024x512, M5/qx=qy=3,
adaptive Cash--Karp 5(4), absolute tolerance 1e-12, Ta=0; no dissipation or
breaking model. Duration 1100.8 s (approximately 80Tp), output every .2 s,
5505 samples per probe. Same five probes. Initial full eta/psi are retained;
production outputs probe eta only, with no dense evolving full-field output.

One continuous production run is used; there is no separate dt study or
duplicate short run. `past_10s` becomes true only when every phase has logged
at least 10 physical seconds. Passing it would exceed the .12 failure time,
but would not establish full-record stability.

The unchanged HOS controller verifies initial probes, binary hash, finite
records, all expected sample times and zero exit codes. Aggregate HOS RSS
guard is 16 GiB; launch peak was approximately 7.6 GiB, with approximately
605 GiB host memory available before launch. Preparation and GL guards are
32/128 GiB, with minimum host available memory 64 GiB. No completion-time
forecast is inferred from startup.

After all four HOS phases complete, the frozen GL time-series comparison
runs under `medium/gl-time-comparison-v1`, using GL16, 3.75/1.875-degree
direction bins, identical reference/prediction output bands, fixed 10--70Tp
scoring window, and the existing one-third amplitude eligibility rule.
No MF12 eta33 comparison, fitted amplitude, phase or time shift is added.
The script then sends one completion/failure email through the existing
remote recipient configuration. `accepted_by_local_mail` means local MTA
acceptance, not confirmed inbox delivery.

## Source and actual launch

New maintained files: `prepare_jonswap_rescaled.m`, `run_jonswap_rescaled.py`.
Local deployment snapshot is in ignored
`artifacts/hos_ocean/jonswap-kphs006-20260927-v2/`.
`launch-source-manifest.json` records base commit
`f5b0a57f2f37677f182f99c794d7270fb1bb0c5d` and exact source SHA256 values.
The remote `inherited-source-manifest.json` separately records old inputs,
controller, export, GL code, dependency archive and email adaptation.
The HOS binary SHA256 is unchanged:
`93d64252f732449fad9d2171979c91e13b0b9de1d17a73d7587fa00721efcbb6`.

The target directory was checked absent, then created explicitly. The two
new scripts and manifest were uploaded as a tar snapshot; hashes and Python
syntax passed before execution. The actual launch in the new run directory:

```sh
screen -DmS hos-jonswap-kphs006-80tp-v2 bash -c \
  'exec python3 -u run_jonswap_rescaled.py > runtime.log 2>&1'
```

MATLAB Code Analyzer and all initial-field assertions passed before HOS.
The preceding `...-v1` attempt stopped before initial-field generation or
HOS on two Code Analyzer array-growth notices in the statistics helper.
The helper was corrected; v1 and its preparation-failure mail record remain
untouched. That event is not a .06 hydrodynamic failure.

For current state, read both JSON status files and the live processes/logs;
do not infer current progress or completion from this launch checkpoint.
