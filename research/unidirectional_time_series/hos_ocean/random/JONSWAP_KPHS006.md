# JONSWAP kp Hs / 2 = 0.06 continuation

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
