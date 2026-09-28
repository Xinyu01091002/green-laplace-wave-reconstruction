# Existing unidirectional HOS eta33 check, 2026-09-27

The single-direction space-time GL eta33 trial was checked against the
existing Akp=.02 and .12 HOS four-phase records. After extracting harmonics
on the complete available 50Tp record and cropping to the unchanged
32--48Tp comparison interval, the FFT trial's relative L2 errors are:

| Akp | Main-group error against HOS | Full 16Tp comparison error | Error against analytic GL reference |
| ---: | ---: | ---: | ---: |
| .02 | 0.32621% | 0.60925% | 0.24940% |
| .12 | 8.63221% | 8.64640% | 0.27662% |

These results use every nonzero temporal input bin through 4omega_p, no
old kh cutoff, L/h=32pi, and Nx=4096. The analytic GL reference on identical
inputs has HOS main/full errors of 0.31727%/0.59003% and 8.72908%/8.74195%.
Thus the large-steepness discrepancy is not resolved by this numerical
executor. HOS phase sectors remain distinct from pure perturbation orders.
There was no new HOS integration, MF12 third-order calculation, taper,
gain/phase/time fitting, or modification of an existing production run.

## Data and execution provenance

Existing source:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1`.
The preserved setup is kph=1, kp=.0279/m, h=35.8422939068 m,
Tp=13.7619970876 s, Alpha=1, four global phases, Akp=.02/.12, 50Tp integration,
and Tp/40 output. The source records have 2001 samples; the comparison
interval contains 641 samples. The original HOS setup and its differences
from OW3D are documented in `hos_ocean/timeseries/README.md` and `RESULTS.md`.

Isolated postprocessing run:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-eta33-resolvent-20260927-v1`.
Execution used MATLAB R2026a as lxy, with `-singleCompThread`. The initial
check completed with exit 0 at 18:54:53 UTC; the separate extraction-window
diagnostic completed with exit 0 at 19:02:25 UTC. Timings from this remote
single-thread check are not interchangeable with earlier local CPU scaling.

Deployment uses a SHA-256 source manifest. The source `comparison.mat`
hashes are:

- Akp=.02: `57a854bd65d8a6d951075853cddfddfba71e3e8fe8193f9e438de7a277e20c61`.
- Akp=.12: `bb2b660df2148d03bac7d5d30b8dc802055a1334f8a3716359993346926023e3`.

All eight full `probes.dat` hashes are retained in `window_status.json`.
Raw probes and numerical field MAT files remain remote. Only compact
metrics, summaries, status/provenance and figures were collected locally.

## Initial check using the old 16Tp extraction window

`check_hos_eta33_resolvent.m` reproduced the saved first-phase extraction,
then extracted the third phase sector with the established Hilbert/four-phase
sign convention. It evaluated two explicitly separate input supports:

1. The historical cubic rule kh>.5 and all triple frequencies below the
   native Nyquist: 98 parents. This preserves the old comparison definition.
2. All nonzero positive temporal bins through 4omega_p: 64 parents. There
   is no kh or initial-spectrum min/max cutoff. Cubic output remains within
   the native temporal band.

The main-group window is fixed from the existing first-harmonic envelope
peak +/-2Tp and reused for both supports and the later diagnostic. At
L/h=100.531, the FFT executor's relative complex errors against analytic GL
are 0.24641%/0.27229% for historical low/high inputs and 0.25034%/0.27683%
for full-band low/high inputs. Smaller domains were also retained.

The old short-window extraction gave:

| Akp | Support | Analytic GL vs HOS, main | FFT vs HOS, main | FFT vs HOS, full |
| ---: | --- | ---: | ---: | ---: |
| .02 | historical cubic | 0.52009% | 0.37369% | 34.7174% |
| .02 | through 4omega_p | 0.39713% | 0.40219% | 34.7174% |
| .12 | historical cubic | 11.7630% | 11.6020% | 11.6620% |
| .12 | through 4omega_p | 8.73290% | 8.63537% | 8.71760% |

The unexpectedly large low-steepness full-record residual required a
separate extraction-window check; it was not discarded or hidden by
reporting only the main group. The original metrics and figures remain
in the isolated run's `results/` directory.

## Extraction-window diagnostic

`check_hos_eta33_record_window.m` reads the original eight full HOS probe
files and performs exactly the same Hilbert/four-phase extraction on all
2001 samples first. It then crops to the same 641 comparison times and
recomputes the first-order input and GL predictions through 4omega_p.
The comparison interval and main-group mask do not change. No taper,
demeaning, output-band filter, fitted alignment or amplitude correction
is introduced.

For Akp=.02, 99.9885% of the old analytic-reference residual energy was
outside the main-group window. The longer extraction changes the first
harmonic by only 0.020701% relative L2, but changes the very small third
phase sector by 34.7087%. The analytic GL full error drops from 34.7177%
to 0.59003%, and the FFT trial full error becomes 0.60925%. The figures
localize the removed discrepancy to the ends of the cropped record.
This supports finite-window harmonic-extraction effects as the dominant
cause of that low-steepness full-record discrepancy.

For Akp=.12 the first-harmonic change is 0.021674%, the third-sector change
is 1.08002%, and only 1.8237% of the old residual energy is outside the
main group. The analytic full error changes from 8.81301% to 8.74195%.
Consequently the high-steepness residual is materially different from the
low-steepness endpoint issue; no further cause is assigned without evidence.

The complete-record extraction is the better supported route for continuing
this HOS comparison. It does not establish that 50Tp eliminates every
finite-record effect, nor that a generic short measured record has the same
available context. Keep the initial and extended-extraction results separate.

## Artifacts and continuation

Local compact results:
`artifacts/unidirectional_time_series/hos-eta33-resolvent-20260927-v1/`.
Its `results/metrics.csv` and `results/summary.json` describe the initial
check. `window_diagnostic/summary.json` and `akp002.png/.pdf`, `akp012.png/.pdf`
describe the complete-record extraction. `manifest.json`, `window_manifest.json`,
status files and the executed window-script snapshot retain provenance.
The local window-script change after execution only annotates a false
allocation warning for the already preallocated four-column raw array.

The five new MATLAB files across this HOS check and the directional trial
passed Code Analyzer. The three deployment/collection Python scripts passed
syntax compilation. Source mathematical identities are covered by the prior
cubic checks and the new five directional identities.

Directional continuation is documented in
`../directional_wave_data/DIRECTIONAL_TIME_FFT_20260927.md`: a small
frequency-direction fixture converges on the tested larger domains. The
subsequent real directional HOS trial is now documented in
`../directional_wave_data/DIRECTIONAL_HOS_ETA33_TRIAL_20260927.md`; it did
not pass validation. The single-direction results in this report are separate.
