# Hilbert four-phase spectrum audit, 2026-09-26

No HOS simulation was rerun. Existing final fields were processed in MATLAB.

Authoritative local implementation inspected:
`C:/Research/VWA/VWA time series/unidirectional/paperplot_VWA_2_4order_OW3D.m`,
coefficient matrix at lines 9--13, spatial quadrature at lines 217--218.
Temporal counterpart: `paperplot_VWA_time_series.m`, lines 214--229,
uses +imag(hilbert(time record)) instead of the spatial -imag(hilbert(field)).

The original matrix applied to [E,-Hx(E)] agrees with the previous spatial
positive-k implementation at all six saved snapshots. Maximum relative L2
of eta differences across harmonics: 7.053e-12; psi: 9.237e-12. Tiny DC/Nyquist
and roundoff terms account for the finite discrepancy; no first/third exchange
or material Hilbert-sign error was found. Earlier descriptions as ordinary
four-point DFT were incomplete: this is a Hilbert-quadrature four-phase method.
Do not use generic phase-DFT aliasing alone as a demonstrated cause of the
5--9% discrepancy or as evidence that extra phases are required.

A workflow difference remains: the original plot calls frequency_filtering
on the separated first/second/third harmonics (lines 50,53,54). The accompanying
implementation located under `VWA time series/directional/` uses smooth tanh
passbands nominally near .5--3, 1--4.5, and 2--5.5 kp. Its historical runtime
MATLAB path resolution has not been certified. The HOS comparison used raw
separated target fields and a fixed input-support projection instead, so it
has not reproduced the entire filtered legacy workflow. The effect of those
filters on reconstruction error has not been measured in this audit.

Delivered plot: `artifacts/hos_ocean/deepwater-results/fourphase-spectrum-audit/fourphase_spatial_spectra_final.png`
(and PDF). Upper row Akp=.02, lower .12; columns t=0,3Tp,20Tp.
This is a spatial wavenumber spectrum, not a measured temporal frequency
spectrum. Positive-mode amplitude is 2*abs(FFT(eta_h))/N in metres (DC and
Nyquist are not doubled). No bandpass or window, no individual curve
normalization, and no k-to-frequency dispersion remapping. Full positive
spectra are in six CSV files; the figure displays 0<=k/kp<=6 on a log amplitude
axis. Operator parity is recorded in operator-audit.json.

HOS runtime inventory (GNU time, excludes MATLAB and orchestration):

| Scope | Completed processes | Sum of process wall seconds | User + system CPU seconds |
| --- | ---: | ---: | ---: |
| Deep-water experiment including tolerance audits | 28 | 1629.02 | 1628.72 |
| Earlier deployment/directional/precision pilots | 5 | 539.29 | 539.19 |
| Total | 33 | 2168.31 | 2167.91 |

The deep-water HOS batches ran concurrently. Summing each actual batch's
slowest process gives approximately 355.54 s of active HOS waiting, excluding
MATLAB stages, scheduling overhead and gaps between batches; this is an
estimate, not the end-to-end campaign duration. Raw fields stayed remote.
See `hos-runtime-inventory.csv` and `hos-total-runtime.json` in the local
artifact folder. This audit added zero HOS runtime.