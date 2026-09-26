# HOS replication of the established GL time-series test

Completed 2026-09-26. Both four-phase families ran to 50Tp (688.09985438 s).
The selected 32--48Tp window contains 641 samples at Tp/40. This is the
practical near-matched replication requested by the user, not exact HOS/OW3D
solver parity. See README.md for physical settings and intentional differences.

## Main result

Raw main-group relative L2 errors in percent, against each solver's own
Hilbert-four-phase second harmonic. The group window is selected from the
shared first-harmonic input envelope, +/-2Tp, without fitting or alignment.

| Akp | Method | HOS reference | Existing OW3D reference |
| --- | --- | ---: | ---: |
| .02 | GL6 | 0.263498% | 0.259207% |
| .02 | GL12 | 0.108555% | 0.103794% |
| .02 | MF12 | 0.108186% | 0.103431% |
| .12 | GL6 | 3.973633% | 3.934372% |
| .12 | GL12 | 3.887767% | 3.848632% |
| .12 | MF12 | 3.895686% | 3.856581% |

HOS reproduces the reconstruction behavior seen in the established OW3D
test: very close low-steepness recovery and a few-percent discrepancy at
higher steepness shared by GL and MF12. This supports using HOS as the
working numerical reference for this single-direction eta22 time-series
workflow without waiting for new OW3D runs. It is not a validation of every
harmonic, depth, directional field or inverse reconstruction problem.

No GL formula, quadrature rule, or reference-led spectral cutoff was adjusted.
Both data sources retained 155 temporal parent bins under the original
kh>=.3 and twice-parent-frequency-below-Nyquist rules. The reconstructed
first-harmonic input was the observed time series, not the initial spatial
spectrum. First-input projection relative L2 is 1.6748e-4 / 2.1063e-4 for
low/high HOS, similar to the stored OW3D values.

The same processing function also reproduced all six saved OW3D predictions:
relative Frobenius differences 1.92894e-15 (low), 2.31816e-15 (high).
The MF12 linear synthesis agreed with the common first-harmonic input to
relative L2 <1e-11. Every saved HOS probe sample was finite; the physical
sample grid and required analysis count passed explicit checks.

## Filtering

Raw predictions and references both use no additional output bandpass.
A supplemental FFT sum-support mask was applied identically to the reference
and all predictions. The mask is determined only from retained input bins.
For HOS GL6 the main-group relative L2 fractions changed from
0.0026349796705 to 0.0026349796075 (low), and from 0.0397363347513 to
0.0397363320282 (high). Thus common output filtering does not materially
change this result. Filtered and raw full/main-window metrics are all saved.
This is not a claim to have duplicated every historical smooth filter.

## Measured time and memory

All production HOS processes used one thread. Four phases ran concurrently
within each low/high family; the families were processed sequentially.
Wall time and maximum RSS from GNU time:

| Stage | Low Akp=.02 | High Akp=.12 |
| --- | --- | --- |
| MF12 four-phase preparation, incl. MATLAB startup | 28.75 s, 1.593 GiB | 27.51 s, 1.613 GiB |
| HOS per phase through 50Tp | 292.35--299.45 s, <=18.406 MiB | 383.25--388.98 s, <=18.406 MiB |
| Analysis, old-OW3D replay, metrics and per-case plots | 40.41 s, 1.406 GiB | 43.32 s, 1.390 GiB |
| GL6 reconstruction call only | 0.84458 s | 0.81860 s |
| GL12 reconstruction call only | 0.67805 s | 0.68925 s |
| MF12 reconstruction call only | 0.49671 s | 0.51606 s |

The individual reconstruction call timings are single measurements, not a
cold/warm performance benchmark or a speedup claim.

This run added 9 HOS executions: 8 production phases and one 0.33 s IO smoke
case. Sum of process wall times: 2726.69 s = 45 min 26.69 s. Concurrent HOS
batch waiting is approximately 299.45+388.98+0.33 = 688.76 s = 11 min 28.76 s,
excluding preparation, analysis, build time, scheduling and operation gaps.
Including the previously recorded 33 HOS executions, cumulative process wall
time is 4895.00 s = 1 h 21 min 35 s across 42 HOS processes. These sums exclude
MATLAB and must not be mistaken for end-to-end elapsed project time.

## Differences and limits

HOS uses periodic boundaries, no31 initialization, and its native nearby
probe at x=10259.08756893 m. Existing OW3D uses walls, enabled third-order
primary corrections at initialization, and x=10259.18359375 m. Their spatial
grids/domain lengths also differ slightly. The user explicitly accepted
close settings rather than exact replication. Weak/free-wave boundary
influence was not separately isolated; this test supports the reconstruction
comparison, not equality of the two boundary-value problems.

A single relative integration tolerance 1e-10 was used. The earlier spatial
study motivated adequate precision; no new tolerance sweep was performed.
The short native-field/probe check verified eta at a grid point to about
1.2e-15 m. The probe reader's header was corrected from ZONE to VARIABLES
before production, with the failed orchestration record preserved.

## Artifacts

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1`.
Raw probes, source snapshots, native smoke fields and comparison MAT arrays
remain remote. Local compact results:
`artifacts/hos_ocean/timeseries-results/`.

- timeseries_summary_final.png / .pdf: main-group comparison, both solvers.
- akp002/metrics.csv and akp012/metrics.csv: all methods, windows and filters.
- akp002/report.json and akp012/report.json: settings, method timings and parity.
- performance.csv and summary.json: stage RSS/time and HOS runtime totals.
- build.json and hos-io.patch: exact probe-capable executable provenance.

Code and results were added locally; no Git commit or push was performed.