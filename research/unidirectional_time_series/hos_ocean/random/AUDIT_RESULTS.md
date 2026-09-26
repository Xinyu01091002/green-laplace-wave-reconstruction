# Random-wave record and directional-discretization audit

Completed 2026-09-26 using only the existing high focusing-label random HOS
records. No HOS rerun, no low family, no amplitude rescaling, and no new email.
Original center 7.5-degree prediction reproduced exactly (relative L2 zero).
All scores use the unchanged absolute interval 41.285991--178.714009 s.
Reference and prediction receive the same sum-band Fourier projection; raw
metrics are also retained. Harmonic sectors are not exact perturbation orders.

## Directional resolution

Full 220 s records, no taper; relative L2 against HOS second-harmonic sector:

| y offset (m) | 7.5 deg (%) | 3.75 deg (%) | 1.875 deg (%) |
| --- | ---: | ---: | ---: |
| 0 | 4.6191 | 2.9730 | 2.6827 |
| -52.7823 | 5.5281 | 4.2010 | 3.9632 |
| +52.7823 | 4.2278 | 2.8795 | 2.3946 |
| -70.3762 | 5.6951 | 4.7432 | 4.5277 |
| +70.3762 | 4.2981 | 3.0996 | 2.4183 |

The 3.75 -> 1.875 degree waveform changes are 1.30--2.16%, so the five-probe
calculation is not claimed fully converged. A further center-only 0.9375-degree
calculation gives error 2.6926%, versus 2.6827% at 1.875 degrees, and waveform
change 0.5369% relative to the 1.875-degree prediction. This supports practical
center stability with respect to angle, not a proof for all probes or all errors.

## Record-length and endpoint effects

Center, 1.875 degrees, same scoring interval:

| Record end (s) | Untapered error (%) | 10 s taper error (%) | 20 s taper error (%) |
| --- | ---: | ---: | ---: |
| 220 | 2.6827 | 4.5033 | 2.4017 |
| 210 | 3.9208 | 3.1576 | 3.0316 |
| 200 | 6.8972 | 2.5322 | 2.0540 |

The untapered predictions change by 4.2900% / 8.6971% on shortening to 210 /
200 s, compared with the full 220 s prediction. With a 20 s endpoint taper,
the corresponding changes relative to that taper's own full-record prediction
are 3.7368% / 2.2281%. A 10 s taper instead gives 5.6489% / 6.4052%; there is
no general guarantee that tapering improves this finite-record reconstruction.

The taper is a half cosine on each end: multiply first-order phase inputs and
initial directional prior by w, and the quadratic reference by w^2. The entire
scoring interval has w=1 for every tested record. This is a quadratic window
diagnostic, not an exact covariance of the frequency-dependent GL operator.
It must not be promoted as a fitted physical correction or best-error selection.

Keeping the full-record Hilbert separation then truncating the separated input
does not remove record sensitivity: the 210 s error is 2.8090%, but the 200 s
error is 9.0238%. Native FFT bin spacing, support projection, and directional
allocation remain coupled. Thus it would be incorrect to blame the Hilbert
transform alone or interpret the residual directly as intrinsic GL error.

## Practical interpretation

Coarse direction bins explain part of the earlier 4--6% discrepancy. At the
center the angle-refined score is about 2.7%, but temporal record processing
still shifts predictions by several percent. Keep these effects separate from
nonlinear physical/model error. Before treating normalized high/low random
results quantitatively, use a longer observation record with interior margins
and a declared consistent temporal analysis. These longer simulations are a
recommendation, not launched by this audit. One seed is not sea-state statistics.

## Evidence and resources

Remote root: hos-random-fourphase-20260926-v1 under the established runs root.
Successful main audit: record-direction-audit-v2; center refinement:
record-direction-center-finer-v1. Raw audit.mat files remain remote. Local plot:
artifacts/hos_ocean/random-audit/audit.png. Compact metrics accompany this note.
Only CPU40 was used; eight retained OW3D jobs were left running.

Main successful audit wall 79.45 s, peak RSS 1497944 KiB (1.43 GiB); center
refinement 26.28 s, peak RSS 2522664 KiB (2.41 GiB). The earlier v1 audit stopped
at 49.24 s because a proposed 30 s taper entered the fixed scoring interval
for the 200 s prefix. Its log is retained in v1. The successful audit uses
10/20 s tapers, both strictly outside the scored interval. No v1 partial result
is substituted for the completed v2 results. Source hashes are saved remotely.
