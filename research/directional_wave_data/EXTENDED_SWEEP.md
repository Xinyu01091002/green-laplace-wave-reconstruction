# Farther probes, depth, bandwidth, steepness and directional eta20

This user-requested extension freezes the candidate formula, direction
allocation, 15/7.5-degree angular grids, eta22 GL8 and psi22 GL2+2 before
reading high-order scores. The shared-scale eta20 diagnostics use ranks
6/12/16. No high-order target selects a location, rank, formula or weighting.
Public executors and historical source data remain unchanged.

## Existing OW3D cases

The local `VWA Unidirectinal/Directional` source provides:

| Family | kh | spread label (degrees) | Akp | Scope |
|---|---:|---:|---:|---|
| test1 | 1 | 25 | .02 | center; x offsets +/-3 lambda_p; y offsets +/-1 and +/-1.5 lambda_p |
| test1 | 2 | 25 | .02 | center |
| test1 | 5 | 25 | .02 | center |
| test1 | 1 | 5 | .02 | center; narrower *angular* spread |
| test6 | 1 | 25 | .12 | center |
| test1 | 5 | 15 | .02 | center; paired angular family for deep steep case |
| test6 | 5 | 15 | .12 | center |

This is 13 predeclared probe/case records. Grid-nearest physical coordinates
are stored, not interpolated. All test1 files at steps 0:10:900 and test6
files at 0:10:1200 are read and hashed, giving 4 s sampling. The source grids,
domains and record lengths differ between test1 and test6, so the comparison
is not claimed as an otherwise-identical single-parameter numerical experiment.
Each probe uses its own eta1 record; the prior directional structure uses
only that case's initial first phase-sector spectrum. Every family has an
independent initial eta/psi polarization audit before the global convention
conversion. Failures are recorded rather than silently removed or retuned.

Actual relative frequency bandwidth is measured from the initial complex
spatial coefficients as energy-weighted std(omega)/mean(omega). The directory
spread is not renamed frequency bandwidth. The inspected historical folders
do not establish a matched, independently varied radial-bandwidth OW3D family.
Other local directionally spread folders with explicit Alpha labels either
have only Alpha=8, confounded parameter changes, or lack four phases; no claim
of an isolated OW3D frequency-bandwidth validation is made.

## Directional eta20

The existing shared-scale GL graph is transferred with full vector dot
products and |k_i-k_j|. Its scale uses vector wave-number variance weighted
by the declared joint-input complex amplitudes. The published physical
assembly factor is retained. Strict K=0 pairs are excluded, but omega_i=omega_j
with different directions and nonzero K are **retained** in the raw prediction.
Their stationary contribution is separately recorded in each audit.

The observed phase-zero record is the average of four phases and may include
higher harmonics. Two explicitly oscillatory comparisons project *both*
prediction and reference onto 0<|omega|<0.5 omega_p and 0<|omega|<3 omega_p.
On these 4 s records, 3 omega_p exceeds the sampled Nyquist, so the latter
retains all available nonzero temporal FFT bins. It cannot recover unresolved
high frequencies or guarantee removal of aliased fourth harmonics. Removing
temporal DC also removes legitimate nonzero-spatial-K stationary contributions;
this comparison is therefore not a complete eta20 or mean-flow certification.
Raw candidate/reference records are preserved and no offset is fitted.

Wolfram checks the vector difference forcing, pair symmetry and the distinction
between zero temporal frequency and nonzero spatial K. MATLAB compares all
three eta20 ranks against the original spatial diagnostic, including equal-
frequency/different-direction pairs: max absolute difference 1.0842e-19 m.

## Bounded computational implementation

For uniform temporal bins, pair coefficients are accumulated at sum/difference
indices and synthesized by FFT. Integer wrap reproduces the *same discrete
sample values* as direct exponentials, including output frequencies above the
sampled Nyquist. It is not an alias-free continuous-time claim or padding.
A separate test compares direct and accumulated eta22/psi22/eta20, including
above-Nyquist sums: max relative difference 1.0801e-14. This avoids repeating
large exponential matrices for every probe; it changes no physical kernel.

## Frequency-bandwidth benchmark (synthetic, separate)

To vary true frequency bandwidth independently, six synthetic cases combine
kh=1/5 with nominal Gaussian energy widths beta=.08/.16/.32. Frequencies are
fixed temporal bins with omega/omega_p=.4:.05:1.6; seven known directions span
-60 to +60 degrees with a prescribed Gaussian angular amplitude weight.
The true directional phases are specified, not inferred from reference fields.
A fixed, frequency-dependent amplitude/phase multiplier changes the initial
record while preserving the relative direction composition. The joint-input
constraint must recover those known evolved coefficients exactly.

Compare GL eta22/psi22 and nonzero-spatial-K eta20 against independent spectral
MF12 **order two only**, on identical known inputs and in physical units. This
tests the assumed allocation and GL kernel approximation across bandwidth;
it is not OW3D evidence, a simulation of nonlinear evolution, or proof that
real nonlinear waves preserve their directional composition. The measured
post-modulation relative frequency widths, rather than nominal beta alone,
are reported. No candidate is tuned or chosen from these reference values.

## Reproduction

```matlab
addpath('research/directional_wave_data');
test_directional_difference_time;
run_directional_extended_sweep('C:/Research/VWA/VWA Unidirectinal/Directional');
run_directional_bandwidth_benchmark('C:/Research/spectral domain implementation of wave interaction theory');
```

Symbolic check: `verify_directional_difference.wl`.
OW3D results: `results/directional_sweep/`; synthetic results:
`results/directional_bandwidth/`. Per-probe fields include raw eta20, stationary
terms, both oscillatory bands, full and main-group metrics, conditioning and
window-completeness status. All plots use only saved OW3D samples.

## Executed OW3D outcomes

All 13 planned probe/case records completed with finite outputs, complete
input-defined main windows and no execution failures. That does not mean all
accuracy tests passed. The primary table keeps the predeclared 7.5-degree
direction grid; eta20 below means rank 16, oscillatory band 3*omega_p capped
by the available sampled Nyquist, excluding temporal DC from *both* sides.
Numbers are main-group relative L2 in percent.

| Center case (kh, spread, Akp) | eta22 | psi22 | eta20 oscillatory |
|---|---:|---:|---:|
| 1, 25 deg, .02 | 0.267 | 0.140 | 1.873 |
| 2, 25 deg, .02 | 0.286 | 0.156 | 3.677 |
| 5, 25 deg, .02 | 1.445 | 0.822 | 14.649 |
| 1, 5 deg, .02 | 0.350 | 0.102 | 2.088 |
| 1, 25 deg, .12 | 4.431 | 2.668 | 4.460 |
| 5, 15 deg, .02 | 3.405 | 1.911 | 56.627 |
| 5, 15 deg, .12 | 2.374 | 1.309 | 29.059 |

For the high-steepness kh=1 case, narrowing the oscillatory eta20 band to
0.5*omega_p changes the error from 4.460% to 1.061%; both are retained. For
kh=5/spread15 the corresponding errors remain large: 54.653% at Akp=.02 and
26.941% at Akp=.12. Thus the deep-water discrepancy is not removed by that
particular temporal band choice. A raw four-phase average is not independently
certified as a pure bound difference-frequency field; input/model error,
unresolved/aliased harmonic content, adjustment waves and mean constraints
have not been separated. Do not diagnose this as a proven GL-kernel failure.

At kh=1/spread25/Akp=.02, farther probes give:

| Offset requested from center | eta22 (%) | psi22 (%) | eta20 oscillatory (%) |
|---|---:|---:|---:|
| -3 lambda in x | 0.274 | 0.289 | 3.869 |
| +3 lambda in x | 0.448 | 0.306 | 3.715 |
| -1 lambda in y | 20.979 | 210.965 | 6.876 |
| +1 lambda in y | 43.219 | 201.847 | 7.216 |
| -1.5 lambda in y | 82.581 | 71.967 | 10.152 |
| +1.5 lambda in y | 78.027 | 51.844 | 11.196 |

Observed-energy-weighted directional cancellation is about 1.84 for x offsets,
33.21 at |dy|=1 lambda, and 120.90 at |dy|=1.5 lambda. The far lateral signals
are much weaker than the center signal; `far_absolute_metrics.csv` preserves
their absolute amplitudes/errors alongside relative scores. These failures
must not be suppressed or clipped by a posteriori denominator regularization.

After inspecting the primary failure, a separately labelled convergence
diagnostic halves angular cell width from 7.5 to 3.75 degrees for both
|dy|=1 lambda points. No formula, GL rank, input band or target alignment changes.
For dy=-1/+1, eta22 errors become 16.458%/9.542%, and psi22 errors remain
104.847%/107.905%. The candidate-to-candidate change is still large. This does
not establish angular convergence or uniquely attribute the failure to the
allocation model; it demonstrates substantial discretization sensitivity in
a cancellation-sensitive configuration. The primary 7.5-degree scores remain
the main results, not replaced with the better refinement scores.

## Executed synthetic frequency-bandwidth outcomes

Measured relative frequency standard deviations are .08140, .16984 and .30749.
These cases use known directions and a constructed common frequency-dependent
modulation, so successful allocation is an implementation test under its own
assumption. The independent spectral MF12 order-two comparison tests the GL
kernel approximation; it is not evidence of real directional-evolution recovery.

| kh | Relative frequency width | eta22 GL8 L2 (%) | psi22 GL2+2 L2 (%) | eta20 GL16 L2 (%) |
|---|---:|---:|---:|---:|
| 1 | .08140 | 0.0123 | 0.0268 | 3.9929 |
| 1 | .16984 | 0.3183 | 1.0301 | 4.7227 |
| 1 | .30749 | 2.1420 | 6.1640 | 9.3516 |
| 5 | .08140 | 0.000675 | 0.00440 | 0.01323 |
| 5 | .16984 | 0.00362 | 0.01813 | 0.09109 |
| 5 | .30749 | 0.01768 | 0.12547 | 0.53970 |

These fixed-rank errors show a frequency-bandwidth dependence, particularly
at kh=1. They do not license fitting the scale or rank to an OW3D waveform.
Depth effects in the kernel benchmark and errors against raw OW3D phase
sectors are different questions and must not be conflated.

## Current interpretation

The joint-input method has useful evidence near the packet center and along
the main propagation direction. It is not yet robust at strongly cancelling
far lateral probes. Directional eta20 remains a diagnostic: some shallow/
intermediate cases are encouraging, but deep-water OW3D discrepancies remain
unresolved. Before widening claims, investigate directional allocation
conditioning and angular convergence, and separately establish the quality
of the eta20 observable. No new OW3D campaign, MF12 third-order calculation,
Stokes correction or reference-tuned repair was performed.

All new/modified MATLAB files passed Code Analyzer. The three Wolfram
difference-sector identities passed, and no released `src/`, `tests/` or
`symbolic/` file was changed. Additional refinement outputs live in
`kh1_s25_a002/y_m1_refinement` and `y_p1_refinement`; they are supplementary
diagnostics and are not substituted into the fixed primary summary.
