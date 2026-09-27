# Third harmonic and R second-order difference time-series pilot

Current execution update: the user explicitly rejected and requested removal
of the added interpolation/low-rank/cache route. Its active local source,
deployment/tests, experiment source bundles, deployed copies and saved operator
plan have been deleted. Retained reports are withdrawn-route provenance, not
evidence for the original GL algorithm's performance. Both earlier direct-triple
controllers remain stopped, with `superseded-by-fft.json` as historical stop
records; do not restart them. R4 completed results below are unchanged.
See SECOND_ORDER_IMPLEMENTATION.md for the actual second-order call path.

User requested third harmonic and second-order subharmonic, explicitly using
R for 2-minus. This pilot reuses the completed low-steepness GPU HOS record
(kpHs/2=.02, kph=1, depth 35.8423 m), the same observed first-harmonic sector,
initial complex directional spectrum, and original 10--70Tp scoring window.
Only the center probe is selected for the initial higher-order trial.

Remote run on 93:
`/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-higher-harmonics-20260927T145332Z`.
Controller PID at launch 2666977; MATLAB PID 2668029, CPU40, one thread.
Read status.json/progress.json/matlab.log; do not duplicate the running job.
Private run-local scientific adapter sources remain in ignored
`artifacts/hos_gpu/higher-harmonics-v1` and the hashed remote snapshot, not
vendored into public runtime code. Local evidence directory is
`artifacts/hos_gpu/higher-20260927T145332Z`.

## R second-order difference

Use the existing pure R4 Neumann pair kernel, transferred to arbitrary
directional temporal-bin parents. The algebra is evaluated directly by pair
summation with a stable polynomial in s^2/A. This is the same truncated R4,
not the exact denominator and not shared-scale GL. The physical assembly is
twice the existing spatial R half-field, as used by the original HOS initial
condition. R6 is a refinement diagnostic, not selected by matching HOS.

Direction widths 3.75 and 1.875 degrees use 4768 and 9536 active parents.
Exactly zero spatial difference K is excluded. Same-frequency/nonzero-K
pairs are retained, including their temporal-DC contribution in raw output.
The observed phase-average reference and raw predictions are saved with their
means unchanged. The oscillatory subharmonic comparison applies the identical
fixed 0 < abs(omega) < .5*omega_p mask to both; it does not validate the omitted
DC/mean-flow sector. Four-phase averaging also contains higher phase-zero
harmonics, so the raw reference is not pure second perturbation order.

First completed R results, fixed 10--70Tp window:

| Method | Direction width | Raw L2 | Oscillatory subharmonic L2 |
|---|---:|---:|---:|
| R4 | 3.75 deg | 13.9036% | 15.1596% |
| R4 | 1.875 deg | 11.8796% | 12.4312% |
| R6 | 1.875 deg | 11.9705% | 12.5567% |

HOS phase-average temporal mean is -1.64949e-4 m; R4 at 1.875 degrees gives
-6.06111e-5 m. Neither mean is fitted or removed in the raw records.

## GL third harmonic

The directional ordered-triple adapter transfers the public no-Stokes eta33
Euler forcing graph. Inner eta2/flat phi2 are generated from first-order
parents internally, never supplied from HOS. The frozen inner/outer tail
scales are 2nu_p-nu(2qp) and 3nu_p-nu(3qp). Both ranks 6 and 8 are tested.
The same 149 temporal parent bins satisfy kh>.5 and the temporal triple sum
Nyquist condition. There is no parent-energy cutoff or spatial-grid rounding.
Output triples beyond the original HOS spatial representable band are excluded
and counted; parents are retained.

Third-order direct cost is cubic in active directional parents. This initial
pilot therefore uses 30 and 15 degree direction bins (596/1192 parents), NOT
the earlier 1.875 degree resolution. Both direction and quadrature refinements
must be inspected; no fine-direction third-order validation is claimed.

HOS third-sector reference is
`(q0-q180+HilbertImag(q90)-HilbertImag(q270))/4`.
The HOS initial state had eta11+eta20+eta22 but no eta33. Consequently its
third-phase sector can contain startup-generated free components. Report raw
errors and the common parent-derived triple-frequency-band errors separately;
do not fit a phase/gain or remove transients to improve agreement.

## Transfer checks

Five Wolfram algebra checks passed (directional pair folds, unidirectional
direct limits, R polynomial equivalence). Across 12 checks each at two times,
two positions and three ranks/orders, temporal adapters matched existing
spatial routines with maximum absolute differences 6.7763e-20 m for R and
4.5297e-17 m for GL33. An independent three-parent directional test compared
temporal-bin FFT aggregation against direct exponential synthesis: relative
errors 3.1798e-15 (R4) and 6.4233e-15 (GL33). These certify implementation
transfer on fixtures, not physical accuracy of the HOS reconstruction.

`metrics.csv` and intermediate plots update as components finish. R_partial.mat
and higher_harmonics.mat stay remote. report.json marks full pilot completion
and reports the direction/rank differences. Preserve large errors and all
partial results; a completed computation need not be an accurate approximation.

## First third-harmonic result and reference diagnostic

At 30 degree directional resolution, GL33 rank 6 gives raw/common-triple-band
L2 errors 92.1935%/80.2993%; rank 8 gives 92.1448%/80.1669%. The 15 degree
cases continue in the original controller. These are coarse-direction trial
results, not validated third-harmonic reconstruction.

`audit_third_reference.m` synthesizes a PURE LINEAR field from the original
initial complex spectrum, forms four phases, and applies exactly the same
finite-record Hilbert third-sector extraction. Its true third harmonic is zero.
Nevertheless, the extracted spurious sector has interior RMS 5.5805e-4 m,
compared with 8.0677e-4 m for the HOS extracted sector (ratio 0.6917).
Within the common triple-frequency band these RMS values are 2.7816e-4 and
5.2438e-4 m (ratio 0.5305). This establishes that finite-record separation
leakage is material at the third-harmonic scale. The evolved HOS first harmonic
differs from the linear control; these ratios are not exact contamination
fractions. The control is NOT subtracted or fitted to the HOS reference.
The original raw references, errors and plots remain unchanged.

## User-requested zero-mean and endpoint diagnostic

The user explicitly requested setting the means to zero. Separate diagnostics
are under remote `demean-check-v1` and `fixed-taper-check-v1`, with compact
metrics/report/plots mirrored in the local evidence directory. These do not
overwrite the original scientific records or alter HOS initial conditions.
They use the completed 30-degree, rank-8 GL prediction.

Subtracting each HOS phase's full-record temporal mean makes the extracted
third-sector mean 2.19e-18 m. The original 10--70Tp raw error is 92.1448%;
full-record mean removal gives 92.9876%; independent scoring-window mean
removal gives 92.1448%. Selected nonzero GL input coefficients change by only
1.72e-16 relative. Common-triple-band error stays 80.1669%. The endpoint jump
is unchanged. The full-record third-sector mean (2.75e-4 m) is substantially
different from its interior scoring-window mean (1.78e-6 m), so subtracting
the full-record mean does not eliminate the interior time-varying leakage.

A separate fixed cosine taper covers the first/last 10% of the record
(110.08 s each) and is identically one throughout the original scoring window.
Apply the same taper, demeaning and four-phase extraction to HOS and to the
analytic third-only GL phases. The error is then 69.8615%. In a pure-linear
control, spurious third-sector RMS drops from 5.58e-4 to 6.50e-5 m (11.65%
remaining). Processing changes the GL prediction within the scoring window
by 1.00e-7 relative. This supports material endpoint leakage but does not
establish accurate/converged GL third-order reconstruction. It is a processing
diagnostic: GL was not recomputed from a newly windowed first-harmonic input.

Initial absence of eta33 does not prevent HOS generating third-order content
through nonlinear evolution. It is a possible source of transient/free
components to diagnose, not grounds for rejecting the comparison. This
low-family steepness is kpHs/2=.02, not monochromatic akp. A larger-steepness
comparison requires a new nonlinear HOS integration, not rescaling its evolved
records. No new-steepness simulation has been launched by these diagnostics.

## User-requested first-harmonic input through 4omega_p

The user rejected narrow edge-padding tests and explicitly requested the full
first harmonic, or retaining all components through 4omega_p. The latter is
now running under `first-band-4wp-v1` in the same remote parent directory.
Controller launch PID 2751775; MATLAB PID 2751777, CPU41. At 15:32 UTC it had
entered `all_nonzero_to_4wp` with 1280 active directional parents: all 320
nonzero temporal frequency bins through 4omega_p and four active direction
bins at the unchanged 30-degree resolution. The old initial-band baseline
had 149 frequencies and 596 directional parents.

The only input frequency exclusions in this new variant are DC and omega
above 4omega_p. There is no original-spectrum min/max cutoff or kh>.5 input
deletion. A run-local copy of the same GL triple kernel relaxes its prior
kh>.5 assertion to q>0; its numerical formula is unchanged. Consequently
low-kh quadrature convergence is unverified, and the directional allocation
outside the initial physical band uses the finite-record prior's FFT tails.
Both limitations must be reported, not hidden by silently dropping bins.
The HOS representable OUTPUT wave-number projection remains unchanged.

Direction resolution, GL rank 8, raw HOS reference, fixed common output band,
fixed taper diagnostic and scoring window are identical across the new and
old input-band cases. Original in-band directional coefficients are checked
unchanged. Report input projection error, third-prediction change, conditioning
and HOS errors independently; never select the input band by its HOS score.

The preceding `first-band-check-v1` narrow extension pilot was intentionally
stopped on the user's steering; its `superseded.json` identifies this reason.
Its controller may show a nonzero MATLAB exit from the intentional stop;
do not call that a scientific failure. Partial data are retained.

This adapter explicitly sums ordered triples and has cubic parent-count cost.
That is a property of this reference implementation, not a statement that a
fixed-grid/fixed-rank FFT Green--Laplace implementation needs cubic work when
more spectrum entries are nonzero. No fast time-series GL performance claim
has been established by these diagnostic runs.
