# Nonzero eta20 and positive-sum eta33 time-series trials

The user authorized commit/push of the current research branch and continuation
with eta20 and eta33, explicitly excluding MF12 third-order comparisons.
The existing branch through `2ff32d4` was pushed before these trials. Public
main and the published spatial executors are unchanged.

## Fixed input and comparison definitions

Use the saved boundary-corrected Alpha=1, kh=1, Akp=0.02/0.12 raw four-phase
records verified in `ALPHA1_STEEPNESS.md`. Construct the same first phase sector
as the eta22 trial. For this eta33 experiment, the common parent domain is
kh>0.5 (the released unified order-three domain) with three times the parent
frequency below the native temporal Nyquist. This explicit shared support is
used for *both* eta20 and eta33 and for every method. No energy cap, spatial
wave-number rounding, interpolation, reference-fitting or gain is used.
It differs from the preceding eta22-only support; retained energy is recorded.

The same input envelope maximum defines the displayed +/-2Tp main-group window.
Each output is synthesized on the original temporal grid; the full time record
is saved. Comparisons remain in physical metres with no fitted alignment.

For eta20, the observed record is the average of the four phases. Four phases
also admit higher harmonic content. Following the existing VWA time-series
script's stated subharmonic band, apply `0 < abs(omega) < 0.5*omega_p` identically
to reference and every prediction. Excluding the zero-frequency FFT bin is a
declared sector projection, not a fitted offset. The original observed mean
and unfiltered candidate/reference fields are saved. This trial does not
predict the strict-zero mean, tank-volume constraint or full mean-flow sector.
The chosen band is not a claim to contain all broadband difference interactions.

For eta33, use the existing temporal Hilbert-four-phase combination
`(q0-q180+Ht(q90)-Ht(q270))/4`. No reference-frequency filter is selected to
improve scores. Both eta20 and eta33 references remain phase/harmonic sectors,
not pure perturbation orders, particularly at Akp=0.12.

## GL transfer and validation

`gl_eta20_time_pairs` transfers the frozen shared-scale GL6/12/16 diagnostic in
`diagnostics/eta20/eta20_green_laplace_shared.m`. It uses the same input-derived
scale sqrt(2*energy-weighted variance of kh) and true wave-number/frequency
differences. Strictly diagonal pairs are excluded. The existing published
`figure_eta20_gl_rank.m` explicitly assembles `2*candidate_half`; the temporal
adapter preserves this physical factor, including u=A/2 for the positive
Fourier coefficient. It is not calibrated from MF12 or OW3D. Shared-scale GL
eta20 remains an accuracy diagnostic, not production or a uniform weak-detuning
guarantee. Neumann R2/R4/R6 are not silently substituted or called GL ranks.

`gl_eta33_time_triples` transfers the pair fields and cubic forcing of
`src/internal/gl_no_stokes_eta33.m`, including its prescribed inner/outer scales
2*nu_p-nu(2qp) and 3*nu_p-nu(3qp). It retains both sinh and cosh Laplace sums,
rather than replacing cosh by Omega times a resolvent approximation. Inner
eta2 and flat phi2 are generated internally from eta1. Triple outputs carry
physical A_a*A_b*A_c/h^2. No Stokes-diagonal or angular repair is introduced.
Ranks 4/6/8 are fixed in advance; eta33 involves no MF12 execution.

Uniform temporal-bin triples are accumulated at their integer sum frequencies,
then synthesized exactly with the negative temporal exponential. No spatial
FFT grid is invented. This remains an O(M^3) reference algorithm, not a new
fast GL time-series algorithm or timing claim.

Wolfram `verify_eta20_eta33_transfer.wl` checks five transfer identities, all
passing. MATLAB checks nonzero complex amplitudes, two times, three positions
and all requested ranks against the existing spatial implementations:

- eta20 physical assembly: max absolute difference 1.0842e-19 m;
- eta33: max absolute difference 4.8773e-17 m;
- cubic uniform-bin aggregation versus direct exponential synthesis:
  relative L2 6.4894e-15.

These certify algebraic transfer/implementation parity on the fixture, not
independent higher-order physical accuracy over arbitrary time records.

## Other methods

- eta20: spectral MF12 **order two only**, negative-pair coefficients with
  strict-zero terms excluded; Walker Eq. (14), with no fitted scale. Walker's
  formula is `-kp/(2*sinh(2*kp*h))*abs(z)^2` before the common band projection.
- eta33: the VWA temporal product from `paperplot_VWA_time_series.m`,
  `Re[z * L_(k B33) z * L_k z]`, and Walker Appendix A's peak-coefficient
  `kp^2*B33(kp*h)*Re[z^3]`. The alternative VWA expression in some later
  comparison scripts is not selected according to error. No MF12 order-three
  coefficients are requested, computed or plotted.
- An independent VWA eta20 temporal approximation has not been certified here;
  the old VWA plotting script itself used Dalzell for eta20. No new VWA eta20
  formula is invented for this plot.

## Reproduction

```powershell
wolframscript -file research/unidirectional_time_series/verify_eta20_eta33_transfer.wl
matlab -batch "restoredefaultpath; addpath('research/unidirectional_time_series'); test_eta20_eta33_time; run_eta20_eta33_pilot('C:/Research/spectral domain implementation of wave interaction theory');"
```

Outputs: `results/unidirectional_time_series/eta20_eta33_boundary_alpha1_akp002/`
and `..._akp012/`. Each contains raw fields, method/sector metadata, full and
main-group metrics, and separate eta20/eta33 PNG and PDF figures. Source MAT
hashes and local commit provenance accompany the saved results.

## Executed results

Both cases use 98 common parents, minimum kh=0.510634. Retained input energy
is 99.9998996% at Akp=0.02 and 99.9999757% at Akp=0.12. Input projection L2
is 0.10025% and 0.04939%. This is explicitly a different common support from
the earlier eta22 trial, required by the declared cubic domain/Nyquist gate.
The eta20 projection retains eight positive and eight negative temporal bins,
with cutoff 0.2282803 rad/s. Observed DC excluded from this sector is
-0.000657712 m / -0.023957419 m; the strict mean is not reconstructed.

Raw relative L2 against the corresponding OW3D projected phase sector in
the input-defined main-group window:

| Sector | Method | Akp=0.02 (%) | Akp=0.12 (%) |
|---|---|---:|---:|
| eta20 | GL6 diagnostic | 7.9154 | 9.6120 |
| eta20 | GL12 diagnostic | 2.6232 | 4.7325 |
| eta20 | GL16 diagnostic | 2.2463 | 4.3139 |
| eta20 | spectral MF12, order 2 | 2.1407 | 4.1349 |
| eta20 | Walker Eq.14, unscaled | 82.8395 | 83.5765 |
| eta33 | GL4 | 14.7162 | 21.6283 |
| eta33 | GL6 | 5.3011 | 13.6108 |
| eta33 | GL8 | 2.1583 | 11.4886 |
| eta33 | VWA temporal product | 18.0446 | 7.6450 |
| eta33 | Walker | 9.4893 | 7.4643 |

The eta20 diagnostic approaches the independent quadratic comparison on this
declared band. It does not certify the excluded mean, other difference bands,
or other depths. Eta33 shows a substantial prescribed-rank dependence; a low
rank must not be reported as a converged GL result. The remaining large-
steepness discrepancy has not been separated into GL quadrature error and
the difference between cubic reconstruction and OW3D harmonic content.
VWA's smaller large-steepness error than its small-steepness error is not
proof of improved asymptotic accuracy; approximation errors can cancel
higher-order input/output contributions. No compensation is fitted here.

## User-requested eta20 band extension to 3 omega_p

`reproject_eta20_band(3)` uses the saved, unfiltered eta20 candidates and OW3D
four-phase-average record. It applies the identical mask
`0 < abs(omega) < 3*omega_p` to all methods, leaving the first-order input,
GL ranks/scales, interaction kernels and main-group window unchanged.
There is no new physical calculation or three-order MF12 execution.
The strict-zero mode remains excluded. Original 0.5*omega_p files remain intact;
new output directories append `_eta20_band3` and include the original-versus-new
main-group error table. This changes the output comparison bandwidth only.
The enlarged OW3D phase-sector band should not be called a strict perturbation-
order separation merely because its frequency limit has been increased.

The 3*omega_p experiment passed finite-value and Code Analyzer checks. It
retains 96 nonzero FFT bins rather than 16. Main-group relative L2 (%):

| Method | Akp=0.02, band 0.5 -> 3 | Akp=0.12, band 0.5 -> 3 |
|---|---:|---:|
| GL6 diagnostic | 7.9154 -> 7.8149 | 9.6120 -> 9.6965 |
| GL12 diagnostic | 2.6232 -> 2.5803 | 4.7325 -> 4.9904 |
| GL16 diagnostic | 2.2463 -> 2.2047 | 4.3139 -> 4.5906 |
| Spectral MF12, order 2 | 2.1407 -> 2.0988 | 4.1349 -> 4.4190 |
| Walker Eq.14 | 82.8395 -> 82.8925 | 83.5765 -> 83.6290 |

The wider output band does not substantially change these relative-error
conclusions. This observation does not establish independence from all other
filter choices or remove the phase-sector versus perturbation-order distinction.
