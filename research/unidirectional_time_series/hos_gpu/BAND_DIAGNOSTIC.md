# Fixed-support precision diagnosis

The user requested testing the explanation for the large-grid FP32 RHS error.
The physical input, h, Lx, Ly, gravity and predeclared tolerances remain fixed.
No hardware settings are changed. The original default unprojected path is
retained. This is an instantaneous order-two diagnosis, not an M=5 validation.

## Predeclared experiment

Input modes are (0,0), +/-(1,2), +/-(3,-1) in 2D, with the same complex
eta/psi coefficients as rhs_check.cu. The fixed enclosing input rectangle is
|mx|<=3, |my|<=2. Every quadratic sum lies in |mx|<=6, |my|<=4.
These rectangles follow from input support and polynomial degree; they are
not selected from reference errors and retain every physical test mode.

Four variants run with identical physical samples:

0. Original path, no projection.
1. Input projection: project eta and psi to the input rectangle before
   constructing G0psi and derivatives. Keep original uploaded states intact.
2. Product projection only: project eta*G0psi, eta*psi_x, eta*psi_y to the
   quadratic rectangle immediately before their spectral derivatives.
3. Both projections.

psi_t's pointwise quadratic terms are not post-filtered in any variant.
The two switches default to false; production is not silently changed.

## Measurements

The independent reference uses signed sparse Fourier convolution, collecting
coefficients once and evaluating with separable double-precision phase tables.
It uses no FFT or GPU operator. Full-grid reference fields are computed for
W=G0psi, A=G0(eta G0psi), Bx=dx(eta dx psi), By=dy(eta dy psi), and the final
eta_t=W-A-Bx-By and psi_t. This separates intermediate errors from the final
subtraction error. Reference conjugate symmetry is checked explicitly.

Before each optional projection, copy the input/product spectrum for diagnosis.
Measure maximum normalised coefficient and Parseval RMS outside the appropriate
rectangle. Also calculate those residual measures after weighting by the
analytic derivative/G0 multiplier (double-precision diagnostic calculation,
not a separate captured GPU transform). R2C interior x columns have weight 2;
the zero and Nyquist x columns have weight 1. The physical-grid errors are
measured against the actual GPU outputs.

FP32 tolerance remains 2e-5; FP64 remains 2e-11, maximum absolute RHS error.
Every variant reports its own pass/fail. Executable exit status checks variant
3, allowing the expected failing baseline to remain in the same report.
Repeated callback-free evaluations must reproduce callback-enabled outputs
exactly, checking that the diagnosis does not change arithmetic or input state.

Run `python run_p40_checks.py --bands` for a fresh hash-verified P40 snapshot.
It also runs earlier FFT/RHS regression tests. Artifacts are stored outside
tracked source. Large-grid runs compare all four variants at h=.15 m in both
precisions, and at h=1.3,20 m in FP32. CPU reference evaluation is optimised
with a Release build; GPU speed is not benchmarked.

## Observed results — 27 September 2026

Hash-verified snapshot `/root/hos-rhs-check-20260927T004353Z` on jfm92/P40.
Raw log/manifest: `artifacts/hos_gpu/20260927T004353Z/`. Compact values and
source hashes: `band_diagnostic_results.json`. Eight CTest tests passed,
including original unprojected small-grid regressions and projection tests.
The original 256x128 unprojected checks also passed at all three depths.

2048x1024, FP32, h=0.15 m:

| Variant | eta_t max abs error | psi_t max abs error | Original 2e-5 threshold |
|---|---:|---:|---|
| No projection | 4.354397584644884e-4 | 3.276512588024083e-6 | FAIL |
| Input only | 5.231239916636543e-6 | 6.250105397676009e-7 | PASS |
| Product-before-derivative only | 3.0951333857554e-5 | 3.276512588024083e-6 | FAIL |
| Both | 8.408705445828168e-8 | 6.250105397676009e-7 | PASS |

At h=1.3 and 20 m, the both-projected FP32 eta_t errors were respectively
1.7590451248383765e-7 and 1.7244998734566863e-7; psi_t errors were
6.520496145245858e-7 and 6.204991702674079e-7. All pass unchanged thresholds.
At h=0.15 m FP64 both-projected errors were 2.1163626406917047e-16 and
2.0816681711721685e-15. All four FP64 variants passed.

## Mechanism supported by the measurements

Raw psi has out-of-band normalised spectral peak 3.83994e-9 and RMS
1.53662e-8. Applying the G0 multiplier to those measured residuals predicts
an out-of-band RMS 6.02698e-6. The actual physical W error is 3.09523e-5
without input projection and 4.64143e-8 with it.

For A=G0(eta G0psi), the unprojected intermediate physical error reaches
1.47247e-3. With input projection only it is 2.29747e-6; with both it is
1.02389e-8. Bx and By errors likewise fall from 1.07169e-3 / 5.09100e-4
to 3.06179e-8 / 3.28715e-8. In this example the large intermediate errors
partly cancel in eta_t; the evidence does not show that cancellation itself
is the principal source. Product-only projection removes the large quadratic
term errors, but leaves W's input-noise error, explaining its remaining fail.

Thus amplified out-of-band roundoff is the dominant cause in this sparse
test. Input projection alone leaves newly generated product/FFT roundoff;
product projection alone leaves input-derivative noise. Both address the
observed chain without deleting any physical mode in the declared test.

This does not certify a universal narrow cutoff, full broadband data, fifth
order, or long-time FP32 accuracy. The next implementation must obtain bands
from the actual discretisation/input contract and HOS-Ocean's order-specific
rules, rather than reusing this test's 3x2 rectangle. Defaults still preserve
the original unprojected baseline; no existing HOS run was replaced.
