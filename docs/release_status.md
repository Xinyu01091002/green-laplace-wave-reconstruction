# Version 0.2.0 release status

Date: 2026-08-27

## Version 0.2 additions

- Public eta20 GL6/GL12/GL16 diagnostic code matches the source
  implementation exactly on a common deterministic fixture.
- Public R2/R4/R6 ordered-pair code matches the frozen source implementation
  exactly on the same fixture.
- R2/R4/R6 Wolfram inverse, residual, endpoint and compiler gates pass.
- GL12 and GL16 Gauss--Laguerre moment and positivity gates pass.
- Portable eta20 figures reproduce the recorded `k_p d=1,2` metrics.
- The committed surface-potential Figure 5 uses the uncorrected pure-GL
  resolvent consistently: GL6 has `Q=0.0346698`, raw L2 `6.75935%`, on the
  8-wavelength matched domain.
- Fixed-FFT R2/R4/R6 execution is intentionally not claimed as validated;
  the public R evaluator is an ordered-pair accuracy diagnostic.

## Executed gates

- MATLAB R2022b release suite: six checks passed.
- MATLAB Code Analyzer: zero messages across the committed `.m` files.
- Public dependency audit: base MATLAB only.
- Exact Wolfram order-two eta22 freeze: passed.
- Exact Wolfram order-two surface-psi22 freeze: passed.
- Exact Wolfram order-three crossing gate: passed.
- Exact Wolfram nested eta33/psi33 freeze: passed.
- Self-contained paper figure command: completed.
- Optional MF12 order-two example: completed without rescaling or alignment.

The deterministic two-parent MF12 smoke case gave raw relative L2 values of
`4.39934e-4` for eta22 and `6.43164e-5` for surface psi22. These values are a
small API/parity smoke test, not a broadband accuracy claim.

## Reproduction boundary

Figures 1--2 are self-contained. The full figure mode regenerates the
MF12-dependent eta22 and eta20 comparisons. Figure 5 and the mixed-rank
diagnostic require the separately archived matched MF12 fields listed in
`paper/data_manifest.json`. Difference-frequency eta20 remains outside the
released total-field API and is labelled diagnostic throughout.
