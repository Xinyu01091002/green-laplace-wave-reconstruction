# Changelog

## 0.3.0 - 2026-09-23

- Replaced the third-order API and paper executors with the pure GL graph,
  removing Stokes-diagonal repairs in the nested second-order states and
  the third-order elevation and surface-potential responses. This changes
  numerical results from earlier releases; an exact diagonal is not imposed.
- Added the independent no-Stokes fourth-order component interface and the
  third-order Two-Scale/Shared-Scale research graphs. The unified spectral
  API continues to accept orders one through three only.
- Added the fourth-order GL--WIT paper package: historical C++ sources,
  frozen inputs, WIT reference output, GL6/8/10 fields, timing measurements,
  plotting commands and build instructions. Historical timings remain
  distinct from newly executed migration checks.
- Added source-integrity, input-reconstruction, C++/MATLAB parity and
  fourth-order paper reproduction checks, including a dedicated CI workflow.
- Reorganized the public guide around reconstruction, paper reproduction
  and symbolic sources; documented dependencies for each task.

## 0.2.0 - 2026-08-27

- Added portable nonzero difference-frequency eta20 diagnostics.
- Added frozen GL12/GL16 shared-scale rank configurations and figures.
- Added exact Wolfram and JSON interfaces for the Neumann R2/R4/R6 sequence.
- Added an independent eta11-only ordered-pair R-series evaluator; it is
  explicitly validation code rather than a fixed-FFT production route.
- Updated eta22 and surface-potential figures, metrics, centerlines and
  mixed-rank cost--accuracy data.
- Removed superseded GL8-only eta20 reference assets; they remain recoverable
  from version 0.1 history.

## 0.1.0 - 2026-08-24

- Added the MF12-compatible two-stage MATLAB interface.
- Added paired physical `eta` and true surface `psi` outputs through the
  released positive pure-sum order-three route.
- Added prescribed-rank `eta22`, dual-branch `psi22`, and bounded
  `eta33/psi33` Green--Laplace graphs.
- Added an optional same-input MF12 order-two comparison.
- Added exact Wolfram residual and freeze sources with committed JSON
  interfaces.
- Added portable manuscript figure entry points and frozen reference assets.
- Added MATLAB R2022b release tests and GitHub Actions CI.
