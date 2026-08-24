# Version 0.1.0 release status

Date: 2026-08-24

## Executed gates

- MATLAB R2022b release suite: five checks passed.
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

Figures 1--2 are self-contained. Figures 3--4 have portable MF12-dependent
generators and retain the current manuscript reference tables/assets. The
full broadband campaign was not recomputed as part of repository extraction.
Figure 5 requires the separately archived matched MF12 field data. The
difference-frequency Figures 6--7 are reference-only in version 0.1 because
that sector is not part of the released unified API.
