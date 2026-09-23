# Third-order Two-Scale and Shared-Scale GL

These independently frozen research graphs reconstruct positive pure-sum
eta33 and surface psi33 from the analytic eta11 spectrum. They have no
Stokes-diagonal repair. Internally generated lower-order states are included
in their execution counts. Wavenumbers and spectra use the dimensionless
convention of the original freeze; inspect the function help and the compact
fixtures before constructing an input.

Run the independent ordered-triple checks from the repository root:

```matlab
setup_green_laplace
addpath('research/two_scale')
run_ordered_triple_parity
run_psi33_ordered_triple_parity
```

The elevation modes are `two_scale` and `shared_scale`. The surface-potential
graph additionally supports `inner_shared_outer_two_scale`. Both enforce
their declared parent support, strict forward cone and alias-safe grid.
The `surface_correction` variable in the surface implementation is the Taylor
evaluation from the flat trace to the free surface; it is not a diagonal repair.

Exact interfaces live in `symbolic/generated`; their Wolfram generators and
original Euler residual sources are included in `symbolic`. Validation does
not select or retune coefficients. These interfaces remain separate from the
main spectral API and do not inherit historical corrected-version accuracy
or performance claims.
