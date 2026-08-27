# Applicability and limitations

- Strict-forward analytic input support is required: every parent has
  `kx > 0`.
- The MF12-compatible `Ux` and `Uy` argument slots are retained, but version
  0.1 requires `Ux=Uy=0`; a uniform-current phase convention has not been
  certified for the released nonlinear graphs.
- `psi22` uses the declared parent domain `kh >= 0.3`.
- `eta33` and `psi33` require every parent `kh > 0.5`.
- FFT grids must be alias-safe for the requested nonlinear order.
- The order-two elevation rank is prescribed by the user; MF12 data do not
  select or retune it.
- Field metrics are computed in the same physical units without alignment or
  fitted rescaling.
- Difference-frequency, strict-zero, free-wave, resonant and near-resonant
  coupled-evolution sectors are outside the released total-field API.
- The eta20 R2/R4/R6 package uses independent ordered-pair reconstruction to
  validate approximation accuracy. It does not claim validated fixed-FFT
  execution, timing, or production readiness.
- GL12/GL16 eta20 results are post-baseline rank diagnostics on fixed inputs,
  not new default settings or uniform weak-detuning guarantees.
- Fourth- and fifth-order symbolic constructions are not advertised as
  released numerical surface APIs until their same-variable validation and
  dependency graphs are complete.
