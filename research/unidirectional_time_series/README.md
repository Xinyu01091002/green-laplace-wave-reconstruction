# MATLAB time-series reference executors

This directory contains the MATLAB modal FFT--GL executors used to generate
the single-direction Python validation fixtures:

- `gl_directional_time_modal_grid.m` produces `eta22(t)`;
- `gl_directional_time_modal_grid_eta33.m` produces `eta33(t)`;
- `gl_directional_time_modal_grid_eta20.m` produces nonzero-wavenumber
  difference-frequency `eta20(t)`.

They retain filtered fields, pointwise products and FFT operators and do not
enumerate wave pairs or triples. They do not use an auxiliary spatial domain.
The Python fixture generator is `python/tools/generate_matlab_fixtures.m`.
