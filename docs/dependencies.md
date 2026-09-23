# Dependencies

## Required

### MATLAB

The release is tested with MATLAB R2022b using double precision. Only base
MATLAB functionality is used: `fft2`, `ifft2`, JSON import, tables, argument
validation and plotting. No toolbox is required by the public executors or
the release tests.

Inputs must lie on the requested FFT grid. The public frontend fails if a
component is rounded onto a different wavenumber or if the nonlinear output
support reaches a Nyquist boundary.

## Optional

### Fourth-order C++ reproduction

The [fourth-order GL--WIT package](../paper/order4/README.md) requires Linux
or WSL, a C++17 compiler, OpenMP, FFTW3 with its threads library, GNU patch,
and Python 3 for process orchestration. MATLAB R2022b provides independent
numerical validation. FFTW is an external build dependency; its installation
can be selected with `FFTW_PREFIX`.

These dependencies are needed to build and rerun the C++ calculations.
Replotting the preserved fourth-order fields and timings, or recomputing the
GL fields against the preserved WIT reference, needs only base MATLAB.

### MF12 MATLAB implementation

MF12 is needed only for comparison examples. It is kept external so that its
history, license and citation remain intact. `setup_green_laplace` accepts
either the repository root or the directory containing
`mf12_spectral_coefficients.m`.

The comparison helper uses MF12 spectral coefficients and reconstructs only
the positive pure-sum order-two `eta22` and surface `psi22` sectors. It does
not compare a GL sector against a differently defined MF12 total field.

### Wolfram Mathematica

WolframScript is needed only to regenerate files in `symbolic/generated/`.
The committed JSON interfaces allow ordinary MATLAB users to run the
software without Mathematica.

Run the frozen gates from the repository root:

```powershell
wolframscript -file symbolic/wolfram/order2_finite_depth_directional_eta22_pure_gl8_freeze.wl
wolframscript -file symbolic/wolfram/order2_finite_depth_directional_psi22_dual_branch_gl_freeze.wl
wolframscript -file symbolic/wolfram/order3_finite_depth_directional_crossing_stokes_gate.wl
wolframscript -file symbolic/wolfram/order3_finite_depth_directional_nested_green_laplace_freeze.wl
```

Generated run reports are written under the ignored `artifacts/` directory.

The retained third-order symbolic sources include historical Stokes traces
and correction identities. These are mathematical reference records; the
version 0.3.0 MATLAB execution path does not apply those repairs. See the
[migration record](no_stokes_migration.md) before regenerating historical
interfaces.
