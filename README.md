# Green--Laplace wave reconstruction

MATLAB implementation and reproducibility material for Green--Laplace
reconstruction of finite-depth directional bound waves. The public interface
follows the two-stage MF12 pattern: prepare a spectral reconstruction once,
then evaluate the surface elevation and the true free-surface potential on a
regular FFT grid.

The current third-order runtime is the pure Green--Laplace graph with all
Stokes-diagonal repairs removed, including repairs in its nested second-order
states. The paper entry points use that same implementation. This changes
third-order results from earlier releases; an exact monochromatic diagonal
is not imposed. See [the migration record](docs/no_stokes_migration.md).

```matlab
setup_green_laplace

gl = gl_spectral_coefficients( ...
    order,g,h,a,b,kx,ky,Ux,Uy,options);

[eta,psi,X,Y,components,audit] = ...
    gl_spectral_surface(gl,Lx,Ly,Nx,Ny,t);
```

The first four outputs of `gl_spectral_surface` match the ordering of
`mf12_spectral_surface`. `psi` always denotes
`phi(x,z=eta,t)`. Flat bulk-potential traces are internal states and are not
silently returned as surface potential.

## Released outputs

| Order | Components | Sector | Status |
|---:|---|---|---|
| 1 | `eta11`, `psi11` | linear | released |
| 2 | `eta22` | positive pure sum | prescribed GL rank |
| 2 | `psi22` | positive pure sum | dual-branch GL2+2 |
| 3 | `eta33`, `psi33` | positive pure sum | bounded GL4 reconstruction |

Additional preserved GL research implementations are available separately:

- `gl_pure_sum_order4`: experimental no-Stokes `eta44` and surface `psi44`,
  with dimensionless analytic inputs and outputs documented in its help;
- [`research/two_scale`](research/two_scale/README.md): independent
  Two-Scale/Shared-Scale third-order elevation and surface-potential graphs,
  with direct ordered-triple parity checks.

These research interfaces do not extend the supported order range of
`gl_spectral_coefficients` or silently add components to the total field.

The paper's [fourth-order GL--WIT comparison package](paper/order4/README.md)
includes the C++ implementations, frozen field and timing inputs, original
reference outputs, GL6/8/10 fields, measured timing rows, and standalone
reproduction commands. All GL execution paths in that package omit Stokes
corrections. Historical timings and newly rerun validation are distinguished.

Use `gl_supported_sectors` for the machine-readable applicability table.
Unsupported difference-frequency, strict-zero, free-wave, resonant and
near-resonant sectors fail explicitly; they are never filled with zeros.

Version 0.2 also includes a separate diagnostic package for nonzero
difference-frequency `eta20`:

- shared-scale GL6/GL12/GL16 rank diagnostics;
- exact frozen Neumann R2/R4/R6 formulas;
- independent ordered-pair field validation for the R sequence.

These diagnostics are not silently added to the total field returned by
`gl_spectral_surface`. R2/R4/R6 are not Green--Laplace ranks, and their
fixed-FFT production implementation is not claimed as validated.

## Dependencies

| Dependency | Purpose | Required |
|---|---|---|
| MATLAB R2022b | GL execution, tests and figures | yes |
| Base MATLAB | FFT, JSON, tables and plotting | yes |
| External MF12 MATLAB repository | GL--MF12 examples | optional |
| Wolfram Mathematica / WolframScript | regenerate exact frozen interfaces | optional |

The release tests use base MATLAB only. No Python package or MATLAB toolbox is
required. The current local release gate was executed with MATLAB R2022b;
newer releases are not claimed until the same checks run there.

## Quick start

Run the self-contained example:

```matlab
setup_green_laplace
run("examples/run_minimal_example.m")
```

Run the release checks:

```matlab
setup_green_laplace
addpath("tests")
run_release_tests
```

## MF12 comparison

MF12 is an external optional dependency and is not copied into this
repository. Configure it once:

```matlab
setup_green_laplace("MF12Root","C:/path/to/mf12-repository")
```

or set `MF12_ROOT` before running:

```matlab
run("examples/run_mf12_order2_comparison.m")
```

The comparison isolates positive pure-sum `eta22` and surface `psi22` from
the same parent components. It reports raw relative L2, raw relative Linf,
norm ratios and signed extrema. No gain, offset, phase, time shift, spatial
shift or fitted rescaling is applied.

## Reproducibility

- `symbolic/residuals/` contains the exact Euler residual sources.
- `symbolic/wolfram/` freezes the Green--Laplace interfaces.
- `symbolic/generated/` contains the machine-readable interfaces consumed by
  MATLAB.
- `paper/` contains the portable figure entry points and frozen manuscript
  reference assets.
- `diagnostics/eta20/` contains the explicitly non-production GL-rank and
  Neumann R-series difference-frequency studies.
- `symbolic/SOURCE_MANIFEST.yml` records the extraction source and hashes.

Regenerate the self-contained paper figure set with:

```matlab
setup_green_laplace
addpath("paper")
run_all_figures("self-contained")
```

The `"full"` mode adds the MF12-dependent second-order field calculations.
Data-heavy manuscript comparisons remain identified separately in
`paper/README.md`; the repository does not pretend that a locked raster is an
end-to-end numerical reproduction.

## Scope

The released nonlinear graphs accept first-order surface elevation only.
They assume strict-forward analytic support and alias-safe FFT grids. The
order-two surface-potential graph requires parent `kh >= 0.3`; order three
requires every parent `kh > 0.5`. The MF12-compatible `Ux`, `Uy` slots are
present, but the released total-field API requires both to be zero. Evidence outside these
declared domains is not a release claim.

This software is a fixed-order computational reformulation of regular bound
wave interaction kernels. It is not a solver for strict-zero modes, mean
flow, resonant primary-harmonic corrections or coupled nonlinear evolution.

## License and citation

Source code is released under the MIT License. Cite the repository using
`CITATION.cff`; cite the Green--Laplace paper and MF12 theory separately when
using the associated methods or comparison implementation.
