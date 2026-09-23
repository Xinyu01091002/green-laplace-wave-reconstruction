# Green--Laplace wave reconstruction

MATLAB implementations and paper reproduction material for Green--Laplace
reconstruction of finite-depth directional bound waves from first-order
surface elevation. The repository also includes the C++ implementations,
inputs and results used for the paper's fourth-order comparison with Wave
Interaction Theory (WIT).

| I want to... | Start here |
|---|---|
| Run a GL reconstruction in MATLAB | [Quick start](#quick-start), then the [minimal example](examples/run_minimal_example.m) |
| Reproduce the paper's results | [Paper guide](paper/README.md); for fourth-order fields and timings, use the [GL--WIT package](paper/order4/README.md) |
| Inspect the derivation and frozen formulas | [Euler residuals](symbolic/residuals), [Wolfram generators](symbolic/wolfram), and [exported interfaces](symbolic/generated) |

## Quick start

The example requires **MATLAB R2022b and base MATLAB only**. No additional
toolbox, Mathematica, MF12 installation or C++ build is needed. R2022b is the
tested release; newer MATLAB releases have not been separately certified here.

Clone the repository, or download and extract its ZIP from GitHub:

```bash
git clone https://github.com/Xinyu01091002/green-laplace-wave-reconstruction.git
cd green-laplace-wave-reconstruction
```

Set MATLAB's current folder to the repository root, then run:

```matlab
setup_green_laplace
run('examples/run_minimal_example.m')
```

The example reconstructs the released first- to third-order positive pure-sum
components, displays the third-order elevation and surface potential, and
saves fields and a figure under `results/minimal_example/`.

To check the installation:

```matlab
addpath('tests')
run_release_tests
```

## Available outputs and interfaces

**The unified spectral API supports orders one to three. Fourth-order
components and the paper's fourth-order comparison have separate entry points.**

| Order | Components | Interface | Scope |
|---:|---|---|---|
| 1 | `eta11`, `psi11` | Unified spectral API | Linear |
| 2 | `eta22`, `psi22` | Unified spectral API | Positive pure sum; prescribed GL rank for elevation, dual-branch GL2+2 for surface potential |
| 3 | `eta33`, `psi33` | Unified spectral API | Positive pure sum; bounded GL4 reconstruction |
| 4 | `eta44`, `psi44` | [`gl_pure_sum_order4`](src/gl_pure_sum_order4.m) | Experimental positive pure-sum component interface; dimensionless analytic inputs and outputs |
| 4 | Paper's `eta44` comparisons | [`paper/order4`](paper/order4/README.md) | GL6/8/10 versus WIT on the declared paper cases |

The unified API follows the two-stage MF12 calling pattern:

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

Passing `order=4` to `gl_spectral_coefficients` is not supported. The
fourth-order paper comparison concerns elevation; it is not a general
validation of all fourth-order surface-potential cases.

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

The [separate nonzero difference-frequency `eta20` diagnostics](diagnostics/eta20/README.md)
include:

- shared-scale GL6/GL12/GL16 rank diagnostics;
- exact frozen Neumann R2/R4/R6 formulas;
- independent ordered-pair field validation for the R sequence.

These diagnostics are not silently added to the total field returned by
`gl_spectral_surface`. R2/R4/R6 are not Green--Laplace ranks, and their
fixed-FFT production implementation is not claimed as validated.

## Dependencies by task

| Task | Dependencies |
|---|---|
| Run the example, release tests, or preserved fourth-order plots | MATLAB R2022b; base MATLAB only |
| Recompute fourth-order GL fields against the preserved WIT reference | MATLAB R2022b; base MATLAB only |
| Build and rerun the fourth-order C++ comparison | Linux or WSL, C++17 compiler, OpenMP, FFTW3 with its threads library, GNU patch, Python 3 for process orchestration |
| Run GL--MF12 comparison examples | MATLAB plus an external MF12 MATLAB implementation |
| Regenerate exact symbolic interfaces | Wolfram Mathematica / WolframScript |

The release tests use base MATLAB only. No Python package or MATLAB toolbox is
required. The current local release gate was executed with MATLAB R2022b;
newer releases are not claimed until the same checks run there.

MATLAB also provides independent validation of C++ results. See
[dependency details](docs/dependencies.md) and the
[fourth-order build instructions](paper/order4/README.md#build-and-check-c).

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

For the fourth-order comparison, run from the repository root:

```matlab
setup_green_laplace
addpath('paper/order4')
reproduce_order4_fields(false) % plot the preserved GL and WIT fields
reproduce_order4_fields(true)  % recompute GL6/8/10 against the preserved WIT reference
reproduce_order4_runtime      % plot historical measured timing data
```

Outputs go to `results/order4/`. The [fourth-order guide](paper/order4/README.md)
also explains how to compile and run C++ on the supplied inputs. Plotting
historical timing measurements does not rerun the 32-core timing campaign.

The source and data directories are organized as follows:

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
`run_all_figures` does not run the separate fourth-order package. See the
[paper guide](paper/README.md) for all entry points, comparisons requiring
external data, and the status of historical reference assets.

The [MATLAB release checks](https://github.com/Xinyu01091002/green-laplace-wave-reconstruction/actions/workflows/matlab.yml)
and [fourth-order reproduction checks](https://github.com/Xinyu01091002/green-laplace-wave-reconstruction/actions/workflows/order4-paper.yml)
run in GitHub Actions. The latter check frozen inputs, compile C++, recompute
the paper GL fields and compare C++ outputs with MATLAB.

The original extraction manifest is supplemented by the
[no-Stokes migration manifest](docs/no_stokes_source_manifest.json) and the
[fourth-order source manifest](paper/order4/source_manifest.json).

## Scope

The released nonlinear graphs accept first-order surface elevation only.
The unified API requires strict-forward analytic support and alias-safe FFT
grids. Its order-two surface-potential graph requires every parent `kh >= 0.3`;
order three requires every parent `kh > 0.5`. Its `Ux` and `Uy` slots must
both be zero. Separate research and paper interfaces have their own declared
domains; these unified API restrictions must not be assumed to describe
every research case. Evidence outside the declared domains is not a release claim.

This software is a fixed-order computational reformulation of regular bound
wave interaction kernels. It is not a solver for strict-zero modes, mean
flow, resonant primary-harmonic corrections or coupled nonlinear evolution.

## Version 0.3.0

**Version 0.3.0 uses the pure GL graph without Stokes-diagonal repairs.**
The third-order API and paper entry points omit repairs in both the nested
second-order states and the third-order response. Results therefore differ
from earlier corrected versions; an exact monochromatic diagonal is not
imposed. The Taylor terms needed to evaluate potential at the free surface
are retained. See the [changelog](CHANGELOG.md) and
[migration record](docs/no_stokes_migration.md).

## License and citation

Source code is released under the [MIT License](LICENSE). Cite the repository using
[`CITATION.cff`](CITATION.cff); cite the Green--Laplace paper and MF12 theory separately when
using the associated methods or comparison implementation.
