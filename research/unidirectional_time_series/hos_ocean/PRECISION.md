# HOS-Ocean high-precision surface IO

This extension uses upstream v2.1.0 commit
`4deb3b4913d993c4e6ea16f736e5fc5792e14f12`.
The accompanying `precision.patch` changes exactly four lines in two files:

- The external surface eta/psi reader uses list-directed input, so it can
  read both the legacy ES12.5 values and full-precision decimal values.
- The ordinary physical free-surface output uses ES25.16E3 for coordinates,
  eta/psi and frame time (17 significant digits). The time-format label is
  also shared by modal frame headers.

The 67-line external-input header and its fixed-width zone record are retained.
The MATLAB exporter writes that header and uses 18 significant digits for the
eta/psi data. This patch targets the external `Initial surface quantities`
route and ordinary physical output; it does not upgrade modal values, refined
output, probe files, diagnostic files, or certify other restart/import routes.
The HOS equations, zero-mode handling, FFTs, dealiasing and integration are
unchanged. Original upstream sources and the official prebuilt binary remain
available separately.

## Build and controls

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0/precision-v1`.

A private GNU Fortran 13.3 toolchain was unpacked from Ubuntu packages;
no system installation or global environment change was needed. Both builds
use identical compiler, Release flags, FFTW, BLAS/LAPACK and MPI/HDF5 settings
(OFF). CMake's upstream Release flags, including fast-math, are unchanged.
Builds use two compilation workers. Each simulation uses one thread.

- Baseline binary SHA-256:
  `4abb84906c97f9c57d66be947beb6a514cc943acc53ba1db927968447af38c0a`.
- Patched binary SHA-256:
  `21a171f03c32a49ae05b62e4896d5b081c26a895372552d6ba563ba8fd44f9b4`.
- Patch SHA-256:
  `f8bdb88987a95bdd03a548bd7047d101e72ce5c16724a981346b9348d6304d6f`.

`build_precision.py` is the corrected clean-build procedure. During the actual
build its first output-format anchor assertion stopped before modifying
output.f90; the bounded `finish_precision_build.py` fixed the subroutine anchor
and completed the build. Both execution scripts and logs are retained remotely.
The resulting diff was inspected: four replacements, with no equation changes.

The three 2-second pilot cases use the same frozen wavegroup, phase zero,
1024x256 grid, M=5, qx=qy=3, RK tolerance 1e-8, 0.2 s output, no ramp or
breaking/dissipation. The cases distinguish compiler and IO/input effects:

1. `case-baseline-low`: unmodified source build, original quantized input.
2. `case-precise-low`: patched source build, identical quantized input.
3. `case-precise-full`: patched source build, original double-precision input
   serialized with 18 significant digits.

The MATLAB validation checks every field value in all eleven frames and
uses thresholds defined before reading results: raw t=0 relative L2 <1e-12,
coordinate error <1e-9 m, and same-input output agreement within the legacy
ES12.5 half-unit rounding bound plus 1e-11 times the field peak. It does not
fit gain, bias, signs, phase, shifts or alter raw output. The full-versus-low
input trajectory difference is reported separately without an accuracy claim.

## Verified results: 2026-09-25 22:35 UTC

All three solver processes exited with code 0. MATLAB verified every eta/psi
value in the eleven frames from 0 to 2 s and all preregistered checks passed.

| Metric | eta | true surface psi |
| --- | ---: | ---: |
| High-precision input text round-trip relative L2 | 0 | 0 |
| Raw t=0 relative L2 against original doubles | 8.6465986343e-16 | 8.0356447136e-16 |
| Raw t=0 maximum absolute difference | 1.5543122345e-15 m | 2.8421709430e-14 m2/s |
| Full-versus-quantized input relative L2 at t=2 s | 2.6897707472e-6 | 1.3221542504e-6 |

Coordinate maximum absolute differences were 1.81899e-12 m (x) and
9.09495e-13 m (y). Same-input original/patched outputs passed the pointwise
legacy-rounding bound throughout. Their reported relative differences are
approximately 1.5e-6 after the initial frame, reflecting legacy output
quantization; they are not an evolution-accuracy estimate.

| Case | Wall time | Peak RSS |
| --- | ---: | ---: |
| baseline-low | 136.16 s | 910980 KiB |
| precise-low | 136.50 s | 911016 KiB |
| precise-full | 139.83 s | 910968 KiB |

These cases ran concurrently on three single-thread processes alongside the
sixteen OW3D runs. The figures include startup and output and are not a scaling
benchmark or a 220 s runtime prediction. All three HOS pilots have completed.

Remote `precision-validation.json`, per-case `status.json` and `resources.txt`
contain the result. `audit-manifest.json` records executed script/patch/report
hashes. Validation report SHA-256:
`777f63de9505c8275a81f85978a307be31c4135f7e9eab1e1dc33afd8e30aff7`.
A small JSON summary was copied locally; raw fields remain remote. The full
precision strip is `case-precise-full/surface_strip.mat`.

Use `build-precise/sources/HOS-Ocean` for subsequent external-surface trials
with the same private BLAS/LAPACK paths. The identified IO limitation is now
resolved for this route. Long-time stability, convergence and high-order
physical accuracy remain unvalidated; no 220 s HOS production family has
been launched in this precision trial.