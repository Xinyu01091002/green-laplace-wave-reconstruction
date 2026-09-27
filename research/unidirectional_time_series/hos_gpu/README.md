# GPU HOS development

The user explicitly requested GPU development on 27 September 2026, superseding
the earlier initialization-only restriction on solver rewrites for this work.
This directory is isolated from current CPU HOS campaigns and public GL APIs.

## Implemented scope

CUDA/cuFFT infrastructure, finite-depth HOS RHS through degree five,
and a fixed-step RK4 validation integrator. **Not a production HOS replacement**:

- FP32 and FP64 device-resident real and half-spectrum buffers.
- Reusable forward/inverse plans with one shared, explicitly sized workspace.
- Device-side inverse normalization; no host transfer inside the transform API.
- Analytic DC/diagonal Fourier coefficient and round-trip checks.
- Array/workspace byte counts from actual plans, not a whole-solver estimate.
- Device-resident quadratic eta/true-surface-psi right-hand sides, with an
  independent sparse Fourier convolution reference. See SECOND_ORDER.md for
  equations, strict input-support contract and production-solver distinctions.
- Optional fixed-support projections and a four-variant precision diagnosis;
  BAND_DIAGNOSTIC.md records the restored large-grid FP32 agreement for the
  declared sparse inputs, while retaining the original failing baseline.
- Degree-one through degree-five Taylor HOS recursion, independently checked
  against MATLAB Fourier convolution; see HOS_ORDER5.md and order5_results.json.
- Fixed-step RK4 for a deliberately specified retained-band Galerkin system;
  see TIME_INTEGRATION.md. This is not yet HOS-Ocean's adaptive integrator.
- A separate M=5/q=3 MPI-compatible RHS adapter, directly checked against
  unchanged production HOS-Ocean numerical objects. OCEAN_ALIGNMENT.md records
  the normalization, partial dealiasing, Nyquist and linear-splitting checks.
- Production-sized low/high JONSWAP phase-zero RHS comparisons pass in FP64
  and FP32 (PRODUCTION_RHS.md). Raw fields stay remote.
- The original linear-split Cash-Karp attempted step and embedded error are
  matched in FP64, including rejected-step retries (OCEAN_CASH_KARP.md).
- The outer adaptive controller matches original step decisions and short
  continuous evolution, including actual low/high phase-zero JONSWAP states
  through 0.3 s; FP32 tolerance sensitivity is recorded in OCEAN_ADAPTIVE.md.
- Same-host CPU8/GPU timing, the isolated official MPI global-zero correction,
  27.6 s tests and completed 138 s low-case precision comparisons are recorded
  in LONG_RUN_BENCHMARK.md. LONG_RUN_HANDOFF.md records the preserved CPU32
  validation outputs and the latest recovered status.
- The added third-harmonic interpolation/low-rank/cache route was removed
  at the user's request, including deployed source bundles and its saved plan.
  Its retained experimental reports are not original-GL performance evidence.
  SECOND_ORDER_IMPLEMENTATION.md records what the earlier second-order
  time-series comparison actually called; it used pairwise GL-kernel evaluation.

Physical array layout is `[ny][nx]` in C row-major order; the half-spectrum
is `[ny][nx/2+1]`. Forward transforms are unnormalised, inverse transforms
apply `1/(nx*ny)`. Plan use is sequential on the default stream. Inverse
transforms can destroy the input spectrum; HOS state must use separate arrays.
There are no power, clock, driver or running-job modifications.

## Build and check

```sh
cmake -S . -B build -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=86
cmake --build build -j2
ctest --test-dir build --output-on-failure
./build/hos_fft_check float 2048 1024
./build/hos_fft_check double 2048 1024
```

Architecture 86 targets the historically recorded RTX 3080 Ti; verify the
actual GPU before deployment. Missing devices return 77 (CTest skip), never
a numerical pass. No GPU performance or numerical validation is claimed
from compilation alone.

## Next implementation gates

1. Run the FFT checks on a visible GPU, including rectangular and 1D grids.
2. Extend the verified initial-state official M=5/q=3 agreement to evolved
   states and additional phases (the original full-support path is retained).
   Relevant upstream source: `resol_HOS.f90`, `variables_3D.f90`, and
   `fourier_r2c_FFTW3.f90` in the deployed HOS-Ocean v2.1.0 source.
3. Extend short adaptive checks to several periods, all phases, and diagnostics
   of phase/energy/harmonic accuracy. The earlier RK4 remains a separate path.
4. Validate short trajectories, then long-record phase and phase-sector
   accuracy in FP64/FP32. FP32 must not inherit a 1e-12 tolerance blindly.

Current production dimensions are 1024x512, expanded to 2048x1024 for its
q=3 product setting, M=5. Do not replace this with a generic dealiasing rule
without reference comparisons. One phase at a time can reuse all buffers;
historical output is streamed to storage rather than retained on device.
