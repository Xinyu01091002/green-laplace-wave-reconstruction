# Second-order GPU RHS checkpoint — 27 September 2026

Follow-up: BAND_DIAGNOSTIC.md isolates out-of-band roundoff amplification.
Explicit predeclared input/product projections restore large-grid FP32
agreement in these sparse tests. The failures below remain the original
unprojected baseline; they do not imply FP32 is intrinsically unsuitable.

The P40 executes the finite-depth quadratic RHS in FP32 and FP64. FP64
passes the current analytic-convolution tests through 2048x1024. FP32 passes
smaller grids but **fails** the declared tolerance on the large grid; do not
promote the current FP32 path to production.

## Snapshot and execution

Source snapshot: `/root/hos-rhs-check-20260927T003325Z`, jfm92, sm_61,
CUDA 12.9.86, GNU 13.3.0. Local log and manifest:
`artifacts/hos_gpu/20260927T003325Z/{run.log,SHA256SUMS}`.

All five source hashes were verified on transfer. `build` is the tested
default CMake build. A separate `build-release` also compiled successfully,
but its executable was not used for the reported results; no speed claim.

| Source | SHA-256 |
|---|---|
| CMakeLists.txt | 4619880fac702c50fc0a302e331ad9db9ecfce1226e341964bec2cb69d5d47ac |
| device_fft.cuh | b93301b79d9c97cade2e866abcea276aa7c2c534c6e55aade815ce07fb8a9d56 |
| fft_check.cu | 4a31ac179d96b196a1023ada576e9cb3e582af105add7edd9361cc0a142f6af8 |
| rhs_check.cu | 72f3598471f08b63815998543baa8c518d6047a0c620798b5375b8197948f4bd |
| second_order_rhs.cuh | 50411639dc3c7ffd50f887484a3dca1a1e8bc564357df462ff64a4a7f92efe4d |

## Results

Six CTest tests passed: two FFT tests and four RHS test executables
(FP32/FP64 in 1D and 2D). Each RHS executable checks three depths.
Two additional 256x128 RHS invocations passed all three depths.

Below, maxima are across the three tested depths unless otherwise stated.
Errors are absolute instantaneous RHS errors, not relative wave-height errors.

| Grid | Precision | Max eta_t error | Max psi_t error | Result |
|---|---|---:|---:|---|
| 64x1 | FP32 | 4.4102e-7 | 1.1555e-7 | PASS |
| 64x1 | FP64 | 1.1276e-15 | 5.5512e-16 | PASS |
| 64x32 | FP32 | 7.8631e-7 | 1.6961e-7 | PASS |
| 64x32 | FP64 | 2.0956e-15 | 6.6614e-16 | PASS |
| 256x128 | FP32 | 6.2270e-6 | 7.2373e-7 | PASS |
| 256x128 | FP64 | 1.4218e-14 | 1.5544e-15 | PASS |
| 2048x1024 | FP32 | 4.3543975846449534e-4 | 3.2765125880240831e-6 | FAIL at h=0.15 m; later depths not run |
| 2048x1024 | FP64 | 6.63136212608606e-13 | 1.8762769116165146e-14 | PASS all three depths |

FP32 tolerance 2e-5, FP64 tolerance 2e-11; neither was relaxed. Large-grid
invocations were sequential and bounded by 100 seconds per executable.
FP32 exited 1, FP64 exited 0. Successful cases also reproduced the same
result exactly when evaluate() reused the existing state and work buffers.

Large-grid FP64 raw errors at h=0.15,1.3,20 m respectively:

- eta_t: 6.5800143111971465e-13, 6.6152638922289952e-13, 6.63136212608606e-13.
- psi_t: 1.3877787807814457e-14, 1.8762769116165146e-14, 1.8707257964933888e-14.

Declared arrays plus cuFFT workspace at 2048x1024: FP32 92,291,072 bytes
(88.016 MiB), FP64 184,598,528 bytes (176.047 MiB). This is the second-order
RHS only; it excludes CUDA context/library overhead and future M=5/RK arrays.

The increase in FP32 error with working-grid resolution is consistent with
roundoff in physical-to-spectral transforms being amplified by derivatives
and cancellation. This is an interpretation, not an isolated root-cause
experiment. Inputs start as sampled physical arrays; no spectrum projection
or cutoff was introduced to tune these scores. The test does not establish
that all FP32 formulations fail, or determine an acceptable HOS time step.

Next gates: production input-band projection and matching HOS-Ocean masks,
higher-order recursion, independent CPU HOS comparison, then time integration.
There is no M=5, long-time accuracy or acceleration claim at this checkpoint.
No hardware settings were changed.
