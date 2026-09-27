# P40 CUDA FFT smoke test — 27 September 2026

Six numerical checks passed on jfm92, Tesla P40, 24 GiB. This establishes
basic CUDA/cuFFT functionality only, not whole-device health, sustained-load
stability, complete HOS correctness, long-time FP32 adequacy or speedup.

## Reproducibility

- Remote source/build: `/root/hos-fft-check-20260927-v1`.
- SSH: root at 192.168.2.92:22 via root at 60.188.112.99:60093;
  existing workstation key, no key contents copied to remote hosts.
- Driver 570.211.01; nvcc 12.9.86; GNU 13.3.0; target sm_61.
- Three source SHA-256 values match BUILD_CHECK.md (CMakeLists.txt,
  device_fft.cuh, fft_check.cu). CUDA build and link passed.
- `ctest --test-dir build --output-on-failure`: FP32/FP64 64x32 both
  **passed**, neither skipped. Raw results: remote build/Testing/Temporary/LastTest.log.
- Four additional invocations of `build/hos_fft_check`: each precision on
  2048x1024 and 1024x1. Each reported pass=true and process exit 0.

## Recorded errors and allocated bytes

The input is a DC component plus one diagonal cosine. Errors are maximum
absolute errors in this test's dimensionless units. Spectral errors compare
normalised FFT coefficients against their analytic values. Tolerances were
set before execution: FP32 2e-6; FP64 2e-13.

| Grid | Precision | Round-trip error | Spectrum error | Array bytes | FFT workspace bytes |
|---|---|---:|---:|---:|---:|
| 2048x1024 | FP32 | 5.3644180297851562e-7 | 2.9811398097683065e-8 | 16785408 | 8388608 |
| 2048x1024 | FP64 | 8.8817841970012523e-16 | 5.1321196989281621e-17 | 33570816 | 16793600 |
| 1024x1 | FP32 | 2.6822090148925781e-7 | 1.8637792425163663e-8 | 8200 | 4096 |
| 1024x1 | FP64 | 4.4408920985006262e-16 | 2.6220609358837407e-17 | 16400 | 8192 |

Array plus workspace storage for the large transform is about 24.01 MiB
FP32 or 48.03 MiB FP64. This excludes CUDA context/library overhead and is
not the memory footprint of a full HOS solver. No timing benchmark was run.

## Read-only device observations

At 00:04:13 UTC before testing: 29 C, 9.80 W, 5 MiB, 0% utilisation;
no compute processes. At 00:05:13 UTC afterwards: 30 C, 9.61 W, 5 MiB,
0% utilisation. These are endpoint observations, not sampled peak readings.
No matching NVRM/Xid/AER/PCIe Bus Error entries appeared in the current-boot
kernel journal after the pre-test timestamp.

Power limit remained 250 W; reported slowdown/shutdown temperatures remained
92/95 C. ECC remained disabled. No settings, firmware, clocks, power limits,
thermal protection or driver state were changed; no GPU reset occurred.

The preceding read-only NVML check returned code 14 (Corrupted infoROM) from
both validation and configuration-checksum queries. Passing these arithmetic
checks does **not** repair or clear that condition. No repair was attempted.
