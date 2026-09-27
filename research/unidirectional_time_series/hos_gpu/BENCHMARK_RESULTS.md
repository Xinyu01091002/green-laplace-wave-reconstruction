# P40 / 8-thread CPU prototype benchmark — 27 September 2026

All numerical gates passed. The current GPU prototype measured 35.6x FP32
and 17.2x FP64 steady RK4-step acceleration at 2048x1024 against its matched
8-thread FFTW/OpenMP CPU backend on the same jfm92 host. These are NOT ratios
against production HOS-Ocean MPI or proof of equal long-time accuracy.

Hardware: Tesla P40 24 GiB, driver 570.211.01, CUDA 12.9 compiler; CPU reported
Intel Xeon E5-2678 v3 @2.50 GHz under VMware, 32 visible logical CPUs.
Both benchmark processes were bound to CPUs 0..7. CPU FFTW and pointwise
loops use 8 threads, not 8 MPI ranks. See BENCHMARK_METHOD.md for definitions.

Source snapshot `/root/hos-benchmark-20260927T011541Z`; local evidence under
`artifacts/hos_gpu/benchmark-20260927T011541Z/`. Source SHA-256, individual
samples, correctness results and derived ratios are in benchmark_results.json.
Raw output fields remain remote under `measurements/`; compact JSON was fetched.
The earlier 011410Z snapshot failed compilation due to a missing cstring
include; no benchmark results from that snapshot are used.

## Steady timings

M=5, 273 positive-x directional modes plus conjugates; eta RMS .01 m,
Lx=90 m, Ly=60 m, h=1.3 m. Identical prepared inputs, complete polynomial
support, same RK4 step dt=.005 s. Median of three measurements after warm-up.
No transfer or reference work in the steady timing; GPU is synchronised.

| Grid | Precision | CPU RHS (s) | GPU RHS (s) | RHS speedup | CPU RK4 step (s) | GPU RK4 step (s) | Step speedup |
|---|---|---:|---:|---:|---:|---:|---:|
| 512x256 | FP32 | .0848467 | .00364714 | 23.26x | .356235 | .0153243 | 23.25x |
| 512x256 | FP64 | .109104 | .00799824 | 13.64x | .460432 | .0335515 | 13.72x |
| 2048x1024 | FP32 | 1.44914 | .0405681 | 35.72x | 6.08855 | .170924 | 35.62x |
| 2048x1024 | FP64 | 2.12700 | .124362 | 17.10x | 8.96784 | .521580 | 17.19x |

The large-grid CPU step samples span 6.0593--6.0905 s FP32 and
8.9208--9.0440 s FP64; GPU spans .170883--.170924 s and .521201--.521760 s.
Three samples are a bounded timing observation, not a scaling study or a
statistical confidence interval. Both backends still recompute spectral
multipliers, and neither is claimed optimally tuned.

## First evaluation including setup and round trip

Includes operator/RK allocation, planning, input copy/upload, first RHS and
output copy/download. Excludes executable startup and prepared-input generation.

| Grid | Precision | CPU (s) | GPU (s) | Ratio |
|---|---|---:|---:|---:|
| 512x256 | FP32 | 1.22163 | .339913 | 3.59x |
| 512x256 | FP64 | .533929 | .315174 | 1.69x |
| 2048x1024 | FP32 | 5.78312 | .314491 | 18.39x |
| 2048x1024 | FP64 | 3.68089 | .422976 | 8.70x |

CPU FFTW_MEASURE plan selection is charged here, as is cuFFT planning on GPU.
These figures must not be presented as full program startup latency.

## Measured memory

All values MiB (2^20 bytes). RAM is whole-process peak RSS; GPU memory is
device-wide NVML usage sampled every .05 s, including context/library usage.
GPU baseline was 141.9375 MiB (reserved/idle usage included by NVML), with
no compute job seen at preflight. Samples are not an exact instantaneous peak.

| Grid | Precision | CPU process RAM peak | GPU process RAM peak | GPU total VRAM sampled peak | GPU increase above baseline | Declared GPU arrays + FFT workspace, including RK |
|---|---|---:|---:|---:|---:|---:|
| 512x256 | FP32 | 30.88 | 167.22 | 326 | 184.06 | 19.52 |
| 512x256 | FP64 | 50.38 | 167.69 | 342 | 200.06 | 39.03 |
| 2048x1024 | FP32 | 366.04 | 186.98 | 620 | 478.06 | 312.06 |
| 2048x1024 | FP64 | 711.29 | 232.23 | 930 | 788.06 | 624.14 |

Thus this tested prototype fits comfortably in 24 GiB VRAM. This is not an
allocation proof for a production spectrum requiring different padding or
concurrent phase runs. GPU host RAM remains in use; no whole-system memory
reduction is inferred by comparing CPU RAM to GPU VRAM.

## Numerical results

All 24 full-field comparisons passed (six reference/candidate pairs, four
fields each): eta_t, psi_t, eta after one step, psi after one step.
No fitting, rescaling or alignment. Previous eight CTest regressions passed.

At 2048x1024:

| Pair | Largest relative L2 across all four fields |
|---|---:|
| CPU FP32 / GPU FP32 | 7.46351e-7 |
| CPU FP64 / GPU FP64 | 8.88104e-16 |
| CPU FP64 / GPU FP32 | 9.04074e-7 |

The CPU and GPU share physics expressions, so this is backend agreement;
the independent MATLAB-reference evidence remains in HOS_ORDER5.md and
TIME_INTEGRATION.md. A single RK4 step does not establish long-term drift.

Sampled GPU temperature never exceeded 40 C across this measurement matrix.
No hardware settings, power caps, thermal protection or current production
jobs were modified. No inference about the separate InfoROM warning is made.
