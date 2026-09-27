# HOS-Ocean GPU: same-host timing and longer evolution

Measured on 2026-09-27 on jfm-92: Tesla P40, Xeon E5-2678 v3 VM.
Initial conditions are the actual low/high JONSWAP phase-zero states, base
1024 x 512, extended 2048 x 1024, HOS M=5, production q=3.
These observations do not certify other phases, grids, spectra, or long records.

## Same-host two-second timing

CPU uses eight MPI processes, one thread each, pinned to CPUs 0--7.
GPU cases ran sequentially after CPU timing finished. Each entry is one run.
The CPU executable uses the production Fortran objects, with the isolated
global-zero correction described below. CPU state precision is FP64 and its
absolute integration tolerance is 1e-12 throughout.

| Case | CPU FP64 advance (s) | GPU FP64, tol 1e-12 (s) | Same-precision speedup | GPU FP32, tol 1e-9 (s) | Combined speedup |
|---|---:|---:|---:|---:|---:|
| Low | 52.5113 | 5.69849 | 9.215x | 1.91048 | 27.486x |
| High | 199.7153 | 21.92476 | 9.109x | 1.91166 | 104.472x |

CPU and GPU FP64 take identical accepted/rejected counts: low 11/0, high
30/13. GPU FP32 takes 11/0 for both. The large high-case combined ratio includes
the benefit of fewer attempted steps at the looser tolerance; it is not a
same-precision hardware acceleration claim.

The measured boundary is the time-advancement loop with each implementation's
probe/energy diagnostics, excluding initialization. GPU checkpoint copies and
writes are measured separately and subtracted (about 0.03 s per case). CPU
full-field output is disabled. The implementations have different diagnostic
algorithms and log output; these are application-loop ratios, not isolated RHS
or FFT ratios. GPU inputs are preconverted modal data, whereas CPU reads ASCII
initial fields: do not use whole-process times as a matched startup benchmark.

Two-second CPU/GPU FP64 probe relative L2 differences are at most 3.43e-13
(low) and 3.40e-12 (high). Maximum energy-drift differences are 4.83e-15 and
4.62e-10, respectively. Agreement with the implementation does not establish
convergence to the physical solution.

## Memory

| Quantity | FP64 GPU | FP32 GPU |
|---|---:|---:|
| Sampled total device memory | 968 MiB | 650 MiB |
| Device baseline before run | 142 MiB | 142 MiB |
| Increment above baseline | 826 MiB | 508 MiB |
| Host process maximum RSS | about 151 MiB | about 149 MiB |

CPU8 aggregate sampled RSS is about 1.92 GiB. This is a sum of process RSS,
so shared pages can be counted more than once; it is not unique physical RAM.
NVML sampling is 50 ms for the short GPU runs and 100 ms for the first long
matrix. These are observed peaks, not an allocation upper bound. All evolution
state and FFT work buffers remain on the device; scalar control/diagnostics
and occasional checkpoint exports cross to the host. Historical outputs stream
to disk rather than accumulating in device memory.

## Longer runs

The first GPU matrix requested 27.6 physical seconds (two peak periods).
Low FP64 and FP32 both completed with 138 accepted steps and no rejections.
Uncontended advancement times were 70.1809 and 23.6539 s, respectively.
At the end, FP32 relative full-field L2 differences from FP64 were 2.646e-6
for eta and 2.601e-6 for psi; dominant-mode phase differences were below
4.36e-7 rad. Maximum probe relative L2 was 1.431e-6 and maximum difference
in signed relative energy drift was 5.637e-7.

Both precision variants of the high case stopped at the original HOS
slope guard, after the last completed output at 17.2 s. This is an incomplete
27.6 s test. Its elapsed time must not be used as a successful long-run timing.
Through 17.2 s the two GPU precisions have probe relative L2 differences below
1.701e-6. Their latest common full-field checkpoint is 13.8 s, where eta/psi
relative differences are 2.247e-6/1.763e-6. These comparisons do not identify
the cause of the subsequent stop.

The low case was extended to 138 s (ten peak periods). Both precisions completed
with 691 accepted steps and no rejections. At 138 s the FP32 eta/psi full-field
relative L2 differences from FP64 are 1.2863e-5/1.2637e-5. The five probes have
maximum relative L2 difference 1.055e-5 and maximum absolute difference
1.161e-5 m. Maximum difference in relative energy drift is 2.822e-6; maximum
absolute drift is 3.579e-9 for FP64 and 2.824e-6 for FP32. The final dominant
mode phase differences are -1.525e-8 rad (eta) and -8.991e-7 rad (psi).
This is a cross-precision comparison, not an independent CPU validation over
all ten periods or a physical convergence study.

CPU32 validation is still running at this report checkpoint. Its overlap with
the 138 s GPU cases caused severe MPI scheduling contention. Neither those
CPU32 times nor the 138 s GPU times enter the uncontended timing table above.
The CPU32 runs have a 3600 s per-case watchdog; a watchdog stop is an incomplete
test, not a scientific failure. See LONG_RUN_HANDOFF.md for continuation.

## Isolated corrections and provenance

In the deployed official MPI source, solveHOS_lin subtracts
`g_star*a_etark(1,1)` from local `(1,1)` on every rank. On ranks with
`local_y_start != 0`, that is a nonzero global transverse mode, whose linear
gravity is already in the analytical linear propagator. The diagnostic fix
restricts this subtraction to the global spatial zero mode. See
`mpi_global_zero_gravity.patch`. With this single numerical change, MPI8
matches the GPU/global-zero convention to the probe errors reported above.
Earlier original-object checks used the MPI executable at one rank and could
not expose this decomposition issue. “Official agreement” must state rank
count and whether this correction is present.

The timing driver also supplies both required arguments to MPI_BARRIER
(`MPI_COMM_WORLD, Statinfo`). The first driver had inherited a one-argument
call and crashed; that failed run is retained and excluded. Production source,
executables, and jobs were not modified. The CPU energy output stores absolute
relative drift, whereas GPU diagnostics retain signed drift; the comparator
accounts for this format difference without changing either result.

Remote immutable GPU source/executable: `/root/hos-long-20260927T040628Z`.
CPU runtime/results: `/root/hos-official-benchmark-20260927T041158Z`.
Corrected CPU audit: `/root/hos-global-zero-audit-20260927T043227Z`.
Use `collect_long_evidence.py` for compact metrics, hashes and comparisons;
raw modal checkpoints remain remote. No device protection, clock, power,
driver, or ECC setting was changed.
