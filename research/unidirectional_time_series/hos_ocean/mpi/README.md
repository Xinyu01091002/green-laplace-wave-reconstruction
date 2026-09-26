# Directional HOS MPI check: one versus four processes

**Later probe-IO update:** The full-spatial-field benchmark below passed,
but a subsequent production check found the upstream MPI probe path omitted
broadcast of nprobes. Use the corrected v2 production binary documented in
research/directional_wave_data/hos_production/LAUNCH.md for MPI probes.
The old binary and spatial-field timing/validation evidence are preserved.

Completed 2026-09-26. This is a short single-case MPI latency/memory and
implementation-consistency check, not a complete 220 s campaign or a
parallel-scaling study.

## Fixed case and comparisons

Same high-precision frozen focused directional input used in the earlier
serial pilot: 1024x256 unique points, kph=1, Akp=.12, M=5, qx=qy=3,
absolute RK tolerance 1e-8, Ta=0, 2 s evolution, output every .2 s.
One phase only. MPI-1 and MPI-4 execute the identical binary/configuration
sequentially, with rank count the sole numerical-layout change. External
initial eta/psi text is split into contiguous y slabs; decimal values are
copied unchanged. Four processes have 64 y rows each.

MPI is Open MPI 4.1.6 with FFTW-MPI 3.3.10. Ubuntu packages were downloaded
and unpacked under this run's private deps directory, without system
installation. The MPI wrapper's relocated library path and PMIx component
path were configured privately. Existing OW3D runtime environments were not
changed. The source is upstream v2.1.0 commit
`4deb3b4913d993c4e6ea16f736e5fc5792e14f12` plus the already validated IO-only
surface/probe precision patch; no evolution formula was changed.

MPI binary SHA-256:
`b00f8d20a4d6ba82c472bb239e8577f1e6c1b5a1c6320edb3257831f75c2edde`.
Initial full field text SHA-256:
`aa910d9f22195ee7e842f54169d760c912b36ccb1c8e535cd090960b79f4de78`.
Package hashes and build configuration are in build-provenance.json.

## Measured resources

| Metric | MPI 1 | MPI 4 |
| --- | ---: | ---: |
| External elapsed time, including launch/output | 135.6721 s | 67.2499 s |
| Sum of HOS rank user+system CPU time | 135.04 s | 267.25 s |
| Peak sampled job RSS, launcher and descendants | 945260 KiB (0.901 GiB) | 1020632 KiB (0.973 GiB) |
| Peak sampled HOS-rank RSS sum | 933084 KiB | 1003456 KiB |
| Peak sampled job PSS, sharing apportioned | 934344 KiB (0.891 GiB) | 959407 KiB (0.915 GiB) |
| Per-rank GNU time maximum RSS | 932576 KiB | 249400--251148 KiB |

Observed elapsed-time speedup: 2.0174x; efficiency versus four slots: about
50.4%. Aggregate rank CPU time increases about 1.98x. Thus this improves
single-task latency, not CPU-hour efficiency. No 8-process result exists.

Memory was sampled at approximately .25 s intervals over the MPI launcher
process tree. RSS sums may count shared mappings more than once; PSS is
also reported. Sampling can miss brief peaks. Sum of individual per-rank
GNU time high-water marks is reported separately and is not an observed
simultaneous aggregate peak. Linux /proc RSS and getrusage accounting can
also differ slightly. The Python monitor and later MATLAB validation are
outside the reported MPI process tree.

OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and MKL_NUM_THREADS were 1; Open MPI
background threads were present but each rank used approximately one CPU.
Ranks were bound to cores and shared-memory/TCP transports were used on one
host. Builds used two compilation workers. The two trials never overlapped;
maximum extra compute-process count was four, with an 8 GiB sampled-job-RSS
abort guard. Sixteen OW3D processes remained running, available RAM stayed
about 472 GiB. No OW3D process or frozen case was modified.

## Numerical validation

MATLAB joined rank-local files in y order and checked their coordinates,
all values, and all eleven physical times 0:.2:2 s. Predetermined relative
L2 acceptance threshold was 1e-10, with no gain/offset/phase/shift fitting.

Maximum MPI-4 versus MPI-1 relative L2 over all frames is below 4e-14 for
both eta and true surface psi. Both were also compared against the existing
serial build: maximum relative L2 approximately 9.97e-12 for eta and
2.08e-12 for psi. All checks passed and all ranks exited zero.

## Use and limitations

Remote run root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-mpi-check-20260926-v1`.
Executable: `build/sources/HOS-Ocean`.
Private launch environment: environment.json; run_mpi_check.py records the
actual mpirun command and rank.sh supplies individual resource logs.
Use this environment rather than adding the private libraries globally.

Native fields and memory time series remain remote; local compact summary:
`artifacts/hos_ocean/mpi-results/summary.json`.
The baseline and four-process job wall times sum to 202.92 s. These are two
physical trials with five solver ranks, not five independent wave cases.
Long-run near-focus behavior, other process counts, multi-job contention,
and final production accuracy/output settings are not measured here.

With a fixed budget of four CPUs, four independent serial phase runs may
provide better total throughput than running each phase sequentially on four
MPI ranks. Running all four phases with four ranks each would require 16
compute slots and needs a separate concurrency decision. This check does
not authorize or launch such a batch.