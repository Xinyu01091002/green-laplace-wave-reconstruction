# Directional HOS resource assessment, 2026-09-26

Read-only remote checks at about 02:57--03:00 UTC; no new solver launched.

Current host: 48 visible CPUs, 747 GiB total RAM, approximately 472 GiB
available RAM, 830 GiB disk available. Two vmstat interval samples showed
about 31% user, 3% system and 67% idle CPU, with zero swap-in/out. Sixteen
single-thread OW3D processes are active. The reported topology is not a
claim of exclusive physical-core ownership.

Active OW3D campaign: ow3d16-kpd1-akp012-20260925T212822Z. Wall elapsed about
5 h 19 min at inspection. Focused-case logs explicitly report physical time
10 s; a random-case log reports 8 s, versus a 220 s target. These logs are
buffered/coarsely reported, so they are lower bounds on current simulation
time, not precise throughput measurements. Current per-process RSS is about
16.3--16.6 GiB, with observed startup maxima 23.38 GiB. No job was modified.

Existing directional HOS evidence: precision-v1/case-precise-full used
1024x256 unique horizontal points, kph=1, M=5, qx=qy=3, Ta=0, absolute
integration tolerance 1e-8, output every .2 s. The same frozen second-order
focused input was used. Two physical seconds completed in 139.82 wall
seconds, peak RSS 910968 KiB = .869 GiB, one thread. Three concurrent HOS
IO-control pilots had similar times and peaks; this is not a scaling study.

Planning estimates, per phase (not measured complete-run timings):

| Configuration | Target duration | RAM estimate | Wall-time estimate |
| --- | ---: | ---: | ---: |
| Compact 1024x256 directional grid | 220 s | about 1--1.5 GiB budget | about 3--5 h |
| Legacy small-grid 1024x256 geometry | 360 s | about 1--1.5 GiB budget | about 5--8 h |
| Legacy large-grid 2048x512 geometry | 480 s | about 3.5--5 GiB, reserve 6 | roughly 1--2 days |

The compact estimate is anchored to 139.82/2*220=15380.2 s (4.27 h), with
startup/output effects apparent in the 2 s pilot. The large-grid estimate
uses approximately four times the field storage and a four-to-five-fold
FFT work factor, then the longer physical duration. These estimates assume
comparable numerical settings and similar timesteps. They do not certify
long-run peak memory, accuracy, near-focus timestep costs, or speedup versus
OW3D. M, dealiasing and integration tolerance changes require updated budgets.

Legacy successful local references inspected: directional_joint_input's
small grid is 1024x256 over 9000x6750 m; the Akp=.12 amplitude-qualified
test6 is 2048x512 over 22500x9000 m, with a 480 s native OW3D run. They must
not be conflated with the compact 50lambda x 20lambda, 220 s campaign.

Suggested resource policy: first one phase for a 5--10 s physical-time
resource check with final probe output/tolerance, then four independent
single-thread phases. Under the tested compact setup, reserve 1.5 GiB per
process (6 GiB total) and four CPUs, while keeping the 16 OW3D jobs intact.
Set OMP_NUM_THREADS/OPENBLAS_NUM_THREADS/MKL_NUM_THREADS=1. Run MATLAB
preparation and reconstruction separately with -singleCompThread and reuse
frozen inputs where possible; their memory is not included in the HOS budget.
Do not increase to 8 or 16 HOS processes until throughput, RAM, CPU and IO
under four concurrent jobs are observed. Independent-case parallelism is the
first choice; no MPI speedup is established by the serial pilots.

Save center/off-center probe records at adequate temporal sampling and sparse
spatial snapshots. Full double-precision ASCII eta/psi fields on 1024x256 at
.2 s over 220 s are about 15 GB per phase, before other files. Probe output
avoids that disk cost but does not eliminate FFT evolution cost. The enlarged
HOS probe formatting patch is already validated in the unidirectional case.