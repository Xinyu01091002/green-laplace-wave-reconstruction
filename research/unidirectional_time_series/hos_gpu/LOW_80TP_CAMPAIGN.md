# GPU replay of the 93 JONSWAP low-family task

## Verified completion

Live status checked after 2026-09-27 08:49:22 UTC: all four low-family phases
completed 1100.8 s in FP64 at tolerance 1e-12. Each has 5504 accepted steps,
zero rejected steps and 5505 finite probe records with the expected final
time. The controller validated and hashed each Results/probes.dat. No GPU
compute processes remained. The concurrent stage lasted about 3 h 16 min;
phi000 had an earlier solo segment, so this is not a simultaneous-start
four-case benchmark. Peak sampled total device memory was 3,607,101,440 bytes.
Compact completion evidence is saved in
`artifacts/hos_gpu/campaign-20260927T052004Z/completed-status.json`.
This completes low-family HOS evolution and output integrity checks only:
four-phase extraction, GL comparison, and high-family long evolution are not
completed by this run. The execution notes below describe its history.

The user authorized using the GPU implementation for the preceding task on 93
on 2026-09-27. The selected source is the low family of
`hos-jonswap-r4gl-80tp-20260926-v2`, not the older 220 s narrowband experiment.
The CPU campaign remains untouched.

## Current execution: four phases concurrently

The user requested a four-phase concurrent trial after the sequential launch.
`switch_parallel_campaign.py` stopped only the original queue controller,
retained the already running phi000 solver (PID 21872), and started the other
three phases through `run_parallel_campaign.py`. No phase-zero computation was
discarded or restarted. Each GPU host process is pinned to a separate CPU
(31,30,29,28 for phases 0,90,180,270). The new detached controller PID is 22248.

The authoritative live status is now `parallel-status.json`; `status.json`
marks the original sequential controller as superseded. `parallel-adoption.json`
preserves its last state, process identity and the handoff source hash.
`parallel-samples.jsonl` records simulation progress and sampled NVML device
memory every ten seconds for measuring aggregate throughput. All four phases
retain FP64, tolerance 1e-12 and the same 1100.8 s target. A phase failure is
recorded without discarding the independent phases' work. The controller checks
each completed record before writing its compatible probe output. This is
ordinary four-process CUDA execution; no MPS service or hardware setting is
enabled or changed.

Initial concurrent observation: over a 150.03 s wall window after all phases
had progress records, each phase advanced at 0.09065 simulated seconds per
wall second, totaling 0.36260. The prior single-phase average was 0.39289,
so the observed aggregate ratio is 0.9229 (about 7.7% lower). This compares
different phase/time positions, not identical fixed work. Device usage was
3,607,101,440 bytes (3.36 GiB); all phases still had zero rejected steps.
The user-requested four-way execution remains active. Compact evidence:
`artifacts/hos_gpu/campaign-20260927T052004Z/parallel-observation.json`.

The following describes the original sequential setup and its input provenance.

The GPU queue uses the original four phase-shifted initial eta/true-surface-psi
fields, base 1024 x 512, extended 2048 x 1024, M=5/q=3, FP64 and absolute
time tolerance 1e-12. Each phase runs to 1100.8 s with 0.2 s output and 5505
expected records at the original five probes. No input renormalization,
phase regeneration, nonlinear ramp, extra damping or breaking model is added.

Preparation runs the original HOS initialization modules at one MPI rank in a
fresh directory on 93. It joins the eight original ASCII slabs without changing
values, exports the initialized modal state, and does no time integration.
Original slab, module, binary and exported-file hashes are retained. Phase-zero
modal bytes must exactly match the previously validated reference. All four
initial probe comparisons must pass before any full record starts.

The immutable GPU solver is `/root/hos-long-20260927T040628Z/build/hos_run`.
It follows the global-zero convention documented in LONG_RUN_BENCHMARK.md;
do not label it identical to the uncorrected multi-rank official executable.
Its original numerical slope guard and diagnostic attempt budget remain active.
If a phase fails, the queue stops and keeps all partial records and logs.

Inputs on 93:
`/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-campaign-inputs-20260927T052004Z`.
Target queue on 92: `/root/hos-low-80tp-20260927T052004Z`.
Local deployment artifacts: `artifacts/hos_gpu/campaign-20260927T052004Z/`.
Read deployment `location.json` and live remote `status.json` before claiming
launch or completion. The controller is detached with its own session and
redirected logs, so it does not depend on the client tool session staying open.

Launch verified at 05:26 UTC on 2026-09-27: all four 0.2 s preflight runs
completed, and all initial probes passed (maximum absolute difference
1.460e-14 m). Phase-zero initial modal data matched the prior reference
byte for byte. The full phi000 run had reached 7.6 s with 38 accepted steps
and no rejections. This is launch evidence, not a completed 80Tp result.
Controller PID was 21823; the observed first full solver PID was 21872.
Input archive SHA-256:
`45e1cabccb835b81c3fcfda263c6be8833539c8bc3a45ce8625360d4f3c8bb33`.

`status.json` reports the active phase/PID and last completed simulation time.
`controller.log`, `<phase>.log`, `<phase>/steps.csv` and
`<phase>/diagnostics.csv` retain execution evidence. After a complete, finite
record with the expected last time and sample count, the controller writes
`<phase>/Results/probes.dat` for the original four-phase/GL workflow. Native
paths differ from 93, so downstream processing must explicitly select these
GPU records. The GL comparison is not automatically run by this queue.

The high family is not queued for 80Tp: both GPU precisions stopped after the
17.2 s output in earlier testing. Its long-record suitability remains unresolved.
Do not disable the numerical guard or add dissipation merely to finish it.
No hardware, thermal protection, ECC, clock, power or driver settings are changed.
