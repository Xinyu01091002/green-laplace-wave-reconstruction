# MPS replay and matched four-process timing

## Completed 1100.8 s MPS result

MPS full batch completed with measured wall time 9527.727701553 s
(158.79546 min, 2 h 38 min 48 s). Each phase completed 5504 accepted steps,
zero rejections and 5505 diagnostic records. All 20 probe records over the
whole interval, energy-drift records, and all saved eta/psi checkpoints
(including the final 1100.8 s state) have zero measured difference from the
prior completed ordinary campaign. The dedicated MPS service shut down.
Completion evidence: local `artifacts/hos_gpu/mps-20260927T090455Z/complete/`.

| Phase | Original ordinary integration wall (min) | MPS integration wall (min) |
|---|---:|---:|
| 0 | 178.3251 | 158.7776 |
| 90 | 195.5460 | 158.7747 |
| 180 | 195.7849 | 158.7862 |
| 270 | 195.7225 | 158.7564 |

The original ordinary workflow had a phi000 solo head start. Its total interval
from first solver start to the last phase's measured completion reconstructs
to about 12209.92 s (203.50 min), using recorded /proc start_ticks, verified
CLK_TCK=100, and setup+integration times. This excludes small process teardown
overhead and is approximate. Comparing actual workflows suggests 44.70 min
less waiting with MPS, but this is NOT the pure simultaneous-start comparison.

The missing simultaneous-start full-length benchmark is now completed on 92:
`/root/hos-ordinary-full-20260927T123931Z` (controller PID at launch 36875).
All four phases started from zero together, using the SAME batch implementation,
executable, inputs, CPU affinities, output cadence and 1100.8 s duration.
MPS was confirmed stopped and the GPU idle before its launch.

Ordinary full batch: 12257.276572942 s (3 h 24 min 17 s).
MPS full batch: 9527.727701553 s (2 h 38 min 48 s).
Matched full-record speedup: 1.28648477x. Waiting-time reduction: 22.2688%,
or 45 min 29.55 s. All probe, energy and stored field comparisons pass with
zero measured difference; each phase has 5504 accepted and zero rejected
steps. Compact evidence is in
`artifacts/hos_gpu/ordinary-full-20260927T123931Z/complete/`.
These are one complete run per mode, not a multi-run statistical benchmark.
Deployment script: deploy_full_ordinary_comparison.py.

User explicitly requested MPS and a repeat of the same four-phase data on
2026-09-27. Active root on 92: `/root/hos-mps-20260927T090455Z`.
Detached controller PID at launch: 29652. Local metadata:
`artifacts/hos_gpu/mps-20260927T090455Z/location.json`.

## First matched result and full replay launch

All four 0--10 s cases completed in 112.42255 s with ordinary processes and
88.51095 s with MPS: batch speedup 1.27015x, waiting time reduction 21.27%.
Each phase took 50 accepted steps and zero rejections in both modes. All
compared probe records and stored modal checkpoints had zero measured
difference. This is a single matched pilot, not the completed 80Tp result.

MPS server PID 29769 reported all four pilot client PIDs. After the numerical
gate passed, the full replay launched client PIDs 29860--29863, all verified
on that same server. At 09:09:06 UTC each full phase had reached 3.4 s with
17 accepted steps and no rejections. Full completion and timings must be
read from live status.json, not inferred from this launch checkpoint.
Compact pilot evidence is in local status.json/pilot-comparisons.json under
the artifacts directory above; the status file is a snapshot, not live state.

`run_mps_comparison.py` uses the exact same hashed FP64 executable and initial
states as the completed low campaign, tolerance 1e-12, M5/q3, the same grid and
0.2 s diagnostics. Four new processes are launched together, pinned to CPU
28/29/30/31 in BOTH configurations. No previously advanced state is adopted.

Execution order:
1. Ordinary four-process execution, each phase 0--10 s.
2. Launch an isolated legacy MPS service with run-local pipe/log directories.
3. Four MPS clients, each phase 0--10 s. Require the MPS server's client list
   to contain all four solver PIDs; do not accept silent non-MPS execution.
4. Compare full probes, energy drift and stored modal checkpoints with the
   ordinary run at 1e-12 relative field tolerance, plus identical step counts.
5. Only after these checks pass, run all four phases from zero to 1100.8 s
   with MPS and compare against the previous completed ordinary campaign.
6. Shut down only this run's MPS service when the controller closes.

`status.json` is authoritative. `pilot-comparisons.json` and
`full-comparisons.json` retain comparator output. Each batch has independent
phase logs, diagnostics, modal checkpoints and solver summary. Full successful
records also receive compatible Results/probes.dat files.

Timing is four-case batch wall time, including process setup and diagnostics,
ending when all four have exited. It excludes post-run comparison. MPS service
startup is outside the batch timer, but first client/context initialization is
inside. Pilot is one matched trial per mode, not a statistical speed guarantee.
Do not compare a sum of process wall times against one batch wall time.

The preceding completed ordinary long campaign gave phi000 a solo head start
before adoption into the parallel queue. It is useful as a numerical reference
but its parallel-stage elapsed time alone omits work and is not an exactly
matched simultaneous-start full-batch timing. Use the new matched 10 s trial
for the direct ordinary/MPS comparison and report full timing with this caveat.

MPS is software scheduling only. GPU compute mode stays Default. No clock,
power limit, thermal protection, ECC or driver setting is changed. The initial
check requires no existing MPS service, and the GPU was idle before launch.
93's existing GL postprocessor continues independently; see POSTPROCESS_GPU_LOW.md.
