# Continue the already running CPU validation

Update 2026-09-27, before the new 80Tp GPU queue: live inspection found no
remaining official_cpu_matrix or HOS validation processes. The low run log
contains BENCH_ADVANCE_SECONDS=1059.085993199 and BENCH_STEPS=138/0, with
final output 27.6 s. The wrapper did not leave summary.json and the high
directory was never created. Thus the older pending-status text below is
historical; do not claim the high CPU validation is still running. Preserve
the low outputs and recover their compact comparisons separately. The new
LOW_80TP_CAMPAIGN.md queue takes priority on 92; do not launch CPU32 over it.

At the 2026-09-27 checkpoint, the isolated CPU32 matrix is still running on 92.
Do not launch a duplicate. The tool exec session was 36561, with remote command:

```
python3 /root/hos-global-zero-audit-20260927T043227Z/official_cpu_matrix.py \
 --runtime /root/hos-official-benchmark-20260927T041158Z/runtime \
 --exe /root/hos-global-zero-audit-20260927T043227Z/HOS-benchmark-globalzero-fixed \
 --out /root/hos-official-benchmark-20260927T041158Z/corrected-long32 \
 --duration 27.6 --ranks 32
```

The script runs low then high sequentially, with a 3600 s per-case watchdog.
Low was still advancing when the GPU 138 s tests finished. All GPU cases have
now ended. Do not overlap new GPU compute with this CPU32 validation: MPI spin
waiting on a shared pinned core caused extreme slowdown. Times from this
validation are excluded from speed claims. CPU8/GPU two-second timing was
completed sequentially before this overlap and is unaffected.

Connect as root to 192.168.2.92:22 through root@60.188.112.99:60093 using the
existing id_ed25519_cursor key. Keep all hardware/device settings unchanged.
Inspect processes, run.log, Results/probes.dat and summary.json before acting.
Only processes belonging to this isolated matrix may be managed; production
simulations elsewhere must remain untouched.

Once closed, run collect_long_evidence.py to refresh long_run_results.json.
It compares completed CPU diagnostics against GPU FP64 and fetches compact
metrics/hashes only. Large modal data remain remote. Check every exit/status:
the wrapper may finish normally even when the scientific executable failed.
If the high case stops, preserve its log and compare the common completed
prefix. GPU high cases stopped after output 17.2 s; their last common full-field
checkpoint is 13.8 s. Do not claim a full 27.6 s pass or long-run speed ratio.

The official MPI global-zero gravity fix exists only in an isolated diagnostic
copy. Its numerical patch and timing-driver communicator fix are documented in
LONG_RUN_BENCHMARK.md. Do not describe it as byte-identical official source.
Update that report's pending CPU32 paragraph after checking the result.
