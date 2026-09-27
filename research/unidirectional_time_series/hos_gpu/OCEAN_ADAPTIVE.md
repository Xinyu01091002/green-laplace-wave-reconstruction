# HOS-Ocean outer adaptive controller

ocean_adaptive.cuh connects the matched Cash-Karp attempted step to an adaptive
loop for the current constant-depth, no-current, no-breaking, no-relaxation
scope. State remains on GPU; step decisions use one returned error scalar.
Controller time and proposed step are FP64 even when the field backend is FP32.

## Rules retained from the production main

- Clip the actual attempted step to the next output time, keeping the proposed
  step separately.
- Reject when error > tolerance; use (error/tolerance)^(-1/4).
- Accept otherwise, commit the candidate and advance time; use exponent -1/5
  when error exceeds FP64 epsilon, otherwise use factor 4.
- Bound that factor to [0.125,4], then multiply by the 0.98 safety factor.
- Update the ORIGINAL proposed step, not the output-clipped actual step.
- Cap the next proposal by dt_out and dt_lin.
- Preserve the source's final-output comparison using epsilon(0.0), i.e.
  default REAL/FP32 epsilon, distinct from the FP64 tiny used elsewhere.
- Retain the numerical slope check at 10. Nonfinite slopes/errors are explicitly
  rejected/stopped rather than allowed to disappear through a max reduction.

No hardware thermal, power, ECC or clock settings are touched. The slope check
is a numerical condition inherited from the solver, not hardware protection.
The diagnostic attempt budget and zero-progress check guard failed tests;
neither silently clips scientific values or changes the acceptance tolerance.

## Official reference

The reference builder extracts the complete original main-program block from
"Going to next time step" through "dt = h_rk", adds only trace logging and a
test attempt budget, and inserts it into a diagnostic helper. The original
controller and instrumented-controller hashes are recorded. Numerical module
objects, including the original Cash-Karp routine, are linked unchanged.
Output/energy-file I/O is omitted from this bounded controller test; breaking
and relaxation are disabled consistently with the current campaign.

The trace records every attempt's starting time, proposed and clipped steps,
embedded error, decision, resulting time and next proposal. FP64 comparison
requires identical decision counts/order, time/step agreement within
1e-10 + 1e-5*abs(reference), and estimator agreement within
5e-15 + 1e-5*abs(reference). These are diagnostic trace gates, not altered
solver tolerances. Saved state errors are measured at identical output times.

State L2 uses the HOS modal normalization's Parseval weights: 1 for zero/x-
Nyquist columns, 1/2 for interior positive-x columns. For the real fields this
is the physical-grid L2 norm, without fitting or alignment. FP64 state gate:
relative L2<=1e-8 and modal max abs<=2e-12+1e-8*reference max. FP32 comparison
against the FP64 reference: relative L2<=1e-3 and max abs<=2e-6+1e-3*reference max.

## Small-grid checkpoint

CPU reference gpu-rhs-probe-20260927T024702Z; GPU snapshot
/root/hos-rhs-check-20260927T025104Z. Records: adaptive_small_results.json.
Low-mode, near-cutoff and Nyquist cases on 64x32 and 128x64, through t=.57
nondimensional with outputs every .13 and a final shortened interval.
The known explosive broad-spectrum stress input from the fixed-step test is
not reused as a physical trajectory. This selection is explicit, not a removed
failure. Six FP64 trajectories at 1e-12 and six FP32 trajectories at 1e-8 pass.

All FP64 accept/reject sequences match the official controller; maximum
observed time/proposal difference is 2.28e-10 and maximum state relative L2
is 8.78e-16. FP32's worst state error is 8.04e-6. FP32 uses a different
tolerance and is NOT expected to follow the FP64 step sequence.

## Production-size bounded test contract

Read existing low/high JONSWAP phase-zero initial states at 1024x512, M=5/q=3.
Run to physical t=.3 s, with outputs at .2 and .3 s. Keep the original main's
initial dt, dt_out and dt_lin; only the bounded end time is selected for this
diagnostic. The reference uses the original 1e-12 absolute tolerance.
GPU FP64 uses the same tolerance and trace gate. FP32 trials use separately
declared 1e-7, 1e-8 and 1e-9 tolerances; compare states to the FP64 reference,
not step-for-step schedules. This does not establish a production FP32
tolerance or long-time phase/energy accuracy.

Raw fields stay remote. A temporary, IP-restricted two-file LAN transfer is
used for the reference trajectories; only compact errors/traces/provenance
are retained locally. Existing production simulations remain unchanged.

## Production-size result — 27 September 2026

Reference snapshot gpu-production-probe-20260927T025257Z on 93 and GPU snapshot
/root/hos-production-check-20260927T025257Z on 92. All initial eight runs passed.
Since the three FP32 tolerances gave identical steps/results, additional 1e-11
and 1e-12 trials were run against the SAME frozen reference to assess whether
stricter temporal control helped. Their state gates were not changed. All four
additional runs passed. Compact source/data provenance, traces and results:
adaptive_production_results.json.

| State | Backend/tolerance | Accepted | Rejected | Maximum state relative L2 over both outputs and fields |
|---|---|---:|---:|---:|
| Low | FP64, 1e-12 | 2 | 0 | 1.51926e-16 |
| High | FP64, 1e-12 | 5 | 12 | 2.86558e-16 |
| Low | FP32, 1e-7 / 1e-8 / 1e-9 | 2 | 0 | 7.61456e-8 |
| High | FP32, 1e-7 / 1e-8 / 1e-9 | 2 | 0 | 7.86996e-8 |
| Low | FP32, 1e-11 | 2 | 0 | 7.61456e-8 |
| High | FP32, 1e-11 | 3 | 2 | 9.19979e-8 |
| Low | FP32, 1e-12 | 3 | 3 | 8.89755e-8 |
| High | FP32, 1e-12 | 13 | 13 | 1.72772e-7 |

Both FP64 decision sequences match the original controller; the high-case
maximum time/proposal difference was 3.72e-12 nondimensional. FP32 is compared
to the FP64 reference states, not required to reproduce its step decisions.
The looser-tolerance plateau is output-interval limited, so these data do not
identify an optimal FP32 tolerance. The strictest FP32 setting increases work
without improving agreement in this short segment. This is a bounded
observation, not a universal precision/tolerance rule.

The high CPU reference repeated several identical clipped steps near the
last output boundary: error about 1.2607e-12 exceeded tolerance, but updating
the larger proposed h_rk left h_loc at the same remaining interval until
enough reductions accumulated. The GPU reproduces this source behavior;
it was not silently optimised away. The client's initial 600-second wait
expired, but the original remote reference kept running and was confirmed
alive. finish_production_probe.py prepared the GPU while waiting, then fetched
the completion provenance and checked the fully closed output files. No
scientific run was restarted and no tolerance was relaxed. Future adaptive
reference waits allow a larger orchestration budget only.

This evidence covers 0.3 physical seconds of phase zero, not multiple peak
periods, all four phases, energy conservation, or long-time phase drift.
The former prototype speedup numbers are not a benchmark of this controller.
