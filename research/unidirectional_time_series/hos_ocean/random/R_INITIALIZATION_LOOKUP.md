# Correction: non-pairwise R initialization implementations exist

Inspected saved source on 2026-09-26, read-only. The previous statement based
only on this repository's diagnostics/eta20 ordered-pair evaluator was too
broad. It does not describe all maintained research implementations.

Current SPARK source at C:/Users/spet5947/Documents/ChatGPT/SPARK includes:

- src/forward/subharmonic/spark_eta20_r_series.m: native fixed-FFT R2/R4/R6.
- src/forward/subharmonic/spark_psi20_r_series.m: matching true-surface Psi20.
- src/forward/subharmonic/spark_r4_low_mode_pair_repair.m: targeted q<=.1
  output convolution for numerical consistency of the same R4 expansion.
- src/forward/subharmonic/spark_difference_audit.m: no whole-support pair
  scan; reports main_operator_fixed_fft and targeted repair work separately.

Cost described by the current code is C_R Ng log(Ng)+O(Nlow Nactive), not a
mandatory O(Nactive^2) enumeration. Fixed grid/rank determines FFT graph size.
R4 still has local pair contributions at selected low output modes; do not
claim the complete R4 path has zero pair work or that repair removes model
truncation error. Exact Neumann pair diagnostics live in validation only.

Historical FORGE under Archive/GL-cleanup-20260925-025101/workspaces also
contains finite_depth_eta20_directional_neumann_resolvent.m (pair_loops=0),
but that historical route has Stokes-correction terms. It was not imported.
SPARK native eta20 declares pure Neumann forcing and has no such correction
in the inspected field operator. Current SPARK scientific/API conventions
must be respected before any adapter is implemented.

Thus R for eta20/psi20 plus GL for eta22/psi22 is a real initialization
candidate, rather than a missing algorithm. It is not yet timed or validated
on the proposed 11128-parent JONSWAP input. A bounded fixed-input benchmark
should measure runtime/RSS and consistency, compare against independent MF12
on an affordable support, and retain all initialization provenance if later
using the hybrid to assess GL against HOS. No code was copied from SPARK,
no new initializer was implemented and no new propagation was launched by
this source inspection.
