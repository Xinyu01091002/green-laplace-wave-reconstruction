# Direct 11128-parent R4 versus MF12 comparison

Authorized 2026-09-26: compare timing and accuracy at the actual large support,
not only a small subset. No HOS evolution is part of this benchmark.

Remote root: `/home/lxy/green-laplace-unidirectional-time-series-runs/r4-mf12-jonswap11128-20260926-v1`.
Screen: `r4-mf12-jonswap11128`. Status: status.json, per-stage logs/resources.

Frozen input: JONSWAP gamma3.3, the minimal11128-mode99%-variance support;
seed20260925, same phases as design; exact retained-spectrum normalization to
kp Hs/2=.12, Hs=8.6021505376 m. kp=.0279, kph=1, domain50x20 lambda;
grid1024x512. One snapshot at t=0. No frequency/output filtering or fitting.
Strict-zero elevation/mean-potential sector is excluded in both methods.

SPARK exact source snapshot:31962e4529aeddca1f9e3ed4d566e330e23b0ea9.
Five native subharmonic source files are copied only into ignored local runtime
artifacts and remote run storage; they are not imported into public source.
Frozen manifests retain SHA256. Historical Stokes-corrected FORGE is not used.

Input to R4 is the Hermitian MATLAB FFT of real eta11, constructed directly
from the prescribed modal coefficients. Factors2 for eta20 and2*sqrt(g) for
psi20 follow current SPARK forward code/contracts, not a fitted adjustment.
MF12 order-one surface parity is checked before timed second-order work.

Cold MATLAB processes run sequentially on CPU40: prepare, R4 eta, R4 psi,
MF12 order2/difference surface, comparison. Each gets GNU time -v wall/RSS,
and the controller samples process-tree memory. No full Neumann pair diagnostic
is attached to R4. Existing low-output R4 repair and all rejection checks remain.
If either R4 operator rejects its output, preserve the rejection and do not
claim a valid field/error for it. MF12 still runs to produce the reference.

MF12 timing is split into full order2 coefficient generation, sector selection,
and one surface call producing eta20/psi20. Its unmodified coefficient API also
generates22 coefficients: report that scope honestly when comparing against
R4's20-only operator. No third-order coefficients, four-phase repetition, or
second amplitude initialization are timed. This is not the total HOS setup cost.

Resource preflight:605 GiB available, about8 CPU background load. Own-job RSS
guard160 GiB, minimum host available64 GiB, one benchmark process at a time.
No other jobs are stopped. Outputs contain dimensional full-grid L2/Linf/RMS
errors and norm ratios, with no gain/phase/offset alignment.

This file initially records the launched experiment, not its result. See the
subsequent result note and saved remote status for completion and measured data.
