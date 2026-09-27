# Original HOS-Ocean linear-split Cash-Karp step

ocean_cash_karp.cuh implements the six-stage Cash-Karp 5(4) attempted step
from the actual HOS-Ocean runge_kutta.f90, using the q=3 nonlinear RHS.
It does not yet implement the outer adaptive step-size controller.

State uses HOS-Ocean's normalised/doubled-positive-x modal amplitudes. The
linear dispersion rotation acts on (g/omega*eta, psi), with g/omega=1 and
omega=0 at the spatial zero mode. The nonlinear psi RHS includes the mean
gravity term just as solveHOS_lin does. Every stage is evaluated on GPU;
the state, stage slopes and candidate remain on device. A preallocated
two-stage reduction returns only one error-estimate scalar to the host.

The stage/final rotation and embedded error calculation are ported separately.
In particular, the error routine's original time-dependent rotation
omega*(t0+c_j*h) is preserved, rather than substituted with another integrating-
factor convention. Tables are explicitly zero-initialised before assigning
the original nonzero Cash-Karp coefficients. No error formula is fitted.

attempt() leaves the starting state intact; accept() explicitly commits the
candidate. This separates later acceptance/rejection logic from the attempted
step and permits rejected attempts without host state copies. No automatic
tolerance relaxation or hardware setting changes are made.

## Independent official reference

ocean_ck_probe.f90 calls the unchanged production module's
RK_adapt_2var_3D_in_mo_lin using its fill_butcher_array routine and exports
the initial/candidate modes and original embedded error estimate. It runs
fixed attempted steps, without a production time loop. The same four low,
near-cutoff, broadband and Nyquist input cases as OCEAN_ALIGNMENT.md are
used on two grids, with three trials each:

- t0=0, h=.01;
- t0=0, h=.2;
- t0=.37, h=.2.

All times here are nondimensional. The nonzero-t0 trial checks the error
rotation. This totals 24 attempted steps. GPU checks are FP64 at this gate;
successful FP32 spatial-RHS checks do not certify FP32 temporal error estimates.

Predeclared state gate: max modal absolute difference <=2e-12 + 1e-8*reference
maximum for eta and psi. Error-estimate gate: absolute difference
`<=5e-15 + 1e-5*abs(reference estimate)`. Also require the same accept/reject classification
at the existing absolute threshold 1e-12. That classification test is not a
test of the outer controller's next-step proposal or a complete adaptive run.

## Verified results — 27 September 2026

All 24 attempted-step comparisons passed, with the same 7 accept and 17 reject
classifications as the original solver at 1e-12. Accepted candidates' maximum
eta/psi relative L2 difference was 1.33e-16. Across all candidates, including
rejected large-step candidates, the maximum state relative L2 difference was
7.60e-14. The largest relative discrepancy in the embedded error was about
1.08e-7, occurring for a very small estimate; the predeclared mixed absolute/
relative estimate gate passed throughout.

Two deliberately large h=.2 attempts on the 128x64 broad-spectrum case produced
explosive candidates with error estimates around 6e118 in BOTH implementations.
They were classified as REJECT, never accepted or counted as successful wave
evolution. Their very large absolute estimator difference is not hidden by
a small-number summary; complete values are retained in the JSON record.

Eight rejected trial-2 attempts were followed, in the same GPU solver object
without re-uploading initial state, by the smaller trial-1 attempt. All matched
the original small-step reference. This includes recovery after the explosive
candidate and checks that attempt() does not contaminate the retained state.
One smaller broad-spectrum attempt remains correctly classified as rejected;
matching a reference is not equivalent to accepting its time step.

Reference: gpu-rhs-probe-20260927T022301Z on 93. Final GPU snapshot:
/root/hos-rhs-check-20260927T023327Z on 92. Full records and source hashes:
ocean_cash_karp_results.json. The earlier 022936Z run passed the original 24
checks; the final run adds rejection/retry checks. A preceding upload timed
out before compilation; it supplies no numerical evidence. Subsequent uploads
use a compressed, hash-verified snapshot. Eight previous CTest regressions pass.

This gate checks FP64 attempted steps on small grids, not FP32 temporal error
control or an actual JONSWAP trajectory. The production-sized JONSWAP RHS
precision check is separately recorded in PRODUCTION_RHS.md. Next steps are
the outer adaptive controller and bounded trajectory comparisons, before any
replacement of existing production runs or revised speedup claim.
