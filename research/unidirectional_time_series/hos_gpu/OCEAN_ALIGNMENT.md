# Direct HOS-Ocean MPI right-hand-side alignment — 27 September 2026

The double-precision CUDA M=5/q=3 adapter passed all 48 full-grid field
checks against an instrumented executable linked to unchanged HOS-Ocean
numerical module objects. This is direct solver-code evidence, additional
to the earlier independent MATLAB-convolution and shared CPU-backend tests.
It is not yet a time-integration or production-size validation.

## Reference authority

Upstream commit: 4deb3b4913d993c4e6ea16f736e5fc5792e14f12 (v2.1.0).
The active JONSWAP campaign records binary SHA-256
93d64252f732449fad9d2171979c91e13b0b9de1d17a73d7587fa00721efcbb6.
It matches `hos-directional-fourphase-20260926-v2/build/sources/HOS-Ocean`.

The probe copies that build's main program, inserts a call after modal-space
initialization, exports test fields and exits before time integration. All
numerical module object files are linked unchanged; only the main object is
replaced. Compiler, module and object hashes are retained in the provenance.
It runs the MPI build with ONE rank. It does not modify or relaunch existing
production jobs, install software or change any hardware settings.

An initial reference used the preceding MPI-check build. The two builds have
identical RHS, FFT, filter, variable and modal-initialization source hashes,
but their full binaries differ (production has a probe-output fix). To remove
that ambiguity, all eight references were regenerated with production-build
objects. Every resulting reference binary is byte-for-byte identical to the
already GPU-tested reference. The SHA-256 equality is recorded, so no duplicate
GPU numerical run was needed. First exporter attempt used mixed Fortran
integer widths in its header; it was superseded before GPU comparisons by
an explicitly int32 header. No values from that first export are used here.

References: `artifacts/hos_gpu/ocean-20260927T013931Z/` locally, and
`/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-rhs-probe-20260927T013931Z`
on 93. GPU snapshot `/root/hos-rhs-check-20260927T013705Z` on 92.
Compact metrics, source/object hashes and reference identities are in
ocean_alignment_results.json; the raw GPU log is under the matching local
artifacts timestamp. New diagnostic source: ocean_probe.f90 and ocean_q3.cuh.

## Important differences resolved

1. **Linear splitting.** solveHOS_lin primarily returns nonlinear terms;
   the nonzero-mode linear evolution is handled separately by the integrator.
   Its mean psi equation already includes -g*eta_mean. The probe explicitly
   removes that mean term to export a pure nonlinear reference, then restores
   G0*psi and -g*eta for the full RHS. The CUDA adapter uses the same two
   definitions. No sign, gain, phase or reference fitting is performed.
2. **Partial dealiasing.** q=3 means a 2Nx by 2Ny working grid and
   order_max=q-1=2. Called odd-degree products (3,5) are projected back to
   |mx|<=Nx/2, |my|<=Ny/2; degree-two/four products are not automatically
   projected. Eta powers, vertical-velocity orders, potentials and W products
   follow their actual call sites in resol_HOS.f90, not the earlier complete
   polynomial-support prototype.
3. **FFT normalization.** HOS-Ocean uses doubled positive-x modal amplitudes
   except the zero and base Nyquist columns. The adapter converts that policy
   to cuFFT's ordinary R2C convention through explicit expansion/reduction.
4. **Nyquist and MPI.** MPI and serial reduce_C paths are not identical.
   This adapter targets the MPI path actually used by the campaign, including
   its y-Nyquist doubling on reduction. It does not silently substitute the
   serial zero-Nyquist policy. Expanded input fields themselves are checked
   against official FFT routines before RHS comparisons.
5. **Units and zero modes.** Test dumps retain the code's nondimensional
   xlen_star, ylen_star, depth_star and g_star. Mean eta and mean psi are
   included. No mean removal, inferred dimensional rescaling or correction.

The new q3 adapter is separate from the earlier full-support HOSRHS class.
Previous speed measurements describe that earlier class and must not be
relabelled as performance measurements of the q3 adapter.

## Tests and results

Base grids 64x32 and 128x64, expanded to 128x64 and 256x128. Both use M=5,
q=3, constant depth, no current, absorption, ramping or breaking. Four inputs
per grid: low modes, near-cutoff modes, a full interior directional spectrum,
and explicit x/y Nyquist modes. The interior spectrum has respectively
961 and 3969 positive-x modes plus conjugates. These are deterministic
operator tests, not a physical evolution/convergence claim.

Each of eight cases checks all cells in expanded eta/psi, nonlinear eta/psi
RHS and full eta/psi RHS: 48 comparisons, all passed. Acceptance fixed before
testing: maximum absolute error <=2e-12 + 1e-8*reference maximum.
The previous eight GPU CTest regressions also passed.

| Field | Worst relative L2 over all eight cases |
|---|---:|
| Expanded eta | 4.37e-16 |
| Expanded psi | 3.71e-16 |
| Nonlinear eta RHS | 1.69e-11 |
| Nonlinear psi RHS | 6.61e-15 |
| Full eta RHS | 3.68e-15 |
| Full psi RHS | 1.02e-15 |

The largest relative nonlinear-eta error occurs in a very small near-cutoff
response; its maximum absolute error is about 6.05e-18. Thus checking only
the full RHS would have hidden that more sensitive comparison; both passed.

Remaining gates: production-size/actual-input checks, multi-rank arithmetic
comparison if needed, FP32 on this q3 adapter, then matching linear evolution,
Cash-Karp stages and adaptive acceptance before long-time comparisons.
No revised acceleration ratio or long-time accuracy claim is made here.
