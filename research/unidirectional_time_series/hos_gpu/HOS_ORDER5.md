# Constant-depth HOS right-hand side through degree five

This implements Taylor HOS for a periodic, flat-bed free surface in eta and
true surface potential psi. It is a GPU right-hand-side evaluator, not yet
a time integrator or a drop-in replacement for the active CPU campaign.
Supported truncation degree is M=1..5. FP32 and FP64 share the algorithm.

## Recursion

Coordinates and Fourier sign are as in SECOND_ORDER.md. Let D_j be the
vertical derivative at still water of the flat-bed harmonic extension:

    D_j(k) = |k|^j                    for even j,
             |k|^j tanh(|k| h)       for odd j.

At k=0 all derivatives j>=1 are zero. Denote the degree-n still-water
potential by phi_n and the degree-n surface vertical velocity by W_n:

    phi_1 = psi
    phi_n = -sum_(j=1..n-1) eta^j/j! D_j phi_(n-j)
    W_n   =  sum_(j=0..n-1) eta^j/j! D_(j+1) phi_(n-j).

Each occurrence of eta or psi counts as degree one. These are homogeneous
components of the instantaneous nonlinear operator, not harmonic phase
sectors and not free/bound components inferred from an observed record.

Write Q=|grad eta|^2. Degree-n boundary-condition terms are:

    E_n = W_n - delta_(n,2) grad eta . grad psi + Q W_(n-2)
    P_n = -delta_(n,1) g eta - delta_(n,2) |grad psi|^2/2
          + 1/2 sum_(p+q=n) W_p W_q
          + Q/2 sum_(p+q=n-2) W_p W_q.

Only positive W indices are included. eta_t=sum E_n and psi_t=sum P_n,
through n=M. Required Taylor terms for true surface potential are retained;
there is no empirical Stokes adjustment, amplitude fitting or mean removal.

## Spectral contract and implementation

Input fields must have declared rectangular support (bx,by). Degree-n
products are projected to (n*bx,n*by) before their spectrum is reused or
their contribution is returned. All degree-M outputs must fit strictly below
working-grid Nyquist: 2*M*bx<Nx and 2*M*by<Ny (Ny=1 requires by=0).
Under that contract, projections remove numerical residue but no physical
polynomial output. Input eta and psi are also projected to their declared
band. Do not declare a narrow band for a broader observed input.

The production HOS-Ocean M=5/q=3 calculation uses a different retained-band
and intermediate-dealiasing policy. This full-support validation path does
NOT reproduce that policy or accept arbitrary production spectra merely
because Nx/Ny match. A port of those masks and direct CPU HOS comparisons
remain necessary before production replacement.

All phi spectra, eta powers, W components, gradients, work arrays and RHS
outputs are allocated once on the GPU. Ordinary evaluate() has no host
transfer or device allocation; optional observers copy fields for validation.
Repeated evaluations preserve uploaded state. No hardware settings are changed.

## Independent reference and checks

`make_hos_reference.m` uses MATLAB direct discrete Fourier convolution
(`conv2` on coefficient arrays), not FFTs or CUDA kernels. Its coefficient
grid covers the complete degree-five support. All requested product degrees
fit, so the `same` crop removes no physical term. It checks conjugate symmetry,
the linear limit, and the independent second-order DNO identity before export.

Three finite depths (0.15,1.3,20 m), the existing two non-collinear conjugate
mode pairs and nonzero mean eta/psi are used. Reference binary files contain
normalised complex Fourier coefficients for eta, psi and each W_n,E_n,P_n,
with a version-specific fixed header described by make_hos_reference.m.

`hos_order_check` compares each field separately at every degree 1..5:

- Small 64x32 grid: all 2048 physical points.
- Large 2048x1024 grid: 257 fixed deterministic points for reference error,
  plus a finite-value check over every output cell. Large-grid errors are
  **sampled errors**, not full-grid maximum error certificates.
- Per-field acceptance: max absolute error <= atol + rtol*reference max.
  FP32 atol=1e-10, rtol=5e-4; FP64 atol=1e-12, rtol=2e-9.
- Total RHS uses the same rtol and 5*atol. Tolerances are fixed before runs.
- Callback-free repeated evaluation must reproduce the same outputs exactly.

Run MATLAB with a fresh output directory, then:

    python run_p40_checks.py --orders <reference-directory>

The runner transfers and verifies source/reference hashes, uses a fresh
remote snapshot, builds for P40 sm_61, runs earlier regression tests, then
checks all depths in both precisions and a large-grid h=1.3 case in both.
No physical trajectory, long-time precision or speedup is certified here.

## Verified checkpoint

Snapshot `/root/hos-rhs-check-20260927T005438Z`, 27 September 2026. Source and
reference SHA-256 values plus all 120 per-degree field checks are retained
in order5_results.json. The local raw log is under the matching timestamp
in artifacts/hos_gpu. All checks passed, as did the eight earlier CTest
regressions and the intermediate-grid second-order checks.

Six small-grid cases cover all three depths in both precisions, with every
grid point checked. Two large-grid cases cover h=1.3 in both precisions at
the predeclared 257 points; every output cell was checked for finiteness.

Large-grid FP32: total eta_t/psi_t sampled maximum absolute errors
2.69263e-7 / 6.44133e-7. Fifth-degree eta/psi sampled relative L2 errors
1.77931e-5 / 5.09635e-6. Large-grid FP64 totals: 3.05311e-16 / 1.08247e-15;
fifth-degree eta/psi sampled relative L2: 2.63605e-14 / 3.06016e-15.

At 2048x1024, declared allocations (arrays plus FFT workspace) are
276,889,600 bytes FP32 (264.063 MiB) and 553,795,584 bytes FP64 (528.141 MiB).
These were obtained from the actual constructed evaluator, excluding CUDA
context/library overhead and any subsequent RK stage buffers. No speed claim.
