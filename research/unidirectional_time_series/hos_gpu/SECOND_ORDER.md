# Finite-depth second-order GPU right-hand side

This is an incremental HOS implementation gate, not a replacement for the
current M=5 HOS-Ocean production runs. It evaluates the complete terms through
degree two in the surface elevation eta and true surface potential psi.
No time integration is implemented at this gate.

## Equations and conventions

Horizontal domain [0,Lx)x[0,Ly), periodic; z=0 at still water, bed z=-h.
Use exp(i k.x) spatial modes, g>0, h>0, physical units (metres, seconds).
psi is phi evaluated on z=eta, not the potential at z=0. Define
D(k)=|k|tanh(|k|h) and G0=F^-1 D F. The order-two expansion is

    eta_t = G0 psi - G0(eta G0 psi) - div(eta grad psi)
    psi_t = -g eta - (grad psi).(grad psi)/2 + (G0 psi)^2/2.

This follows the finite-depth Dirichlet-to-Neumann expansion at first order
in eta and the dynamic free-surface boundary condition through degree two.
Mean elevation and the spatially uniform potential are retained, including
the mean dynamic RHS; no zero-mode correction or amplitude fitting is used.

The independent sparse convolution reference in rhs_check.cu evaluates:

    eta_t[K] = D(K) psi[K]
             + sum_(p+q=K) (K.q - D(K)D(q)) eta[p] psi[q]
    psi_t[K] = -g eta[K]
             + 1/2 sum_(p+q=K) (p.q + D(p)D(q)) psi[p] psi[q].

It sums Fourier coefficients directly at each physical point, with no FFT
or GPU kernel reuse. Complex conjugate signed modes make the input real.
Two non-collinear modes have different complex eta/psi coefficients; their
sums, differences, self-interactions, zero modes, signs and finite-depth
multipliers are therefore exercised. Depths 0.15, 1.3 and 20 m at fixed
Lx=9, Ly=7 m cover distinct depth regimes. The 1D variant uses Ny=1.

## Support and memory contract

The caller supplies fields already on the working/expanded grid. Declared
input support |mx|<=bx, |my|<=by must obey 4*bx<Nx and 4*by<Ny, with by=0
for Ny=1. These strict inequalities keep every quadratic output below the
Nyquist frequencies. The constructor checks the declared band/grid relation;
the default path does NOT inspect or filter input arrays to enforce the declared support.
An arbitrary sampled waveform cannot be passed under a false band declaration.
No aliasing-free claim is made outside this contract.

The optional diagnostic switches `evaluate(project_input, project_products)`
now enforce the declared input rectangle before differentiation and the
twice-wide quadratic rectangle before differentiating products, respectively.
Defaults remain false/false to preserve the original failing baseline.
For this declared band-limited problem these projections remove no physical
input or quadratic output modes. They are not a general licence to discard
nonzero modes in real input data. The dynamic RHS pointwise squares are not
post-filtered in this diagnostic.

This gate computes all representable quadratic outputs, without reducing
to the production base-grid retained band. It is NOT a claim of reproducing
HOS-Ocean's M=5/q=3 intermediate product masks. Higher orders require explicit
matching of those masks and CPU HOS comparisons before production use.

eta, psi, G0psi, horizontal derivatives, RHS outputs, spectrum and reusable
FFT storage stay on the GPU during evaluate(). Upload/download are explicit
boundary operations. An additional reusable eta-projection buffer is allocated
at construction; original eta is preserved across repeated calls. The reported allocated bytes cover the declared arrays
and cuFFT workspace, excluding CUDA context and library overhead. No plan or
device allocation occurs inside evaluate(). With diagnostic callbacks absent,
no host transfers occur there. Explicit diagnostic callbacks download snapshots
to measure term errors and spectral residuals; their execution is not timed as
a production benchmark.

Tests repeat evaluation on the same state to detect destructive spectral
transforms or accidental accumulation. Tolerances fixed before testing:
FP32 maximum absolute error <2e-5 and FP64 <2e-11 for each RHS field.
These tolerances concern instantaneous RHS arithmetic, not wave evolution.
