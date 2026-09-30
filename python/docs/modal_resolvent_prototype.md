# One-dimensional modal-resolvent prototype

The prototype evaluates single-direction `eta22(t)` and `eta33(t)` directly
on a uniform dimensionless-wavenumber grid. It uses baseband temporal slots,
FFT convolution of filtered modal fields, and the finite-node Green--Laplace
resolvent. It does not create an auxiliary spatial domain or enumerate input
wave pairs or triples.

## On-grid identity

A two-mode case was constructed with

```text
q2 = 2 q1,
omega2 / omega1 = 3 / 2,
```

so both modes lie exactly on the temporal DFT bins and the `q`-modal grid.
At the probe and `t=0`, the modal result agrees with the spatial FFT--GL graph:

| Component | Relative difference |
|---|---:|
| `eta11` | `3.40e-16` |
| `eta22` | `5.63e-16` |
| `eta33` | `1.97e-15` |

For 256 modal points and GL4, 25 repeated calls gave a median time of
`0.00290 s` on the local Python 3.12/NumPy environment.

## Off-grid refinement

The fixed-phase four-component fixture has 384 output samples and uses GL6.
Changes below are measured against the 65536-point modal result.

| Modal points | Single-call time | q projection RMS | `eta22` change | `eta33` change |
|---:|---:|---:|---:|---:|
| 16384 | 0.675 s | `5.19e-4` | 0.221% | 6.688% |
| 32768 | 1.444 s | `1.31e-4` | 0.023% | 0.083% |
| 65536 | 3.136 s | `7.29e-5` | reference | reference |

The 32768-point result meets a 0.1% successive-grid target for both components
on this fixture. These are local single-call measurements, not general timing
claims.

## MATLAB modal parity

For the same 512-point modal discretisation, the Python and MATLAB direct
modal executors agree to `3.10e-16` for `eta22` and `1.06e-15` for `eta33`.

