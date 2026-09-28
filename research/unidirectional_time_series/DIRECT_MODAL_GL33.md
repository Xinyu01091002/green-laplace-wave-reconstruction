# Direct joint-modal Green--Laplace eta33 prototype

## Scope and graph

This prototype extends the Cartesian joint-modal time-series representation
to positive pure-sum eta33. It uses the same peak-centred observed
first-harmonic wave-group band and initial conditional direction allocation as
the eta22 trial. It does not use parent-pair/triple enumeration, MF12 third
order, a Stokes correction, fitted coefficients or a polynomial resolvent.

The inner nested-GL state is generated as

```text
E2   = (-A_Q S2 Sd + C2 Sk)/4,
Phi2 = i (C2 Sd-S2 Sk)/4,
```

where `S2,C2` are the rank-8 Gauss--Laguerre modal responses at the frozen
inner scale. The required eta22/Phi22 spatial and temporal derivatives are
reconstructed on the same joint grid. They are combined with the first-order
state to form the frozen directional cubic forcing `Fd,Fk`. The outer result
is

```text
E3 = A_Q S3 Fd-i C3 Fk,
eta33 = E3/h^2,
```

with a separate frozen outer scale and rank-8 Gauss--Laguerre response.

The archived first run kept absolute FFT-bin labels, giving 79 internal
frequency slots for bins 7--26 and 76 for bins 6--25. The maintained executor
now removes the common carrier exactly: every contiguous 20-bin input uses
`3*20-2=58` baseband slots. Thus `128` and `256` denote the number of
Cartesian modes in each horizontal direction, not the temporal sample count.
Current principal work arrays are approximately `128 x 128 x 58` and
`256 x 256 x 58` complex values. The returned record still has the original
1101 time samples.

## Frozen-graph consistency

An on-grid single-direction, single-mode fixture matched the frozen spatial
`gl_no_stokes_eta33` graph at the probe and at `t=0` to `1.10e-15` relative.
Both executors used GL8 for the inner and outer resolvents and reported zero
pair and triple loops.

## HOS focused-wave-group result

Run `hos-modal-grid-wavegroup-gl33-20260928-v1` used the completed HOS Hilbert
third phase sector at five probes. That sector is not certified exact
perturbation order three: the HOS initial state contains 11+20+22 but no
31/33, and free transient content is not removed.

The 256 x 256 results are:

| Probe | Retained bins | 128-to-256 main-window change | HOS third-sector relative L2 | q projection RMS | dispersion residual RMS | GL time (s) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 7--26 | 1.940% | 9.315% | 2.57% | 2.79% | 14.92 |
| 2 | 7--26 | 2.011% | 6.926% | 2.83% | 3.15% | 14.79 |
| 3 | 7--26 | 2.011% | 7.004% | 2.83% | 3.15% | 14.92 |
| 4 | 6--25 | 1.932% | 5.848% | 2.76% | 3.00% | 13.41 |
| 5 | 6--25 | 1.932% | 5.925% | 2.76% | 3.00% | 13.23 |

The mean HOS-sector relative L2 error is 7.004%; the maximum is 9.315% at the
centre probe. Total 256-grid compute time for five independent calls was
71.27 s. The inner invalid Laplace-domain source fraction was at most of order
`1e-16` on the 256 grid and the outer fraction was of order `1e-32`.

The result is materially closer to this HOS sector than the older
initial-spectrum-only, time-snapshot spatial GL33 baseline (about 15.8--16.1%
at the five probes). The comparison is informative but not a single-variable
benchmark: the new route also uses the observed first harmonic and a selected
wave-group band. Grid convergence is weaker than eta22, so the 256-grid
curves remain a current approximation rather than a final spatial limit.
