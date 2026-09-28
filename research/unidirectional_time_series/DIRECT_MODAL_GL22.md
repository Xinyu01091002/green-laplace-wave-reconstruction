# Direct joint-modal Green--Laplace eta22 prototype

## Scope

This prototype evaluates the positive pure-sum second-order elevation at a
fixed probe. It uses both the observed first-harmonic temporal coefficients
and the conditional directions supplied by the initial directional spectrum.
It does not evaluate eta20, mean flow, signed sectors or third order.

The only approximation of the modal resolvent is the prescribed
Gauss--Laguerre rank. No Chebyshev, fitted kernel or parent-pair evaluator is
used. The first prototype projects non-Cartesian directional wavevectors onto
a Cartesian modal grid. This projection is a separate numerical approximation
and is audited independently from the Green--Laplace rank.

## Joint first-order representation

At a probe chosen as the spatial origin, write the analytic first-order field

```text
v(x,t) = sum_alpha C_alpha exp(i q_alpha.x-i nu_alpha t),
q_alpha = h k_alpha,   nu_alpha = omega_alpha sqrt(h/g).
```

For fixed depth, `nu^2=q tanh(q)` is monotone in the magnitude `q`. Temporal
frequency therefore determines the free-wave magnitude, but not direction.
The initial spectrum supplies the conditional direction distribution, while
the observed first harmonic supplies the temporal complex coefficient.

The directional points normally do not lie on a Cartesian `(qx,qy)` lattice.
The prototype deposits each point onto its four neighbouring Cartesian modes
with bilinear weights whose sum is one. Consequently the spatial-origin
coefficient at every temporal bin is conserved exactly. The deposited modes
do not satisfy the free-wave dispersion relation exactly; the weighted
wavevector error and dispersion residual are returned in the audit.

## Modal Green--Laplace resolvent

Let

```text
A_Q = Q tanh(Q),   a_Q = sqrt(A_Q),   s = nu_1+nu_2.
```

The exact positive-sum response is

```text
eta22_hat(Q,s) = [-A_Q Sd_hat(Q,s)+s Sk_hat(Q,s)]
                 / [4 h (s^2-A_Q)].
```

Its Green--Laplace representation is

```text
1/(s^2-a_Q^2) = integral_0^inf exp(-s tau)
                 sinh(a_Q tau)/a_Q dtau,

s/(s^2-a_Q^2) = integral_0^inf exp(-s tau)
                 cosh(a_Q tau) dtau.
```

For Gauss--Laguerre nodes `x_j`, weights `w_j`, scale `lambda`, and
`tau_j=x_j/lambda`, define `beta_j=w_j exp(x_j)/lambda`. The finite-rank
modal transfers are

```text
Gd_J(Q,s) = -A_Q sum_j beta_j exp(-s tau_j)
                    sinh(a_Q tau_j)/a_Q,

Gk_J(Q,s) = sum_j beta_j exp(-s tau_j) cosh(a_Q tau_j).
```

They are evaluated in the balanced form

```text
ep_j = w_j/lambda exp[x_j-(s-a_Q) tau_j],
em_j = w_j/lambda exp[x_j-(s+a_Q) tau_j],

Gd_J = sum_j -a_Q (ep_j-em_j)/2,
Gk_J = sum_j       (ep_j+em_j)/2.
```

## Moving the node damping to the output frequency

For every positive-sum contribution,

```text
exp(-tau_j nu_1) exp(-tau_j nu_2)
  = exp[-tau_j (nu_1+nu_2)]
  = exp(-tau_j s).
```

Therefore the eight undamped first-order filtered fields and the two source
fields need to be formed only once. After transforming the sources to joint
output coordinates `(Qx,Qy,s)`, the Gauss--Laguerre node damping is applied as
an output-frequency multiplier. This is algebraically identical to damping
each parent at every node, but removes the repeated field transforms from the
rank loop.

The source fields are

```text
Sd = 2 v v_nu2 + v_nu^2 - hx^2 - hy^2,
Sk = 2 v v_q2_over_nu + 2 (hx jx + hy jy).
```

The joint convolution is evaluated by transforms and pointwise products. No
parent pairs or triples are enumerated. The final modal spectrum is contracted
directly to the probe before the full-record temporal synthesis.

## Declared numerical errors

The prototype reports separately:

1. observed temporal-bin reconstruction error at the probe;
2. Cartesian wavevector projection error;
3. free-wave dispersion residual after wavevector projection;
4. source energy falling outside the convergent Laplace domain `s>a_Q`;
5. GL rank and modal-grid dimensions.

No anti-alias padding, amplitude correction, phase correction or HOS-fitted
coefficient is introduced. The current modal grid is a feasibility device,
not yet a production low-cost solver.

## Initial synthetic smoke result

A non-physical timing fixture used 36 temporal bins, 13 directions, 468 joint
modes, 1101 output samples, and GL8. Its purpose was to exercise the complete
multi-mode array path; it is not HOS validation. The observed coefficient at
the probe was conserved to about `1e-16` at every tested grid.

| Cartesian grid | Time (s) | Weighted q error | Weighted dispersion residual | Source energy outside s>a | Output norm |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 32 x 32 | 0.131 | 37.66% | 108.9% | 7.8851% | 0.009554 |
| 64 x 64 | 0.177 | 18.61% | 40.95% | 1.5430% | 0.008684 |
| 128 x 128 | 0.567 | 9.21% | 16.39% | 0.31533% | 0.008135 |
| 256 x 256 | 2.611 | 4.58% | 7.41% | 0.044078% | 0.0080 |

The 128-to-256 output change was 11.09% in relative L2. Uniform Cartesian
gridding is therefore not converged for this deliberately broad
`kh approximately 0.05--3.86` stress test. This does not determine its
accuracy for the energetic range of a focused wave group.

## HOS focused-wave-group result

Run `hos-modal-grid-wavegroup-gl22-20260928-v1` applied the same GL8 executor
to the five saved HOS probes. Each input used a peak-centred contiguous band,
at least 20 frequencies, and at least 99.9% of the complete first-harmonic
record energy. The centre and near-centre probes selected bins 7--26; the two
farther probes selected bins 6--25. Actual retained fractions were
99.9999836%--99.9999977%, with first-record relative projection errors
0.01524%--0.04055%.

The 256 x 256 results are:

| Probe | 128-to-256 main-window change | HOS second-sector relative L2 | q projection RMS | dispersion residual RMS | GL time (s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.5394% | 4.1507% | 1.71% | 1.84% | 4.478 |
| 2 | 0.5618% | 2.9721% | 1.90% | 2.09% | 4.387 |
| 3 | 0.5618% | 3.0126% | 1.90% | 2.09% | 4.293 |
| 4 | 0.4366% | 2.4512% | 1.78% | 1.97% | 3.506 |
| 5 | 0.4366% | 2.4912% | 1.78% | 1.97% | 3.482 |

The energy-weighted source fraction outside the convergent Laplace domain was
at most `7.90e-17`; the deposited station coefficients reproduced the selected
first-harmonic bins to `2.24e-16` relative or better. The coarser 64 x 64
predictions happened to have smaller HOS errors, but the HOS error increased
under grid refinement while the grid-to-grid change decreased. That behaviour
is treated as numerical/physical error cancellation, not as evidence for
selecting the coarser grid.

These results validate the feasibility of the Cartesian approximation for
the energetic focused-wave-group range. They do not validate the discarded
spectral extremes or establish that the four-phase HOS sector is an exact
perturbation-order reference. A NUFFT remains an optional refinement rather
than a prerequisite for continuing the wave-group calculation.
