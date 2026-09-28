# Direct joint-modal Green--Laplace eta20 prototype

## Scope

This prototype reconstructs the nonzero-spatial-wavenumber second-order
difference sector from the same two authorized inputs as eta22/eta33: the
observed first-harmonic probe record supplies temporal complex coefficients,
and the initial spectrum supplies conditional directions. Strict spatial
`Q=0` is excluded as a separate mean/volume/gauge sector. Temporal zero with
nonzero directional difference `Q` is retained internally.

The comparison projects both HOS and GL identically to

```text
0 < abs(omega) < 0.5 omega_p.
```

Thus the displayed result excludes temporal DC. The HOS reference is the
four-phase even-sector average and can contain higher even harmonics; it is not
certified pure perturbation-order eta20.

## Modal graph

For the analytic first-order joint spectrum, the difference sources are built
by products of filtered fields with conjugated filtered fields. Their joint
spectra are functions of output wavevector difference `Q` and signed temporal
difference `sigma`. The shared-scale GL response is evaluated as

```text
ep_j = w_j/delta_q exp[x_j-(a_Q-sigma) x_j/delta_q],
em_j = w_j/delta_q exp[x_j-(a_Q+sigma) x_j/delta_q],

Md = sum_j -a_Q (ep_j+em_j)/4,
Mk = sum_j  i   (ep_j-em_j)/4,

eta20_hat = (Md Fd_hat+Mk Fk_hat)/(2h).
```

The factor two relative to the low-level frozen spatial helper is the declared
physical assembly `2*candidate_half`, not a fitted gain. A two-mode on-grid
fixture matches the frozen shared-GL6 spatial graph to `3.79e-16` relative.
Neither executor uses a parent-pair loop.

For 20 contiguous input frequencies, exact baseband difference convolution
uses `2*20-1=39` internal temporal slots.

## HOS focused-wave-group result

Run `hos-modal-grid-wavegroup-gl20-20260928-v2` uses GL16 and 20 input
frequencies. The 512 x 512 results are:

| Probe | 256-to-512 main-window change | HOS even-sector subharmonic L2 | Candidate/reference norm | GL time (s) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.81% | 5.30% | 0.95 | 10.55 |
| 2 | 0.88% | 2.44% | 0.99 | 10.62 |
| 3 | 0.88% | 2.38% | 0.99 | 10.42 |
| 4 | 1.59% | 2.50% | 1.01 | 10.44 |
| 5 | 1.59% | 2.39% | 1.01 | 10.70 |

The waveforms reproduce the broad negative set-down pulse without fitted
alignment. The centre discrepancy is mainly an amplitude difference near the
minimum; the four lateral probes are within about 2.4--2.5% in the selected
main window.

## Band, grid and rank checks

Expanding the first-harmonic input from 20 to 36 frequencies changes the final
main-window eta20 by 0.22% at probes 1 and 4. Twenty input frequencies are
therefore sufficient for this focused wave group.

The 256-to-512 change is 0.81--1.59%, so eta20 is less spatially converged than
eta22. The 512 grid is retained for the displayed final curves.

At 256 x 256, changes relative to GL16 are:

| Probe group | GL6 change | GL12 change |
| --- | ---: | ---: |
| Centre | 13.42% | 2.67% |
| Near-centre pair | 9.01% | 1.80% |
| Outer pair | 8.97% | 2.01% |

GL6 is not adequate. GL12 is closer but not rank-converged at a one-percent
target; GL16 is the retained result. At the outer probes GL12 happens to have
a smaller HOS residual than GL16, which is treated as numerical/physical error
cancellation rather than a reason to select GL12.

## Strict-zero audit

The source energy at strict spatial `Q=0` grows from about 5% at the centre to
about 20% at the outer probes as the grid is refined. It is intentionally
excluded and is not counted as a failed nonzero resolvent. On the 512 grid,
the genuinely nonzero-`Q` source outside `a_Q>abs(sigma)` is only
`4.6e-6--1.1e-4` of total source energy. The v1 and corrected-audit v2 runs
have identical predictions and reported waveform metrics.

This result supports nonzero eta20 reconstruction for the present wave group.
It does not solve the strict mean, mass/volume constraint, return flow or every
signed low-frequency sector.
