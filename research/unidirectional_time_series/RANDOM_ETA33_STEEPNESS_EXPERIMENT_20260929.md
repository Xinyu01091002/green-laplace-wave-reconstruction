# Unidirectional random-phase eta33 steepness experiment, 2026-09-29

## Design and completed execution

This experiment replaces the coherent focusing phases of the earlier wave-group
family by independent uniform positive-mode phases while preserving the same
unidirectional semi-Gaussian modal magnitudes. Three reproducible MT19937
realizations use seeds `20260925`, `20260926` and `20260927`. Each realization
uses the same random phases at every amplitude.

The nominal amplitude ladder is

```text
Akp = 0.02, 0.04, ..., 0.18
```

at fixed `kph=1`. No realization is renormalized to a target Hs. Consequently
the actual random-sea steepness is much smaller than the nominal focusing
label:

| Nominal Akp | Actual kp Hs / 2 |
|---:|---:|
| 0.02 | 0.003243 |
| 0.04 | 0.006486 |
| 0.06 | 0.009729 |
| 0.08 | 0.012973 |
| 0.10 | 0.016216 |
| 0.12 | 0.019459 |
| 0.14 | 0.022702 |
| 0.16 | 0.025945 |
| 0.18 | 0.029188 |

For each seed and amplitude, four global phases were initialized with MF12
`11+20+22+33`, no `31`, `muStar=0`, and linear frequencies. HOS used order 5,
50 Tp duration and Tp/40 probe output. A highest-amplitude smoke test passed,
then all 108 production HOS phase runs completed with exit code zero. At most
24 single-thread HOS processes ran concurrently. Solver SHA-256:
`13c4bae40855f3d3a29da12d5567559c4a41e415c27bc5f67d5b927db6c2f060`.

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-eta33-random-steepness-20260929-v1`.

The first and third Hilbert--four-phase sectors were extracted on the complete
50 Tp records before cropping to the fixed 10--40 Tp scoring interval. GL8
used only the non-enumerating `gl_unidirectional_time_series` implementation.
No gain, phase, time or spatial alignment was fitted.

## Direct comparison is not an eta33 validation

The raw GL8--HOS-third-sector errors are:

| Nominal Akp | Actual kp Hs/2 | Median error | Three-seed range | Median GL/HOS norm |
|---:|---:|---:|---:|---:|
| 0.02 | 0.003243 | 100.01% | 100.00--100.04% | 0.0159 |
| 0.04 | 0.006486 | 99.96% | 99.81--100.29% | 0.0634 |
| 0.06 | 0.009729 | 99.81% | 98.91--101.03% | 0.1418 |
| 0.08 | 0.012973 | 99.40% | 97.20--102.74% | 0.2494 |
| 0.10 | 0.016216 | 98.57% | 95.36--105.90% | 0.3829 |
| 0.12 | 0.019459 | 97.17% | 94.24--110.80% | 0.5364 |
| 0.14 | 0.022702 | 95.05% | 94.13--117.38% | 0.7016 |
| 0.16 | 0.025945 | 94.88% | 92.20--125.08% | 0.7545 |
| 0.18 | 0.029188 | 96.28% | 88.68--133.08% | 0.7820 |

The trend is seed-dependent and is not a monotone steepness law. At low
steepness the GL output is tiny relative to the HOS third phase sector and the
normalized waveform inner product is approximately zero. The median inner
product rises from `0.007` to `0.438` over the ladder, but remains realization
dependent.

The decisive diagnostic is amplitude scaling. Across the nine amplitudes,
the GL `eta33` norm scales with powers `2.97`, `3.00` and `3.01` for the three
seeds, as a cubic quantity should. The raw HOS third-sector norm scales only
with powers `1.06`, `1.08` and `1.54`. Thus the HOS third phase sector is not a
pure cubic `eta33` reference for these long random-wave records. Treating its
raw relative error as GL physical error would be incorrect.

## Amplitude-order regression

The nine amplitudes permit a pointwise odd-power regression

```text
q(a,t) = a q1(t) + a^3 q3(t) + a^5 q5(t) + a^7 q7(t),
a = nominal_Akp / 0.18.
```

This removes the dominant linear-in-amplitude leakage before comparing the
cubic coefficients. Four-term fit residuals are `1e-5--1e-4`, and changing to
the three-term basis changes the extracted HOS cubic coefficient by `4.2--8.8%`.
The resulting cubic comparisons are:

| Seed | Cubic relative error | Cosine | GL/HOS cubic norm |
|---:|---:|---:|---:|
| 20260925 | 200.61% | 0.287 | 2.049 |
| 20260926 | 45.63% | 0.890 | 0.882 |
| 20260927 | 86.03% | 0.564 | 0.804 |

The regression is numerically consistent but does not make the HOS cubic
coefficient equal to positive pure-sum `eta33`. It contains the complete cubic
response admitted by the nonlinear evolution, including main-harmonic and
mixed-sign sectors. In a broadband random wave these sectors overlap in
temporal frequency, whereas the current GL result contains only positive
pure-sum `eta33`. This explains the strong seed dependence and prevents the
current HOS coefficient from serving as a unique `eta33` oracle.

## Numerical bounds

The fixed `4 omega_p` GL input retains an observed first-sector projection with
relative errors `2.75--4.63%`. The final `L/h=1600` successive-domain changes
are `0.198--0.315%`, so the predeclared `5e-4` domain gate is not met. These
errors are material for a precision claim and are reported as failures of the
declared numerical gates. They are nevertheless much smaller than the raw
`89--133%` sector mismatch and do not alter the conclusion that the reference
sector is not pure `eta33`.

No input frequencies were silently removed to force the gate, and no result is
called converged. A future precision run would need a broader representable
input band and larger auxiliary domains, but should only be undertaken after a
sector-matched physical reference is defined.

## Conclusion

The experiment successfully extends the nominal amplitude scale through
`Akp=.18`, but the actual random-sea range is only
`kp Hs/2=0.00324--0.02919`. All simulations are stable. The experiment does
not validate or falsify GL `eta33` accuracy because the available long-record
HOS third phase sector is not the same physical component.

The next valid comparison must isolate positive pure-sum `eta33` itself, for
example from a perturbation-order spatial reference or a sector-complete
third-order decomposition. Merely increasing the number of random seeds or
fitting the raw third phase sector more accurately will not close that physical
definition gap.

## Evidence

Local ignored outputs:

- `artifacts/unidirectional_time_series/hos-eta33-random-steepness-20260929-v1/analysis/`
- `artifacts/unidirectional_time_series/hos-eta33-random-steepness-20260929-v1/amplitude_order/`

Tracked reproduction entries:

- `prepare_hos_random_eta33_steepness.m`
- `deploy_hos_random_eta33_steepness.py`
- `analyze_hos_random_eta33_steepness.m`
- `analyze_hos_random_eta33_amplitude_order.m`
- `plot_hos_random_eta33_results.m`
- `collect_hos_random_eta33_steepness.py`
