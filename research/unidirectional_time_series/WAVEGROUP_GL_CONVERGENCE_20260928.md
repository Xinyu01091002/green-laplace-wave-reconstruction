# Focused-wave-group GL convergence, 2026-09-28

## Scope

This check isolates frequency-band, Cartesian-grid and Gauss--Laguerre-rank
changes for the focused HOS wave group. Only the centre probe 1 and outer
probe 4 are used. HOS errors are reported after every run but do not select
the numerical settings.

The executors use exact baseband relocation. For 20 contiguous input bins,
the eta22 and eta33 work-frequency counts fall from the old absolute-bin
counts 53/79 to 39/58. Three-frequency complex fixtures reproduce the old
paths to `2.10e-16` for eta22 and `9.23e-16` for eta33.

## Frequency-band convergence

All band runs use a 256 x 256 Cartesian grid, GL8 and the same spatial modal
limits set by the 36-bin input. Changes are measured against 36 bins in the
main wave-group window.

| Probe | Input count | Retained bins | eta22 change | eta33 change |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 20 | 7--26 | 0.119% | 0.518% |
| 1 | 24 | 6--29 | 0.100% | 0.505% |
| 1 | 28 | 2--29 | 0.158% | 0.306% |
| 4 | 20 | 6--25 | 0.091% | 0.337% |
| 4 | 24 | 5--28 | 0.098% | 0.350% |
| 4 | 28 | 2--29 | 0.101% | 0.157% |

For this focused wave group, 20 frequencies are sufficient at the present
accuracy target. This is an observed final-output result, not a universal
fixed frequency count.

## Cartesian-grid convergence

These runs use the 20-bin wave-group input, GL8 and identical physical modal
limits at both resolutions. Changes are measured against 512 x 512.

| Probe | Grid | eta22 change | eta33 change | eta22 time (s) | eta33 time (s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 256 x 256 | 0.109% | 0.482% | 2.85 | 10.97 |
| 1 | 512 x 512 | 0 | 0 | 12.06 | 42.01 |
| 4 | 256 x 256 | 0.115% | 0.443% | 2.11 | 9.09 |
| 4 | 512 x 512 | 0 | 0 | 9.42 | 41.73 |

The HOS residual increases under refinement: at probe 1 eta22 changes from
4.15% to 4.25% and eta33 from 9.31% to 9.73%; at probe 4 eta22 changes from
2.45% to 2.55% and eta33 from 5.85% to 6.23%. The coarser residual therefore
contains modest numerical/physical error cancellation. Numerical settings are
selected by grid-to-grid change, not by the smaller HOS residual.

## Green--Laplace rank convergence

These runs use the 20-bin input and 256 x 256 grid. Changes are measured
against GL10.

| Probe | Rank | eta22 change | eta33 change |
| ---: | ---: | ---: | ---: |
| 1 | GL6 | 0.142% | 2.474% |
| 1 | GL8 | 0.021% | 0.560% |
| 4 | GL6 | 0.166% | 2.345% |
| 4 | GL8 | 0.026% | 0.527% |

GL8 is resolved for eta22 at this scale. GL8 is within about 0.6% of GL10
for eta33 and is adequate for a one-percent working target; GL10 is the more
conservative setting for a final high-accuracy eta33 curve.

## Current numerical choice

For continued focused-wave-group development:

```text
input frequencies: 20, selected adaptively from the observed first harmonic
eta22: 256 x 256, GL8, baseband work count 39
eta33: 256 x 256, GL8 for routine work or GL10 for final curves,
       baseband work count 58
```

The individual band/grid/rank changes are diagnostics, not independent error
bounds and should not be added as if statistically independent. The remaining
HOS discrepancy includes physical-model differences and the fact that the HOS
third phase sector is not certified pure perturbation order three.
