# Single-direction eta33 steepness experiment, 2026-09-29

## Question and fixed design

This experiment tests whether the discrepancy between the non-enumerating
Green--Laplace `eta33` time series and the HOS third phase sector increases
with steepness. It does not change water depth, spectrum shape, phase
realization, focusing position, record duration, harmonic extraction or the
comparison window.

The six cases use

```text
Akp = 0.02, 0.04, 0.06, 0.08, 0.10, 0.12
kph = 1
```

with the existing Alpha=1 semi-Gaussian unidirectional spectrum, four global
phases, 50 Tp HOS integration and Tp/40 output. The Hilbert--four-phase first
and third sectors are extracted on the complete 50 Tp record before cropping
to the unchanged 32--48 Tp interval. No gain, phase, time or spatial alignment
is fitted. The HOS phase sector is not labelled a certified perturbation-order
`eta33` field.

The existing Akp=.02/.12 records are reused read-only. Sixteen new HOS runs
provide the four intermediate steepnesses. Every run completed with exit code
zero. The solver SHA-256 is
`13c4bae40855f3d3a29da12d5567559c4a41e415c27bc5f67d5b927db6c2f060`.
The isolated remote root is
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-eta33-steepness-20260929-v1`.

GL postprocessing calls only `gl_unidirectional_time_series`, with the full
declared nonzero input band through `4 omega_p`, no parent-pair/triple
enumeration, no Stokes correction, GL8, and a fixed auxiliary-domain ladder
through `L/h=1600` when required. All cases pass the declared successive-domain
change tolerance `5e-4`.

## Six-point result

Raw, unaligned main-window relative L2 errors are:

| Akp | GL8--HOS error | Full-window error | GL/HOS main norm | Final L/h | Last domain change |
|---:|---:|---:|---:|---:|---:|
| 0.02 | 2.0585% | 2.1177% | 0.9812 | 800 | 4.9346e-4 |
| 0.04 | 2.7572% | 2.7962% | 0.9740 | 800 | 4.9757e-4 |
| 0.06 | 3.9454% | 3.9725% | 0.9619 | 1600 | 1.8862e-4 |
| 0.08 | 5.6037% | 5.6229% | 0.9453 | 1600 | 1.9091e-4 |
| 0.10 | 7.7106% | 7.7248% | 0.9244 | 1600 | 1.9402e-4 |
| 0.12 | 10.2390% | 10.2499% | 0.8994 | 1600 | 1.9807e-4 |

The error is strictly monotone. A pure all-point power fit gives approximately
`error proportional to (Akp)^0.896`, but this hides a nonzero low-steepness
floor. The two-term model

```text
E(Akp) = E0 + c (Akp)^2
```

fits all six points as

```text
E0 = 0.018344815
c  = 5.85366443
R2 = 0.999960300
```

where `E` is a relative fraction. Thus the result is consistent with a roughly
`1.83%` steepness-independent baseline plus an additional relative discrepancy
that grows quadratically with steepness. Because `kph=1` is fixed, this trend
cannot be attributed to changing water depth.

The main-window HOS third-sector norm scales as `(Akp)^2.9899`; the GL `eta33`
norm scales as `(Akp)^2.9448`. Both remain close to cubic overall, while the
GL/HOS norm ratio decreases monotonically from `0.9812` to `0.8994`. This is
consistent with an increasingly important higher-order/model difference, but
does not by itself identify a unique omitted physical sector.

## Fixed-input GL rank diagnostic

The same six inputs and HOS records were evaluated at GL4, GL6 and GL8. No rank
was selected against HOS.

| Akp | GL4 error | GL6 error | GL8 error |
|---:|---:|---:|---:|
| 0.02 | 14.6555% | 5.2308% | 2.0585% |
| 0.04 | 15.2493% | 5.9014% | 2.7572% |
| 0.06 | 16.2367% | 7.0255% | 3.9454% |
| 0.08 | 17.6120% | 8.5977% | 5.6037% |
| 0.10 | 19.3662% | 10.6026% | 7.7106% |
| 0.12 | 21.4873% | 13.0191% | 10.2390% |

GL rank therefore matters materially. GL6 to GL8 reduces the HOS discrepancy
by `2.78--3.17` percentage points, and the current API has no GL10/12 result.
GL8 must not be described as rank-converged. The rank effect provides a
credible explanation for much of the approximately constant low-steepness
floor, whereas the nearly quadratic growth with `Akp` remains after separating
that qualitative baseline.

Historical analytic-resolvent results at the two endpoints had smaller errors
(`0.317%` and `8.729%`), consistent with this interpretation, but no prohibited
interaction enumerator or analytic interaction reference was run in the new
experiment.

## Interpretation

The new evidence supports three bounded conclusions:

1. For this fixed unidirectional `kph=1` family, the `eta33` discrepancy grows
   monotonically as steepness increases and decreases as steepness decreases.
2. The six-point trend is exceptionally well represented by a constant floor
   plus an `O((Akp)^2)` relative contribution. This is compatible with omitted
   fifth-and-higher odd-order content relative to a cubic signal, but is not a
   unique proof of that mechanism.
3. GL quadrature rank is also a material error source. The experiment separates
   it from auxiliary-domain convergence, but does not yet establish the
   infinite-rank limit.

The next discriminating numerical step, if required, is an independently
certified GL10/12 nested-graph extension on these frozen inputs. A water-depth
sweep should come later and retain the same nondimensional spectrum and
steepness ladder.

## Evidence

Local ignored artifacts:

- `artifacts/unidirectional_time_series/hos-eta33-steepness-20260929-v1/analysis/`
- `artifacts/unidirectional_time_series/hos-eta33-steepness-20260929-v1/rank_analysis/`

Tracked reproduction entries:

- `deploy_hos_eta33_steepness.py`
- `collect_hos_eta33_steepness.py`
- `analyze_hos_eta33_steepness.m`
- `analyze_hos_eta33_rank_steepness.m`
- `plot_hos_eta33_steepness_experiment.m`
