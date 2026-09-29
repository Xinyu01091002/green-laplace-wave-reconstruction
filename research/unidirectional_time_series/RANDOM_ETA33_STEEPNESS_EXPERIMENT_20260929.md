# Unidirectional random-phase eta33 steepness experiment, corrected 2026-09-29

## Correction to the first analysis

The first analysis applied an FFT Hilbert transform directly to stationary
random-wave records whose two ends were nonzero and discontinuous under
periodic extension. That analysis produced apparent HOS third harmonics
with amplitude powers `1.06--1.54` and raw GL--HOS errors near 100%. A pure
linear four-phase control produced almost exactly the same spurious third
harmonic: at low amplitude its leakage/HOS-third-harmonic norm ratio was approximately
one and its waveform correlation was approximately one.

Those untapered random-wave results are therefore invalid as `eta33` evidence.
They diagnose finite-record Hilbert leakage, not missing HOS `eta33` and not GL
physical error.

Following the user's correction, the final procedure applies the same
raised-cosine taper to both ends of all four phase records before the Hilbert--
four-phase separation and evaluates only a central interval. A `5 Tp` taper at
each end of the complete `50 Tp` record is retained; the scoring interval is
`15--35 Tp`, at least `10 Tp` away from either tapered edge.

## Frozen random-wave campaign

The experiment uses the same unidirectional semi-Gaussian modal magnitudes as
the focused family, but independent uniform positive-mode phases. Three
reproducible MT19937 realizations use seeds `20260925`, `20260926` and
`20260927`; each seed retains the same phases across the amplitude ladder.

```text
nominal Akp = 0.02, 0.04, ..., 0.18
kph = 1
```

No Hs renormalization is applied. The actual random-sea steepness is:

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

Each seed/amplitude uses four global phases initialized with MF12
`11+20+22+33`, no `31`, `muStar=0`, and linear frequencies. HOS uses order 5,
50 Tp duration and Tp/40 output. A highest-amplitude smoke test passed, and all
108 production HOS phase runs completed with exit code zero. Solver SHA-256:
`13c4bae40855f3d3a29da12d5567559c4a41e415c27bc5f67d5b927db6c2f060`.

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-eta33-random-steepness-20260929-v1`.

## Taper selection and linear control

Taper widths `0, 2, 5, 10 Tp` per edge were tested with a fixed central
`15--35 Tp` diagnostic. The HOS third-harmonic norm powers are:

| Taper per edge | Seed 20260925 | Seed 20260926 | Seed 20260927 |
|---:|---:|---:|---:|
| 0 | 1.155 | 1.193 | 1.607 |
| 2 Tp | 2.565 | 2.676 | 2.933 |
| 5 Tp | 2.975 | 2.977 | 2.999 |
| 10 Tp | 2.999 | 2.998 | 2.987 |

Thus a `5 Tp` taper already restores essentially cubic scaling while changing
less of the record than `10 Tp`. At `5 Tp`, the pure-linear spurious-third norm
is only `0.16--1.52%` of the HOS third-harmonic norm, compared with `65--95%`
without taper. This establishes that the same HOS solver produces a normal
cubic third harmonic for the random waves once the finite-record extraction
is treated correctly.

## Corrected tapered GL8 comparison

GL is recomputed from the tapered first-harmonic record over the complete
50 Tp interval; only 15--35 Tp is scored. The comparison uses the
non-enumerating `gl_unidirectional_time_series` GL8 implementation, no fitted
adjustments and no Stokes correction.

| Nominal Akp | Actual kp Hs/2 | Median error | Three-seed range | Median cosine | Median GL/HOS norm |
|---:|---:|---:|---:|---:|---:|
| 0.02 | 0.003243 | 60.83% | 51.54--71.61% | 0.794 | 0.784 |
| 0.04 | 0.006486 | 53.25% | 40.86--71.55% | 0.847 | 0.835 |
| 0.06 | 0.009729 | 52.70% | 40.03--71.55% | 0.850 | 0.838 |
| 0.08 | 0.012973 | 52.58% | 39.84--71.55% | 0.851 | 0.838 |
| 0.10 | 0.016216 | 52.53% | 39.75--71.55% | 0.851 | 0.838 |
| 0.12 | 0.019459 | 52.49% | 39.69--71.56% | 0.851 | 0.837 |
| 0.14 | 0.022702 | 52.44% | 39.62--71.56% | 0.852 | 0.837 |
| 0.16 | 0.025945 | 52.39% | 39.54--71.57% | 0.852 | 0.836 |
| 0.18 | 0.029188 | 52.33% | 39.46--71.59% | 0.852 | 0.836 |

After `Akp=.04`, each seed's discrepancy is nearly amplitude-independent. The
random-wave result therefore does **not** show the steepness-driven growth seen
for the focused wave group. It shows a strong realization-dependent systematic
difference: approximately `39--40%`, `52--53%`, and `71.5%` for the three
seeds. The HOS third harmonic and GL waveform correlations are approximately `0.92`,
`0.85`, and `0.70` respectively.

## Numerical bounds

The taper reduces the first-input projection error to only
`0.0011--0.0039%`. Thus the earlier `2.75--4.63%` projection error was also a
consequence of using the untapered finite record.

At the final `L/h=1600`, successive-domain changes are `0.359--0.478%`; none
of the 27 cases passes the stricter predeclared `0.05%` domain gate. The
corrected GL8 curves are therefore not called spatially converged. However,
this sub-percent numerical uncertainty cannot explain the stable
`39--72%` seed-dependent differences.

### Fixed-input rank check

At nominal `Akp=.12`, the same tapered inputs and HOS third harmonics were evaluated
with shared inner/outer GL4, GL6 and GL8:

| Seed | GL4 error | GL6 error | GL8 error |
|---:|---:|---:|---:|
| 20260925 | 53.59% | 52.60% | 52.49% |
| 20260926 | 42.28% | 40.00% | 39.69% |
| 20260927 | 72.47% | 71.67% | 71.56% |

The waveform correlations change by less than about `0.0012` from GL4 to
GL8. Rank therefore contributes a few percentage points at most in this
diagnostic and cannot explain the stable `39--72%` seed-dependent difference.
GL8 is not a mathematically certified infinite-rank limit, but increasing rank
alone is unlikely to close the discrepancy.

## Interpretation

The corrected evidence supports the following bounded conclusions:

1. The same HOS code produces a normal cubic third harmonic for focused
   and random waves. The earlier claim that the random HOS third harmonic was not cubic
   was an analysis error caused by untapered finite-record Hilbert leakage.
2. Extending the nominal scale to `Akp=.18` is numerically stable for all three
   seeds, corresponding here to actual `kp Hs/2<=.02919`.
3. After proper tapering, the GL8 discrepancy is largely independent of
   steepness but strongly dependent on the random phase realization.
4. The fixed-input GL4/6/8 diagnostic rules out low GL rank as the dominant
   cause. The remaining difference may include native spatial representation,
   free versus bound third-order evolution, and sensitivity of broadband
   interaction accumulation to the phase realization. It is not evidence of a
   simple high-steepness breakdown.

The next discriminating step is a native-spatial/current-state comparison on
the three frozen seeds, not a further steepness extension.

An exact repeated-field domain-doubling check has also ruled out spatial
periodicity as the cause: changing from `68 lambda_p / 4096` to
`136 lambda_p / 8192` changes the central second and third harmonics by only
`7.73e-12` and `4.62e-10` relative, respectively. See
`RANDOM_PERIODIC_BOUNDARY_CHECK_20260929.md`.

## Evidence

Local ignored outputs:

- `artifacts/unidirectional_time_series/hos-eta33-random-steepness-20260929-v1/taper_diagnostic/`
- `artifacts/unidirectional_time_series/hos-eta33-random-steepness-20260929-v1/tapered_gl/`
- `artifacts/unidirectional_time_series/hos-eta33-random-steepness-20260929-v1/hilbert_leakage/`
- `artifacts/unidirectional_time_series/hos-eta33-random-steepness-20260929-v1/boundary_plot/`
- `artifacts/unidirectional_time_series/hos-eta33-random-steepness-20260929-v1/tapered_rank/`

Tracked reproduction entries include:

- `prepare_hos_random_eta33_steepness.m`
- `deploy_hos_random_eta33_steepness.py`
- `diagnose_random_hilbert_leakage.m`
- `diagnose_random_tapered_harmonics.m`
- `analyze_hos_random_eta33_tapered_gl.m`
- `plot_random_eta33_boundary_effect.m`
- `analyze_hos_random_eta33_tapered_rank.m`
- `collect_hos_random_eta33_steepness.py`
