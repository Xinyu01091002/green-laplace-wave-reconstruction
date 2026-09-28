# JONSWAP kpHs/2=.06 joint-modal checkpoint, 2026-09-28

## Input and scope

This checkpoint reuses the completed four-phase CPU HOS family
`hos-jonswap-kphs006-20tp-20260927-v1/medium`. It has 1377 samples over 20Tp,
five probes, JONSWAP gamma 3.3, 11128 retained parents, seed 20260925 and
`kp*Hs/2=.06`. HOS is not rerun.

The current non-enumerating joint-modal executors use the observed first
harmonic plus the initial directional spectrum. No amplitude, phase, time or
spatial alignment is fitted. HOS outputs remain phase sectors rather than
certified perturbation orders.

The eta20 GL16 result in this checkpoint is retained as a diagnostic only.
The user selected native R4 as the next eta20 route; do not promote the GL
eta20 diagnostic as the maintained recommendation.

## First observed-spectrum run and timing

The first 99.99%-energy attempt stopped explicitly because the selected cubic
baseband exceeded the native temporal record. This demonstrates that the weak
finite-record spectral tail cannot be included blindly. Its composition may
include physical support, nonlinear broadening, phase-sector leakage and
rectangular-window leakage; it is not certified nonphysical.

The completed 99.9%-energy run retained 61--82 contiguous frequencies. Its
first-record projection error is 3.13--3.16%. Raw fixed-window HOS errors are:

| Component | Minimum | Maximum |
| --- | ---: | ---: |
| eta20 GL16 diagnostic | 40.32% | 117.96% |
| eta22 GL8 | 11.57% | 29.33% |
| eta33 GL10 | 76.07% | 90.38% |

Four-thread median call times are:

| Stage | Median seconds |
| --- | ---: |
| observed-input preparation | 0.233 |
| eta20 GL16, 512x512 | 39.03 |
| eta22 GL8, 256x256 | 4.88 |
| eta33 GL10, 256x256 | 19.37 |
| three GL calls combined | 63.54 per probe |

The five-probe numerical loop took 317.12 s, excluding MATLAB startup and
figure export.

## Grid and rank separation

Probe 1 was rerun with fixed observed input and fixed physical modal limits.

| Component | 256-to-512 change | HOS error at 256 | HOS error at 512 |
| --- | ---: | ---: | ---: |
| eta22 GL8 | 6.83% | 29.33% | 28.38% |
| eta33 GL10 | 27.30% | 90.38% | 90.74% |

Thus 256x256 is not adequate for this broad JONSWAP input, especially at
third order. Refinement does not remove the main HOS discrepancy.

At fixed 256x256, changes relative to GL12 are:

| Component | GL6 | GL8 | GL10 |
| --- | ---: | ---: | ---: |
| eta22 | 6.36% | 3.88% | 1.82% |
| eta33 | 12.28% | 6.87% | 3.15% |

Rank error is not negligible, but the eta33 grid change is substantially
larger than the GL10-to-GL12 change.

## Native initial-wavenumber-grid experiment

The initial spectrum already lies on the HOS 1024x512 Fourier lattice. A new
eta22 experiment retains those wavevectors exactly and applies observed
first-harmonic information only through temporal-frequency conditioning.

A common multiplicative complex correction per frequency shell is unstable:
directional cancellation gives energy-weighted condition numbers 21--38 and
produces 42%--463% HOS errors.

A minimum weighted-norm additive correction avoids whole-shell amplification
and exactly restores each retained station bin. It gives:

| Probe | Native-grid eta22 error | Prediction/reference norm |
| ---: | ---: | ---: |
| 1 | 19.57% | 1.006 |
| 2 | 17.42% | 0.926 |
| 3 | 25.09% | 0.935 |
| 4 | 16.34% | 0.929 |
| 5 | 23.08% | 0.927 |

Median runtime is 18.53 s per probe. The native input-support bins 12--49
retain only 99.45--99.73% of the observed first-record energy, giving
5.18--7.42% first-record projection error. Removing spatial wavevector
projection therefore exposes a separate state-estimation limitation.

## Interpretation and next boundary

For a unidirectional wave, known propagation direction plus dispersion maps a
single probe's temporal spectrum to its spatial spectrum. For a directional
focused group, the observed probe plus a trusted initial directional prior can
conditionally close the problem. For a directional random sea, one complex
station coefficient per frequency cannot uniquely determine the many modal
directions and phases that produce it.

The Green--Laplace/R4 forward problem and the directional-state inverse problem
must therefore be kept separate:

```text
state estimation: probe records -> current C(kx,ky)
forward reconstruction: current C(kx,ky) -> bound eta20/eta22/eta33
```

The next eta20 implementation should adapt SPARK's native fixed-FFT R4 graph to
the joint time representation, preserving its targeted low-output repair and
strict-zero boundary. The next directional random-wave step should use several
probes jointly to regularize the current modal state rather than assigning an
independent direction spectrum from each single probe.

## Evidence

Compact outputs are under:

- `artifacts/unidirectional_time_series/hos-jonswap006-joint-modal-v2/`
- `artifacts/unidirectional_time_series/hos-jonswap006-numerical-convergence-v1/`
- `artifacts/unidirectional_time_series/hos-jonswap006-native-eta22-v1/`
- `artifacts/unidirectional_time_series/hos-jonswap006-native-eta22-v2/`

Remote source bundles record the exact executor hashes. Raw HOS fields were not
downloaded.
