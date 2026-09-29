# Random-wave periodic-boundary check, 2026-09-29

## Question

The random-wave third-harmonic discrepancy might conceivably arise from a
finite HOS domain or a reflected wave. The production input, however, uses a
Fourier periodic domain rather than walls. No `relaxation zones` block is
present in `input.yml`; the HOS parser therefore keeps its default `i_rlx=0`,
and the Runge--Kutta code does not call the relaxation-zone operator.

## Exact repeated-field double-domain test

The most discrepant realization, seed `20260927` at nominal `Akp=.12`, was
rerun on an exactly doubled domain:

| Quantity | Original | Doubled |
|---|---:|---:|
| Domain length | 68 lambda_p | 136 lambda_p |
| Fourier nodes | 4096 | 8192 |
| Delta x | unchanged | unchanged |
| HOS order | 5 | 5 |
| Duration | 50 Tp | 50 Tp |

Each original phase field was repeated exactly twice, so only the containing
periodic domain changed. Two doubled-domain probes were placed at `x_p` and
`x_p+L_original`. Four HOS phase runs completed with exit code zero using the
same solver binary SHA-256
`13c4bae40855f3d3a29da12d5567559c4a41e415c27bc5f67d5b927db6c2f060`.

The initial deployment directory `hos-random-double-domain-20260929-v1`
contains only a preserved controller syntax failure and no HOS result. The
successful isolated run is
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-random-double-domain-20260929-v2`.

## Result

Raw phase-record relative differences over the complete record are
`2.70e-13--6.89e-13` between the doubled and original domains. The two probes
inside the doubled domain agree to approximately `7.84e-14--8.14e-14`.

After the established `5 Tp` edge taper and four-phase extraction, central
15--35 Tp differences are:

| Harmonic | Doubled versus original | Probe x_p+L versus x_p |
|---|---:|---:|
| Second harmonic | 7.7347e-12 | 1.3186e-12 |
| Third harmonic | 4.6197e-10 | 4.1070e-11 |

The larger relative number for the third harmonic reflects its much smaller
norm; the absolute waveform differences remain about machine/accumulated
roundoff scale.

## Conclusion

The doubled periodic domain reproduces the original local evolution and the
repeated probe exactly for practical purposes. There is no wall reflection in
the configured HOS problem, and periodic wrap-around/domain length does not
explain the remaining random-wave third-harmonic difference.

The remaining investigation should focus on the current first-order state,
native spatial GL evaluation, and free-versus-bound third-order evolution.

## Evidence

Local ignored outputs:
`artifacts/unidirectional_time_series/hos-random-double-domain-20260929-v2/results/`.

Tracked reproduction entries:

- `prepare_hos_random_double_domain.m`
- `analyze_hos_random_double_domain.m`
- `deploy_hos_random_double_domain.py`
