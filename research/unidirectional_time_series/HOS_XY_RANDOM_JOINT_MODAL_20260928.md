# Akp=.12 off-axis wave-group and random-phase GL checks

## Focused wave group at new x-y positions

The isolated HOS replay `hos-directional-high-xy-probes-20260928-v1`
completed all four global phases with 1101 finite samples at 12 probes. Every
initial-value check passed to `6.1e-15 m` or better. Four probes are same-x
centreline references; eight probes have both nonzero longitudinal and lateral
offsets.

The eight off-axis positions span

```text
dx/lambda_p = +/-0.4883, +/-0.9766
abs(dy/lambda_p) = 0.2344, 0.4688, 0.7813.
```

Measured HOS first-harmonic maximum ratios relative to the corresponding
same-x centreline are 0.429--0.875. Phase-zero total-wave maximum ratios are
0.386--0.854. All points pass the predeclared one-third amplitude gate.

The maintained non-enumerating joint-modal results use 20 adaptive observed
frequencies, eta20 GL16/512^2, eta22 GL8/256^2 and eta33 GL10/256^2. Raw,
unaligned main-window HOS phase-sector errors are:

| Component | Minimum | Maximum |
| --- | ---: | ---: |
| nonzero eta20 | 2.05% | 6.41% |
| positive-sum eta22 | 2.15% | 5.39% |
| positive-sum eta33 | 4.21% | 8.94% |

The spatially displaced wave-group result passes: no component shows a
material collapse relative to the original same-x five-probe evidence. The
largest eta22/eta33 errors occur at the farthest `(+0.9766,+0.7813)` point,
whose measured first-harmonic ratio is still 0.429.

## Existing phase-only random HOS

Because the focused positional gate passed, the existing completed random
family `hos-random-fourphase-20260926-v1/akp012` was reprocessed without
rerunning HOS. It preserves the same 3777 modal amplitudes and randomizes only
their phases with seed 20260925. Its actual linear values are

```text
Hs = 0.4076151725 m
kp*Hs/2 = 0.005686231656.
```

The old 17-frequency GL result was not reused. Adaptive 99.9%-energy selection
retained 31--72 contiguous frequencies across the five probes. The resulting
first-record projection error is 3.09--3.16%, much larger than the focused
wave-group projection error despite meeting the energy threshold.

Raw fixed-window HOS phase-sector errors are:

| Component | Minimum | Maximum | Norm-ratio behaviour |
| --- | ---: | ---: | --- |
| nonzero eta20 | 19.92% | 35.79% | 0.92--1.01 |
| positive-sum eta22 | 3.90% | 10.44% | 1.01--1.08 |
| positive-sum eta33 | 99.97% | 100.02% | 0.01--0.03 |

The random-wave result is therefore mixed rather than a validation pass.
Eta22 remains qualitatively close, though less uniform than the focused wave
group. Eta20 has the right norm scale but substantial waveform disagreement.
The GL eta33 prediction is far smaller than the HOS Hilbert third phase
sector; this may include finite-record phase leakage, free/transient content
or nonlinear evolution and is not by itself proof of a kernel defect.

No amplitude, phase, time or spatial alignment was fitted. HOS references are
reported as phase sectors rather than certified perturbation orders.

## Evidence

Compact outputs are under
`artifacts/unidirectional_time_series/hos-xy-random-joint-modal-v1/`.
Focused HOS probe records and terminal status are under
`artifacts/unidirectional_time_series/hos-directional-high-xy-probes-v1/completed/`.
Raw full fields were not downloaded.
