# Fixed-location test of directional joint-input reconstruction

The user requested time records at several probe positions. Five locations
were fixed before reading any new high-order scores: domain midpoint, +/-one
peak wavelength in x, and +/-half a peak wavelength in y. Targets are snapped
to the nearest original grid point; no spatial interpolation or probe search
is used. The same kh=1, spread-label=25 degrees, Akp=0.02 four-phase dataset
and the same steps 0:10:900 (4 s interval) are retained.

| Probe | Requested offset in lambda_p | Actual x (m) | Actual y (m) |
|---|---|---:|---:|
| center | (0,0) | 4500 | 3375 |
| x_minus | (-1,0) | 4271.484375 | 3375 |
| x_plus | (1,0) | 4728.515625 | 3375 |
| y_minus | (0,-0.5) | 4500 | 3269.53125 |
| y_plus | (0,+0.5) | 4500 | 3480.46875 |

All four phases and 91 saved times are read once for all probes. The initial
spatial spectrum is shared. Each probe has its own observed eta1(t) and its
own initial-spectrum spatial phase at that coordinate. The same conditional
complex-direction allocation is then applied independently, with 15-degree
and 7.5-degree cells. Its parameters, GL ranks, support rules and denominator
handling are unchanged from `JOINT_INPUT_PILOT.md`.

The main-group window at each point is +/-2Tp about that point's observed
first-order envelope maximum. Every candidate and reference uses that same
local window. Different arrival times are neither fitted nor aligned across
probes. Full-window and main-group errors are both saved.

The center extraction and shared initial spectrum must be byte-for-byte
numerically identical to the preceding single-probe extraction. After this
assertion passes, the center's existing result and global initial-convention
audit are reused. Four additional probe reconstructions are actually executed;
the old center result is not overwritten. The observed high-order fields are
used for scoring only, never for position selection or directional weights.

```matlab
addpath('research/directional_wave_data');
run_directional_multiprobe('C:/Research/VWA/VWA Unidirectinal/Directional/test1');
```

Results are under `results/directional_joint_input_multiprobe/`, with a shared
source-file hash list, extracted records, per-probe result subdirectories,
`summary.csv`, `summary.mat`, and PNG/PDF multi-probe figures. The summary
points explicitly to the existing center result. This is still a test at
sparse saved times and does not certify a continuous high-harmonic record.

## Executed results

Center and source-spectrum parity assertions passed. The four new positions
ran with unchanged 444/888-parent angular discretizations and no target-led
parameter adjustment. Main-group relative L2 against OW3D, in percent:

| Probe | eta22, joint 15 deg | eta22, joint 7.5 deg | psi22, joint 15 deg | psi22, joint 7.5 deg |
|---|---:|---:|---:|---:|
| center (reused baseline) | 0.6458 | 0.2668 | 0.4031 | 0.1402 |
| x_minus | 0.6182 | 0.1911 | 0.6021 | 0.2137 |
| x_plus | 0.6366 | 0.2453 | 0.6447 | 0.2903 |
| y_minus | 1.0651 | 0.6168 | 1.6526 | 0.8794 |
| y_plus | 0.9539 | 0.5242 | 1.6532 | 0.8865 |

The direction-blind comparison ranges from 19.28% to 31.84% for eta22 and
13.75% to 28.29% for psi22. Using the initial spectrum alone gives 0.9439%--
2.7914% eta22 and 2.5148%--5.7853% psi22. These controls and full-window scores
are preserved in the per-probe metrics, not discarded after comparison.

Observed-energy-weighted allocation conditioning on the 7.5-degree grid is
about 1.00 at the center, 1.23 at the x-offset probes and 2.03 at the y-offset
probes. This shows increased directional cancellation away from the center;
the very weak tail-bin conditioning caveat from the original trial remains.

Each local main-group metric contains 13 saved samples. These additional
points support spatial robustness in this small-amplitude family, not dense
temporal resolution, arbitrary probe positions, angular convergence or
larger-steepness validity. Both lateral points and both x-offset points were
retained regardless of their errors.

The four modified/new MATLAB files pass Code Analyzer. The reused initial
spectrum, center records, time grid and h/g/kp metadata passed explicit parity
checks against the original single-probe extraction. No GL kernel was changed.
