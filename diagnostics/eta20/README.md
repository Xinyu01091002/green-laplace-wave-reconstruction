# Nonzero difference-frequency eta20 diagnostics

This directory is deliberately outside the released total-field API. It
contains two distinct approximation families for the strictly nonzero
difference-frequency elevation `eta20`.

## Shared-scale Green--Laplace ranks

`eta20_green_laplace_shared` evaluates the frozen shared-scale graph. GL6,
GL12 and GL16 change only the prescribed standard Gauss--Laguerre rank. The
forcing graph and physical scale remain fixed; MF12 fields do not select or
fit the nodes, weights, formula or rank.

For the committed `k_p d=1`, 120-degree crossing and 30-degree spreading
example, the raw relative L2 values are 23.827%, 9.609% and 7.824% for GL6,
GL12 and GL16. At `k_p d=2` they are 7.890%, 2.235% and 1.003%. These are
fixed-input rank diagnostics, not uniform weak-detuning guarantees.

## Neumann R2/R4/R6 sequence

`eta20_neumann_r_series_ordered_pair` evaluates exact frozen R2, R4 and R6
formulas with an independent ordered-pair reconstruction. The Wolfram freeze
certifies the inverse, residual, endpoint and compiler identities without
MF12 or sampled selection.

For the same `k_p d=1` input, raw relative L2 values are 1.766%, 0.948% and
0.407% for R2, R4 and R6. The strict zero-output mode is set to zero and
belongs to the separate mean, flux and gauge sector.

The committed R evaluator is accuracy-validation code:

- `pair_loops = 2`;
- `production_candidate = false`;
- `fixed_fft_validated = false`.

The exact compiler counts for the corresponding proposed fixed graphs are:

| Method | FFT/IFFT | Products | Parent filters |
|---|---:|---:|---:|
| R2 | 53 | 124 | 47 |
| R4 | 87 | 316 | 79 |
| R6 | 121 | 596 | 111 |

These counts do not imply that the high-power fixed-FFT implementations are
numerically validated or timed.

## Run

Set the external MF12 repository and run from the repository root:

```matlab
setup_green_laplace("MF12Root",getenv("MF12_ROOT"));
addpath("diagnostics/eta20");
figure_eta20_gl_rank("results/paper",getenv("MF12_ROOT"),1);
figure_eta20_gl_rank("results/paper",getenv("MF12_ROOT"),2);
figure_eta20_r_series("results/paper",getenv("MF12_ROOT"),1);
```

All comparisons are raw same-unit fields with no gain, alignment, shift or
fitted rescaling.
