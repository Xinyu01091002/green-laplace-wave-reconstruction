# Paper figures

The current third-order entry points use the no-Stokes runtime described in
[`docs/no_stokes_migration.md`](../docs/no_stokes_migration.md). Previously
committed reference assets are version-specific records and were not all
regenerated in this migration. Recompute a figure before citing it as a result
of the current runtime; a historical corrected inner state is not equivalent
to removing only the final third-order correction.

The portable manuscript entry points are:

- `figure_01_scalar_convergence.m`: scalar modal-resolvent convergence;
- `figure_02_input_cases.m`: the three directional input cases;
- `figure_03_eta22_rank_ladder.m`: the GL2/4/6/8 field ladder;
- `figure_04_eta22_mf12_waveform.m`: the representative GL6/MF12 field;
- `figure_05_psi33_mf12_waveform.m`: the ranked GL/MF12 surface-potential
  field, requiring the separately published large-field archive.
- `diagnostics/eta20/figure_eta20_gl_rank.m`: GL6/GL12/GL16
  difference-frequency rank diagnostic at `k_p d=1,2`;
- `diagnostics/eta20/figure_eta20_r_series.m`: GL12 versus the frozen
  R2/R4/R6 ordered-pair accuracy diagnostics;
- `figure_psi33_mixed_rank.m`: non-selecting inner/outer rank allocation.

Each computed figure writes its PNG and available CSV/MAT provenance under
`results/paper/`.

`run_all_figures("self-contained")` regenerates Figures 1--2 and the compact
unified-surface API figure without MF12. `run_all_figures("full")` also runs
the MF12-dependent eta22 and eta20 figures. The surface-potential figure and
mixed-rank diagnostic are added when
`GL_PSI33_DATA_ROOT` points to the external matched-field archive.

The `reference/` directory preserves the regenerated PNG and CSV assets for
provenance. Difference-frequency eta20 remains a diagnostic rather than a
component of `gl_spectral_surface`. R2/R4/R6 are complementary Neumann
approximations, not Green--Laplace ranks. Some third-order figures depend on
large MF12 fields that are not suitable for Git. Their hashes are recorded in
`data_manifest.json`.

The committed surface-potential Figure 5 uses the uncorrected pure-GL
resolvent field. For GL6 it has `Q=0.0346698` and raw relative L2 `6.75935%`
on the declared 8-wavelength domain. Corrected and uncorrected variants must
not be mixed in captions or metrics.

All new figure calculations write to `results/`, which is ignored by Git.
Committed files in `reference/` must not be edited by hand.
