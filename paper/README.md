# Paper figures

The portable manuscript entry points are:

- `figure_01_scalar_convergence.m`: scalar modal-resolvent convergence;
- `figure_02_input_cases.m`: the three directional input cases;
- `figure_03_eta22_rank_ladder.m`: the GL2/4/6/8 field ladder;
- `figure_04_eta22_mf12_waveform.m`: the representative GL6/MF12 field;
- `figure_05_psi33_mf12_waveform.m`: the ranked GL/MF12 surface-potential
  field, requiring the separately published large-field archive.

Each computed figure writes its PNG and available CSV/MAT provenance under
`results/paper/`.

`run_all_figures("self-contained")` regenerates Figures 1--2 and the compact
unified-surface API figure without MF12. `run_all_figures("full")` also runs
the MF12-dependent Figures 3--4. Figure 5 is added when
`GL_PSI33_DATA_ROOT` points to the external matched-field archive.

The `reference/` directory preserves the current manuscript-facing PNG and
CSV assets for provenance. The second-order difference-frequency Figures
6--7 remain reference-only in version 0.1 because that sector is not part of
the released unified surface API. Some field figures depend on large MF12
fields that are not suitable for Git. A reference raster is never presented
as an end-to-end numerical reproduction.

All new figure calculations write to `results/`, which is ignored by Git.
Committed files in `reference/` must not be edited by hand.
