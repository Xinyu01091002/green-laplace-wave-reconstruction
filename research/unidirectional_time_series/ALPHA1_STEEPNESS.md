# Alpha=1: larger steepness and historical VWA time-series checks

The earlier kh=0.5 and kh=1 pilots already used Alpha=1, Akp=0.02. The user's
follow-up is interpreted as an Alpha=1 test at larger steepness. No completed
dense Alpha=1, kh=0.5, Akp>0.02 family was established in the VWA directory
search. A separate existing raw OW3D probe dataset supplies Alpha=1, kh=1,
Akp=0.02 and 0.12, so those two cases are compared within the same campaign.

## Data and execution

Read-only source:
`C:/Users/spet5947/Documents/ChatGPT/ESC-Time Series/outputs/remote_runs/esc-ow3d-full-20260827T015024Z/compact-results/compact`.
This is a different campaign/probe from the earlier VWA-folder pilot.
The campaign's preserved 641-sample windows are used; the handoff notes that
some runs are complete only over the declared window. No broader completion
claim is made. Every CSV SHA-256 must match `probe_manifest.json` before use.
Check physical parameters against the original `vwa_general_mf12_fourphase`
input files and config. This legacy dataset is retained explicitly as a
diagnostic, not the preferred finite-depth validation dataset (see below).

Only raw eta/psi CSV files and metadata are imported. No ESC prediction,
inverse solution, corrected GL executor or high-order reference is used.
The existing Hilbert-four-phase separation, common parent-domain/Nyquist
projection, ranks 6/8/12, MF12, VWA and Walker formulas remain unchanged.
The main-group window is +/-2Tp about the common first-order envelope maximum,
chosen from the input and applied identically to all methods. Both full-window
and main-group metrics are retained. Output directories carry `compact_` to
prevent overwriting the earlier VWA-folder results.

```matlab
run_alpha1_steepness( ...
    'C:/Users/spet5947/Documents/ChatGPT/ESC-Time Series', ...
    'C:/Research/spectral domain implementation of wave interaction theory',true);
```

The third argument selects the already existing boundary-corrected campaign
(true, default). Set false to reproduce the legacy dataset without overwriting
the corrected outputs. No new OW3D campaign is launched.

This is a test with a projected first phase-sector record, not a claim of
exact eta11 isolation. At larger steepness both higher-order primary-harmonic
content and phase aliasing can affect the input, while the observed second
phase sector can contain higher-order contributions. No cutoff is tuned
against the reference to remove that discrepancy.

## What was found in the old VWA code

The existing `paperplot_VWA_time_series.m` uses dt=0.15 and Tp=12 for a case
whose source metadata report dt=0.150125 and Tp=13.761997. This mismatch was
already corrected in our first pilot. It is a reproducibility issue, not
evidence by itself that it caused a particular historical discrepancy.

A stronger inconsistency appears in
`C:/Research/VWA/VWA time series/UWA_data_20260401/matlab/diagnose_timeseries_vs_spatial_vwa.m`.
The diagnostic declares kh=0.5, h=kh/kp and a fixed probe, but calls
`C:/Research/VWA/VWA Unidirectinal/Existing Theory/eta/Linear_focus_envelope.m`.
That function hard-codes h=150 for the group speed, uses `w=sqrt(g*k)`, and
uses `phi=-w*t+k*t*cw`. Consequently its temporal signal follows a moving
coordinate convention, with frequency w-k*cw, whereas the diagnostic's
temporal VWA reconstruction inverts stationary finite-depth dispersion at
h=0.5/kp. This is not a consistent same-physics fixed-probe comparison.
No historical file was edited, and we have not identified whether this
specific script produced the plot recalled by the user.

Finally, VWA and the pairwise GL/MF12 routes do not use the same interaction
approximation. VWA applies a one-parent spectral multiplier before multiplying
by the analytic first-order signal; GL/MF12 retain pair-dependent sums. Alpha
changes spectral shape, so comparable performance on one narrow/small-amplitude
case does not establish broad-spectrum or large-steepness accuracy. Distinguish
this approximation error (VWA versus MF12 on identical input) from the combined
input/model/reference discrepancy (each method versus the OW3D phase sector).

## Legacy finite-depth boundary problem

The first execution reproduced approximately 53% full-window discrepancy
for *all* GL/MF12/VWA routes, already at Akp=0.02. The local report
`ESC-Time Series/docs/22_ow3d_boundary_corrected_bound_harmonic_report.md`
documents that the legacy kh=1 initial packet intersected the closed-tank
wall: the peak packet centre began only about 2.97 wavelengths from it.
The report records an independent rerun with focus_fraction changed from
0.5 to 0.67, giving about 14.53 wavelengths of initial clearance.

That existing raw replacement is
`outputs/remote_runs/esc-ow3d-kd1-boundary-corrected-20260901T013632Z/compact`,
with source inputs under `vwa_general_mf12_fourphase_kd1_boundary_corrected`.
It has its own manifest hashes, which are verified before execution.
The present comparison uses the *same* GL/MF12/VWA/Walker formulas and common
support rule on both campaigns. It imports no Stokes-informed inversion or
other correction discussed elsewhere in that historical report. The new
outputs have `boundary_` names; the old `compact_` outputs are preserved.
This dataset correction is source-supported, not an adjustment of candidate
coefficients based on reference error. It does not establish the cause of
the separate kh=0.5 VWA-folder tail seen in the earlier pilot.

## Results with the unchanged reconstruction

All 16 CSV files (two amplitudes, four phases, two campaigns) passed their
respective manifest SHA-256 checks. All records have 641 finite samples;
four-phase grids, source dt/h/g and probe locations agree within each group.
The boundary-corrected probe is x=10259.18359375 m; dt_output=0.34405 s.
Modified MATLAB files passed Code Analyzer with no messages.

Raw relative L2 against the boundary-corrected OW3D second phase sector,
within +/-2Tp of the common first-order envelope maximum:

| Method | Akp=0.02 (%) | Akp=0.12 (%) |
|---|---:|---:|
| GL6 | 0.259207 | 3.934372 |
| GL8 | 0.122600 | 3.843050 |
| GL12 | 0.103794 | 3.848632 |
| Spectral MF12 | 0.103431 | 3.856581 |
| VWA | 0.899913 | 3.304667 |
| Walker | 11.031488 | 8.222252 |

The legacy Akp=0.12 main-group MF12/VWA errors were 52.5645%/52.9670%.
The legacy Akp=0.02 already had about 53% full-window MF12 error. Preserved
legacy results and the documented independent boundary correction therefore
identify a major dataset defect, rather than evidence of a VWA-only failure.

Boundary-corrected input projection retains 99.9999974073% / 99.9999958057%
of positive-frequency energy for Akp=0.02 / 0.12. GL12-versus-MF12 full-window
relative L2 is 0.00114451% / 0.0327904%. The current no-Stokes GL graph remains
close to the quadratic spectral reference as amplitude increases. No rank,
formula or support was selected against the numerical reference.

The large-steepness VWA-to-OW3D error is slightly smaller than MF12-to-OW3D,
but this does not certify VWA as the more accurate quadratic kernel. The
OW3D phase-sector comparison also contains higher-order input/output content;
approximation error can partially cancel that discrepancy. Its precise
decomposition has not been measured here. The remaining few-percent error
is not addressed by adding historical Stokes or fitted corrections.

The separate saved `alpha1_steepness_contrasts.csv` makes the distinction
explicit: VWA versus MF12 on the identical main-group input is 0.97035%
at Akp=0.02 and 0.86254% at Akp=0.12. GL12 versus MF12 is 0.0011435%
and 0.032791%. Thus the dramatic legacy OW3D mismatch is not explained
by the roughly one-percent VWA-versus-MF12 approximation difference here.
