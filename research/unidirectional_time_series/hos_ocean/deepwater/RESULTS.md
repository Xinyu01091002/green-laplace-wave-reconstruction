# Results: unidirectional deep-water four-phase HOS / spectral MF12

Completed 2026-09-26. Final accepted runs are `tol16/akp002` (absolute
RK tolerance 1e-16) and `tol14/akp012` (1e-14). The physical setup and
frozen source hashes are in README.md. Both four-phase families reached
20Tp = 240.2001413833 s. Every output record was finite and t=0 eta/psi
matched prescribed input to relative L2 <1e-12. No reference alignment,
amplitude fit, phase shift or offset correction was applied.

## Spatial reconstruction errors

Percent relative L2 over the full periodic domain:

| Akp | Time after initialization | second harmonic | third harmonic |
| --- | --- | ---: | ---: |
| 0.02 | 3Tp | 0.009535% | 7.96591% |
| 0.02 | 20Tp | 0.026933% | 4.94667% |
| 0.12 | 3Tp | 0.358064% | 8.66296% |
| 0.12 | 20Tp | 1.032746% | 6.74302% |

At t=0 the second/third reconstruction errors are approximately 7.70e-14 /
1.58e-11 for Akp=.02 and 1.30e-14 / 4.41e-13 for Akp=.12 (fractions, not
percent). This is initialization/separation consistency, not independent
physical validation.

The second harmonic is reconstructed closely, with increased discrepancy
at higher steepness and later evolution. Third-harmonic-sector reconstruction
improves between 3Tp and 20Tp but retains a 5--9% full-domain discrepancy.
The experiments do not demonstrate a unique adjustment-completion time.
Ta=0 throughout, and both snapshots precede the nominal linear focus at 40Tp.

The third record is the four-phase plus positive-spatial-wavenumber projection,
not an independently certified eta33 reference. HOS can generate primary
harmonic corrections despite excluding them at initialization; four-phase
aliasing can mix +3 with -1 phase sectors. Spatial direction projection does
not eliminate every possible mixed-sign or counterpropagating contribution
of a broadband evolving packet. This remains a candidate explanation, not a
proven diagnosis. A bounded eight-phase audit would distinguish the phase
aliasing contribution before blaming MF12 or an initialization transient.

The first-harmonic projection outside the frozen parent support is small:
relative norm at 20Tp is 1.97e-7 (low) and 1.97e-5 (high). A diagnostic of
third-sector content below the minimum positive triple-sum wavenumber found
only about 1e-6 relative norm. Thus that simple disjoint-band check does not
explain the observed third-sector discrepancy; overlap-region phase aliasing
is not ruled out. No diagnostic cutoff was used to alter the reported errors.

## Time integration audit

Absolute tolerance 1e-10 was inadequate for the low-amplitude harmonics:
its second-sector errors were about 11% and 21%. These failed precision
results are retained and are not the final physical result.

Final output differences from the preceding tighter-run pair, relative to
the observed harmonic norm:

| Akp | Tolerance comparison | Time | second | third |
| --- | --- | --- | ---: | ---: |
| .02 | 1e-14 vs 1e-16 | 3Tp | 1.789e-6 | 4.383e-6 |
| .02 | 1e-14 vs 1e-16 | 20Tp | 2.969e-6 | 2.417e-6 |
| .12 | 1e-12 vs 1e-14 | 3Tp | 4.486e-6 | 1.138e-5 |
| .12 | 1e-12 vs 1e-14 | 20Tp | 9.016e-6 | 7.732e-6 |

These are fractions. They support temporal precision for the main conclusions;
spatial-grid, HOS-order and infinite-depth-limit convergence have not been
certified. The experiment uses finite kph=20 with minimum parent kh=5.588.
Long difference modes need not themselves satisfy a deep-water limit.

## Measured resources

Wall time and per-process maximum RSS (GNU time), not estimates:

| Stage | Akp=.02 | Akp=.12 |
| --- | --- | --- |
| Prepare four MF12 phase fields, including MATLAB startup | 31.24 s; 1.587 GiB | 25.29 s; 1.599 GiB |
| Final HOS, each phase, one thread | 214.15--219.07 s; 18.344 MiB | 71.64--72.48 s; 18.348 MiB maximum |
| Final MATLAB reconstruction, metrics and base plots | 43.42 s; 1.624 GiB | 42.99 s; 1.611 GiB |

The four HOS phases in each batch run concurrently. The low final run costs
more because its absolute tolerance is 1e-16 versus 1e-14 for high, not because
low-steepness physics is intrinsically slower. These are 1D 4096-point runs;
do not compare their timing to the earlier 2D pilot as a solver speedup.

Including the initial and tolerance-audit batches, 28 HOS processes completed:
8 at 1e-10, 8 at 1e-12, 8 at 1e-14 and 4 low-steepness runs at 1e-16.
All stage resource files are preserved. `performance.csv` is a file-level
inventory, not a sum of elapsed time or aggregate peak memory. It contains
preserved analysis-log copies; do not double-count them in totals.

## Artifacts

Remote run root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-deepwater-fourphase-20260926-v1`.

Final native fields and MATLAB comparison arrays remain remote. Locally,
small summaries and plots are under `artifacts/hos_ocean/deepwater-results/`:

- `tol16/akp002/comparison_zoom.png`, `comparison.png`, `metrics.csv`, `analysis.json`.
- `tol14/akp012/comparison_zoom.png`, `comparison.png`, `metrics.csv`, `analysis.json`.
- `performance.csv`, `final-summary.json`, `support_and_tolerance_audit_final.csv`.

The zoom plot uses one common coordinate origin from the first-harmonic
envelope for both curves; this is a display coordinate, not fitted alignment.
The full-domain figures and errors are preserved separately. A summary-packaging
script syntax failure was corrected in `final-processing/collect_results_fixed.py`;
it did not modify or rerun the numerical results. Initial and refined source
snapshots, logs, MF12 archive and all native solver outputs remain available.