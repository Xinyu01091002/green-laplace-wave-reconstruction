# GL comparison of the completed GPU low-family records

## Completed result

The lxy-account run completed at 2026-09-27 09:29:53 UTC, wall time
2060.44 s (34.34 min), sampled peak MATLAB RSS 13,211,760 KiB (12.60 GiB).
Local metrics.csv, report.json and rendered full_comparison.png are in the
deployment artifacts directory below. All five probes passed the original
one-third amplitude eligibility rule.

At direction resolution 1.875 degrees, GL16 and the fixed 10--70Tp window:

| y offset (m) | Raw relative L2 (%) | Common sum-band relative L2 (%) |
|---|---:|---:|
| 0 | 4.2147 | 4.2383 |
| -52.7821 | 3.2956 | 3.3010 |
| +52.7821 | 3.5131 | 3.5344 |
| -70.3762 | 3.1934 | 3.1939 |
| +70.3762 | 4.5602 | 4.5712 |

Direction refinement from 3.75 to 1.875 degrees changes predictions by
1.09--1.18% relative L2 in the interior. The fixed 10 s shorter center record
changes first input by 0.645%, GL prediction by 3.760%, and the filtered
reference by 0.458%. Thus record-length/frequency-grid sensitivity remains
material compared with the observed GL--HOS difference. This is not a
converged sub-percent validation or a pure perturbation-order reference.

## Readable GL--HOS time-series detail

`plot_gl_hos_detail.m` loads the completed comparison.mat without recomputing
or altering predictions. It exports an unfiltered center-probe overlay plus
residual, and a five-probe detail figure over the fixed first five scoring
periods (10--15Tp, 137.62--206.43 s). Files `gl_hos_center_detail.png`,
`gl_hos_five_probe_detail.png`, `gl_hos_detail.json` and the center detail CSV/PDF
are retained locally beside the full figure. The center error over this short
display window is 8.883%, while the originally specified 10--70Tp scoring
window gives 4.2147%; do not confuse the two. The comparison is specifically
the reconstructed second-harmonic phase sector, not total free-surface elevation.

The user requested postprocessing of the completed four-phase GPU HOS run.
Source: `/root/hos-low-80tp-20260927T052004Z` on 92, four complete 5505-row
probe files ending at 1100.8 s. Each source file was verified against the
completed controller's SHA-256 before conversion to CSV. Only compact probe
records were transferred over an IP-restricted, temporary LAN endpoint; full
modal fields remain on 92.

Active processing directory on 93:
`/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-low-postprocess-20260927T085331Z`.
Local deployment evidence: `artifacts/hos_gpu/postprocess-20260927T085331Z/`.
Read `status.json` and `matlab.log` there before launching anything; do not
duplicate an active job. The detached controller runs as lxy, MATLAB is pinned
to CPU40 and uses one computational thread. The first root-account launch
failed at MATLAB sign-in; its logs/status remain as `*.root-account-failed`.
No scientific calculation from that launch was used. Later deployments use
the original licensed account directly.

Scientific code is an unchanged copy of the original campaign's
`compare_jonswap_time.m`, SHA-256
`88023c5af9a97cd5529acfdc614bc6d4df0ef5b9392d21b338e7d62e4d237dc9`.
Its archived GL dependencies and the original initial C/kx/ky/omega arrays
are reused. The original CPU campaign and its postprocessing are untouched.
All source and transfer hashes are in remote provenance.json/snapshot.json.

The procedure retains the four-phase/Hilbert first-harmonic input, the
second-harmonic phase-sector reference, original one-third off-axis amplitude
eligibility, original frequency-support rules, GL16, direction bins of
3.75 and 1.875 degrees, and fixed 10--70Tp scoring window. It also reports
unfiltered errors, common-sum-band errors, initial-spectrum-only and
direction-blind baselines, and the 10 s record truncation sensitivity.
No gain, offset, phase or time shift is fitted. Phase sectors are not exact
perturbation orders. The directional distribution is supplied by the initial
spatial spectrum, not inferred uniquely from a single probe.

Outputs under `gl-time-comparison-v1`: metrics.csv, short_center_metrics.csv,
report.json, comparison.mat, full_comparison.png and full_comparison.pdf.
Retrieve compact metrics/report/plots only; comparison.mat stays remote.
Controller completion means these outputs exist and the MATLAB report says
FULL_RECORD_COMPARISON_COMPLETED. It does not imply small scientific errors:
inspect the metrics and rendered plot before drawing conclusions.

92 also has MATLAB R2025b and ample free RAM, as verified during this task.
93 was chosen to reuse the established R2026a environment, source and prior
inputs. This MATLAB code currently runs on the CPU; moving it to the GPU host
alone would not enable GPU execution. No hardware settings are modified.
