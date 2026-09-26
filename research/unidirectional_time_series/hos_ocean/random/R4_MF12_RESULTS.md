# 11128-mode direct R4/MF12 result

Executed 2026-09-26, CPU40 only, sequential fresh MATLAB processes. Main run
finished19:03:27 UTC; diagnostic capture completed afterwards. Same11128
JONSWAP parents, same phases, same1024x512 grid, kp Hs/2=.12, Hs=8.6021505376m.
No HOS launch, no new GL initialization, no fitting or output filter.

## Measured costs

| Work | Operator seconds | Cold process wall (s) | Process peak RSS (GiB) |
| --- | ---: | ---: | ---: |
| R4 eta20 through original rejection | 12.181351 | 24.00 | 1.5522 |
| R4 psi20 through original rejection | 6.080959 | 16.27 | 1.6422 |
| MF12 order2 coefficients + one20 eta/psi surface | 230.8240 | 243.21 | 67.4119 |

MF12 coefficient generation52.2260s, sector selection.3452s, spectral surface
178.2528s. The unmodified coefficient API includes22 coefficients too. R4 is
20-only; these are implementation-path comparisons, not identical kernel-work
counts. R4 combined compute18.2623s versus230.8240s (~12.64x), and two cold R4
processes40.27s versus one cold MF12 process243.21s (~6.04x). GNU time -v reports
the process peaks above; controller process-tree sampled peaks are saved too.
This is one difference snapshot, not four-phase complete HOS initialization.

## Error from captured R4 fields

Both unmodified R4 operators computed their fields but then rejected them at
the unchanged unrepaired-region Hermitian-defect threshold1e-8:

- eta20:3.3320491e-8, spark:Eta20RSeriesStability.
- psi20:1.7856625e-8, spark:Psi20RSeriesStability.

The benchmark therefore preserves rejected_by_original_operator in the
primary reports. To obtain the requested large-input accuracy comparison,
separate diagnostic copies save intermediate fields immediately BEFORE the
same original rejection. Formulas, repair cutoff, threshold and rejection
remain unchanged. Only the save operation was inserted; original and
instrumented SHA256 hashes are retained in diagnostic-source.json.

| Field | Raw full-grid relative L2 | Relative Linf | Absolute RMS | R4/reference norm |
| --- | ---: | ---: | ---: | ---: |
| eta20 | 1.879125% | 1.290832% | .002524642 m | .98736264 |
| psi20 | 1.773910% | 1.994692% | .273513441 m2/s | .98967899 |

These characterize the finite R4 fields at the rejection point, not accepted
production output. The approximately1.8% model discrepancy and approximately
1e-8 numerical consistency indicator are different quantities. A small check
failure is not evidence that the physical field error is catastrophically large;
conversely this comparison does not silently authorize bypassing the check.

Each R4 low-mode repair treated30 output cells,163021 oriented contributions,
with measured diagnostic repair times.279s eta and.187s psi. This is targeted
low-output work, not enumeration of the full123821256 ordered cross pairs.

## Interpretation and next boundary

The large-support experiment supports R4 as a practical20 candidate: materially
less time/memory and approximately2% raw field error on this specific prescribed
spectrum. It does not yet establish a drop-in accepted HOS initializer: the
automatic interface exits at its existing stability check, which needs an
explicit numerical treatment before integration. No threshold was relaxed,
no Stokes correction was added, and no public SPARK source was modified.
Independent MF12 reference is retained for further checks. Low-amplitude
homogeneity is not counted as a separately executed second-amplitude test.

Evidence: remote r4-mf12-jonswap11128-20260926-v1; compact local reports and plot
in artifacts/hos_ocean/r4-mf12-11128/. Raw field MAT files remain remote.
The source benchmark design and timing scope are in R4_MF12_BENCHMARK.md.
