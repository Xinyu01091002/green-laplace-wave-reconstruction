# Partial directional HOS / GL comparison, 2026-09-26

This is a preview from an immutable common prefix of the ongoing production
run, not the final 220 s result. No HOS process was paused, modified or rerun.
The snapshot was taken at 06:13:49 UTC, with 869 common samples spanning
0--173.6 s at .2 s intervals, for all four phases and five probes.

Remote preview:
`hos-directional-fourphase-20260926-v2/partial-comparison-20260926-v1`.
Source runs and native probe records remain remote. snapshot.json records
common-prefix hashes and the original available sample counts.

## Reconstruction

The archived existing allocate_directional_record and gl_directional_sum_time
implementations were used, with GL quadrature J=8 and 15/7.5 degree angular
bins. Directional complex shape is supplied by the frozen first-order input
spectrum, then constrained by the observed first-harmonic time record at each
probe. Higher-order reference records are not used in allocation. The
negative temporal exponent convention is retained. This is not direction
inference from a single time record without prior information.

Temporal Hilbert-four-phase separation is applied along time. Parent bins
are restricted, before scoring, to the original eligible first-order frequency
range, kh>=.3 and resolvable frequency sums. Actual retained parent angular
frequencies are .18075907098--.65073265552 rad/s (bins 5--18). The same parent
band is applied to the observed first harmonic and its directional prior.
The second-harmonic sum-support projection is applied identically to target
and every prediction. Raw output comparisons are also retained; there is no
fitted gain, phase, offset or time shift.

## Preview results

Main-group relative L2 percentages for the 7.5-degree directional allocation,
with the identical output projection on GL and HOS:

| Actual transverse offset | eta22 peak ratio to center | GL eta22 relative L2 |
| --- | ---: | ---: |
| 0 m | 1 | 4.22797% |
| -52.7821 m | .709279 | 3.10274% |
| +52.7821 m | .709279 | 3.14256% |
| -70.3762 m | .541117 | 2.58541% |
| +70.3762 m | .541117 | 2.62474% |

All five points pass the user's one-third amplitude criterion on this common
record. Main-group windows, selected only from each first-harmonic envelope,
are complete: center 81.876--136.924 s and lateral 82.076--137.124 s.
The center raw-output error is 4.22838%, so the common output filter does not
materially create the observed agreement. The source-based parent projection
relative L2 is .1456% at center and .0577--.0886% laterally.

Angular refinement from 15 to 7.5 degrees changes the main-group GL prediction
by .398% at center and .434--.513% laterally. Observed-energy-weighted
allocation conditioning is about 1.00045 at center, 1.19276 at the nearer
lateral pair, and 1.36940 at the farther pair. No division floor or fitted
regularization was used.

## Partial-record sensitivity

A fixed 10 s shorter center prefix (0--163.6 s) was processed with the same
rules. Within the original main-group window:

- first-harmonic input difference: .02891% relative L2;
- GL 7.5-degree prediction difference: .08291% relative L2;
- common-band reference difference: .002922% relative L2.

The corresponding center reconstruction metric is 4.21392% versus 4.22797%
in the longer record, a .01405 percentage-point change. The sensitivity
includes Hilbert endpoints, Fourier-grid changes and the input-derived
frequency support, not Hilbert effects alone. It supports viewing this as a
useful preview but does not replace the final full-record comparison.

## Resources and artifacts

One MATLAB process, constrained to CPU 32 with -singleCompThread, ran beside
the existing jobs. Entire postprocessing including startup, sensitivity check
and plotting: 113.40 s wall, peak RSS 1587820 KiB (about 1.51 GiB). The five
main probe calculations together took 5.578 s by internal MATLAB timers;
the shorter-window check and plotting are outside those internal timers.
No new HOS solve was launched.

Local artifacts: `artifacts/hos_ocean/directional-partial/` contains
partial_comparison.png/.pdf, metrics.csv and report.json. The original 220 s
production remains independent of this preview. Future full-record processing
must preserve this snapshot and identify changed windows/bins explicitly.