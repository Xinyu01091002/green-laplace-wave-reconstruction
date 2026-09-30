# Single-direction random-wave modal comparison

This comparison uses the frozen three-seed random-phase HOS campaign at
`k_p h=1`. All four HOS phase runs were tapered over `5 T_p` at each record
edge. The first-order input was reconstructed over the complete `50 T_p`
record and `eta22/eta33` were scored over `15--35 T_p` without alignment or
amplitude fitting.

The compact HOS source file has SHA-256
`8195c23027e7ad3d2e9f939368abdc69bbd0ad41be061a8f54278c89672e21eb`.
No HOS simulation was rerun.

## Three-seed fixed-input comparison

The table uses nominal `Ak_p=.12`, GL8 and the same tapered `eta1(t)`. Modal
timings use 8192 points. Auxiliary timings are local Python one-shot calls of
the final `L/h=1600`, `Nx=32768` level on the same machine.

| Seed | Modal time | Auxiliary time | Ratio | Modal eta22/HOS | Auxiliary eta22/HOS | Modal eta33/HOS | Auxiliary eta33/HOS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20260925 | 6.91 s | 164.58 s | 23.8 | 0.159% | 0.160% | 52.473% | 52.474% |
| 20260926 | 6.78 s | 88.22 s | 13.0 | 0.253% | 0.252% | 39.637% | 39.636% |
| 20260927 | 6.82 s | 124.50 s | 18.2 | 0.174% | 0.173% | 71.563% | 71.563% |

For 4096 to 8192 modal points, the three-seed changes were
`0.021--0.022%` for `eta22` and `0.037--0.041%` for `eta33`. The saved
auxiliary-domain last-level changes were `0.361--0.476%`.

## Full random-wave matrix

The automatic modal ladder was run on all 27 combinations of three seeds and
nominal `Ak_p=.02,.04,...,.18`. Every case passed the `0.1%` change gate at
8192 points after evaluating 2048, 4096 and 8192 points.

| Quantity | Result |
|---|---:|
| Modal convergence time, median | 12.07 s |
| Modal convergence time, range | 11.96--12.33 s |
| Saved auxiliary time, median | 126.93 s |
| Saved auxiliary time, range | 125.22--130.76 s |
| End-to-end time ratio, median | 10.55 |
| Modal final change, median | 0.0458% |
| Auxiliary final change, median | 0.405% |
| Modal/auxiliary waveform difference, median | 0.765% |
| Modal HOS eta33 error, median | 52.517% |
| Auxiliary HOS eta33 error, median | 52.531% |
| Modal q-projection RMS, median | 0.307% |

The modal route preserves the earlier physical conclusion: after `Ak_p=.04`,
the eta33 discrepancy is nearly amplitude-independent and strongly dependent
on random-phase realization. The three seed families remain approximately
`39--40%`, `52--53%`, and `71.5%`.

The complete 27-case metrics are stored in
[`python/reference/random_hos_modal_ladder_metrics.csv`](../reference/random_hos_modal_ladder_metrics.csv)
with SHA-256
`b00d75d2bdcf3d4623dc3b9f837a6e548c1041f4b15479da4c0f719fbe03b05a`.

## Decision

Within the maintained scope, the modal route has the same HOS accuracy, a
smaller numerical refinement change, and substantially lower runtime. The
single-direction public time API therefore uses the one-dimensional modal
resolvent. The auxiliary-domain executor is removed from the active package.

Directional spatial reconstruction is unchanged and continues to use the
declared physical `x,y` FFT grid.
