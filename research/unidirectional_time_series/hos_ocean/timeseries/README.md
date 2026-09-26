# HOS GL single-point time-series check

Completed for both steepnesses. See [RESULTS.md](RESULTS.md) for errors,
measured runtime/memory, and the supported interpretation.

This follows the user's clarified scope: near-matched OW3D physical settings
are sufficient; do not pursue exact solver-to-solver parity or broad accuracy
sweeps. Reuse the established temporal reconstruction workflow and compare
predictions/reference with identical output filtering.

## Setup

- kph=1, g=9.81, kp=.0279; Alpha=1 semi-Gaussian spectrum; Akp=.02 then .12.
- Same 68 lambda_p nominal domain and .67L focus location as the boundary-
  corrected OW3D family; 4096 unique periodic HOS points versus 4097 wall-
  bounded OW3D nodes. Start 40Tp before nominal focus, end 10Tp after focus.
- Keep the user's no31 choice: independent MF12 11+20+22+33 eta and true
  surface psi at four phases, muStar=0, linear frequencies. Legacy OW3D
  enabled primary third-order corrections; this difference is intentional.
- Native HOS probe at index 2745, approximately x=10259.09 m; legacy probe
  x=10259.18359375 m. No fitted spatial/time shift is applied.
- M=5, qx=5, Ta=0, no breaking/dissipation; relative integration tolerance
  1e-10. One fixed production setting, not another multi-level sweep.
- Save eta at Tp/40. Select 32--48Tp elapsed (focus +/-8Tp), 641 samples.
  Harmonic separation is the verified temporal Hilbert-four-phase formula.

## Preserved reconstruction

The actual committed gl_eta22_time_pairs and mf12_eta22_time implementations
are archived, along with their required JSON formula provenance. GL ranks
6/8/12, MF12, VWA and Walker use identical first-harmonic inputs. The parent
rule is kh>=.3 and twice the parent frequency below the temporal Nyquist,
matching run_ow3d_eta22_pilot.m. No energy threshold is tuned to the target.

Primary raw comparisons use the identity output filter for both prediction
and target, as did the saved OW3D pilot. Supplementary common_sum_band results
apply the exact same native FFT mask to the reference and every prediction:
positive/negative output bins between twice the minimum and maximum retained
parent bin. This is fixed by input support, not reference error. No distinct
bandpass or phase convention is applied to one side only.

Both HOS records and the already saved OW3D records pass through the same
processing function. Recomputed OW3D predictions must match saved predictions
to relative Frobenius norm <1e-10. Full-window and input-envelope +/-2Tp main
window errors are recorded. No gain, phase, offset, or time alignment is fitted.

## Execution and IO

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1`.

The tested v2.1.0 IO patch was extended by one formatting line for probe
output, ES13.5 -> ES25.16E3. No numerical equation was edited. Binary SHA-256:
`13c4bae40855f3d3a29da12d5567559c4a41e415c27bc5f67d5b927db6c2f060`.

The probe-only output avoids repeatedly saving full spatial fields. A short
smoke case wrote both native fields and probe data. Its initial eta values
agreed with each other and the prescribed input within 1.2e-15 m. The initial
reader assumed a ZONE record, but probes.dat starts data after VARIABLES;
this parsing issue was corrected before production, without rerunning the
smoke solver. The original failed controller record is preserved.

The first/higher fields are finite and t=0 phase recovery is checked during
preparation. Every probe sample and requested time is checked during analysis.
Per-process GNU time records wall, user/system CPU and peak RSS for HOS and
MATLAB stages. Raw fields/probes and comparison MAT files remain remote.

Periodic HOS and wall-bounded OW3D are not identical boundary problems.
The nominal packet centers have substantial boundary clearance in the
analysis window, but weak/free-wave boundary influence has not been ruled
out separately. The aim is a practical independent GL reconstruction check,
not certification of solver equivalence or all wave conditions.