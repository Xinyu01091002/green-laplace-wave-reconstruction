# Off-centerline comparisons with the user's one-third amplitude criterion

The user clarified that extremely weak far-lateral signals are not the
primary application: the eta22 time-record wave amplitude should be at least
one third of the centerline value. This is a new user-defined signal-strength
eligibility criterion, introduced after the earlier stress-test results.
It does not retroactively turn those stress tests into successful results.
All previous outputs remain available and no GL coefficient is fitted.

## Definition and steepness

Primary quantity is the maximum absolute *sampled OW3D reference eta22* over
the full saved record. All lateral comparisons here use x=4500 m, with the
same-x centerline at y=3375 m. The centerline peak is 0.01629001170036 m, so
the threshold is 0.00543000390012 m. The full-record peak-to-trough range is
also reported (centerline 0.03033928369611 m), as is the main-group peak ratio.
The lateral points have identical pass/fail classification under peak and
whole-record peak-to-trough definitions. No sub-sample maximum is inferred.

All these near/far/multiple-probe checks use nominal initial Akp=0.02,
kh=1 and direction-spread label 25 degrees. The separate parameter sweep
also contains Akp=0.12, but that larger steepness was not used for these
off-centerline checks. Akp is the source-case label, not a newly fitted local
steepness. Each record still has 4 s sampling.

On-centerline x-offset probes are reported separately. Their same-x centerline
comparator is themselves; the CSV additionally records their ratio to the
focus-center probe for context. This is not a global-focus-amplitude gate
applied indiscriminately to all centerline positions.

## Recheck of previous lateral probes

At the old nominal +/-0.5 lambda_p locations (actual dy=+/-105.46875 m),
eta22 peak is about 0.0037353 m, or 22.930% of centerline: below one third.
At nominal +/-1 and +/-1.5 wavelengths, full-record peak ratios are about
0.040%--0.043%. The +/-1.5-wavelength main-group ratios are smaller still.
Thus the large relative errors at the far lateral probes lie outside the
newly specified useful-signal region. Even the previous half-wavelength
lateral points, despite their sub-percent reconstruction errors, are excluded.

## Additional fixed lateral probes

Before computing their GL errors, specify offsets +/-0.25 and +/-0.35 lambda_p,
snap to the nearest grid point, and check the amplitude criterion. All four
pass. The algorithm, initial spectrum, angular grids and main-group window
definition are unchanged. No reference-based gain, phase or shift is applied.

| Point | Actual dy (m) | eta22 peak (m) | Peak ratio | Range ratio | eta22 L2 (%) | psi22 L2 (%) |
|---|---:|---:|---:|---:|---:|---:|
| y_m0p25 | -52.734375 | 0.01141364 | 0.70065 | 0.70459 | 0.26543 | 0.22015 |
| y_p0p25 | +52.734375 | 0.01141370 | 0.70066 | 0.70459 | 0.23861 | 0.18916 |
| y_m0p35 | -79.1015625 | 0.00723627 | 0.44422 | 0.44985 | 0.33464 | 0.39411 |
| y_p0p35 | +79.1015625 | 0.00723634 | 0.44422 | 0.44985 | 0.28029 | 0.38513 |

Errors use the fixed main-group window and 7.5-degree joint-input direction
grid. The unchanged center benchmark is eta22 0.26681%, psi22 0.14016%.
The four eligible lateral points are below 0.4% for both variables on these
saved main-group samples. This supports the restricted weak-steepness signal
region, not arbitrary weak-signal probes or larger-steepness validity.

## Artifacts and command

```matlab
addpath('research/directional_wave_data');
run_directional_amplitude_gate('C:/Research/VWA/VWA Unidirectinal/Directional/test1');
```

`results/directional_amplitude_gate/` contains old-point reclassification,
full criteria, new raw extractions/source hashes, per-probe results, summary
CSV/MAT and the amplitude-qualified PNG/PDF figure. Center raw records and
the initial spectrum are checked against the prior extraction before its
convention audit is reused. No previous result is overwritten.
Both new MATLAB functions pass Code Analyzer. The threshold uses the OW3D
reference rather than GL output, and neither pass/fail nor the reported
relative errors rescale any waveform.
