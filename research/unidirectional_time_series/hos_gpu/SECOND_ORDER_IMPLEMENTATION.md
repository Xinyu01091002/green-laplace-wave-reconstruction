# Actual second-order time-series implementation and route removal

The user requested deletion of the added interpolation/low-rank route and an
explanation of the earlier second-order calculation. This records actual
calls rather than claiming the native spatial FFT-GL was used end to end.

## What the second-harmonic comparison actually did

1. Obtain the first-harmonic time record from HOS four-phase/Hilbert processing.
   The observed second-harmonic record is retained solely as a comparison.
2. Fourier-project that first record onto the original eligible initial-spectrum
   frequency range. For the center probe, 149 positive bins span
   .273926--1.118533 rad/s. They retain 99.8065% of positive-frequency energy;
   waveform relative L2 change is 4.3999%. This is an input restriction, not
   use of the untouched full first-harmonic record. The later 4omega_p trial
   concerned third order and did not retroactively change these second-order
   results.
3. At each retained frequency, infer k from linear dispersion. Propagate the
   INITIAL complex spatial spectrum to the probe analytically, group its
   contributions by direction, and obtain a finite-record temporal directional
   prior. Use `allocate_directional_record.m` to set
   A(f,d) = observed(f) * prior(f,d) / sum_d prior(f,d).
   This preserves each prior complex directional shape and matches the observed
   frequency coefficient. It is an assumption, not unique direction inference.
4. Call the existing `research/directional_wave_data/gl_directional_sum_time.m`
   with all nonzero frequency-direction components. That function explicitly
   forms pair matrices of output wavevector, sum frequency and forcing terms.
   It evaluates the GL quadrature response for EACH PAIR (eta uses 16 nodes in
   the completed comparison), then multiplies by A_i*A_j.
5. Aggregate pair coefficients at b_i+b_j with `accumarray`, then call `fft`
   to synthesize the negative-time-exponential convention. The result is
   Re sum C_l exp(-i*omega_l*t). The FFT here accelerates final synthesis;
   it does NOT replace the preceding quadratic pair computation.
6. Compare raw and common-output-band predictions to HOS, with no fitted gain,
   offset, phase or time shift. The fixed scoring window is 10--70Tp. At
   1.875-degree resolution there are 9536 active frequency-direction components,
   hence roughly 90.9 million ordered pairs.

Thus this was a pairwise temporal GL-kernel reference implementation, not a
direct call to the original fast native spatial FFT-GL operator. It did not use
the later interpolation, Chebyshev, TT-SVD or operator-cache route. The reported
second-harmonic errors apply to this particular projected input and directional
prior, not to an unmodified full-record original-GL execution.

The second-order DIFFERENCE trial uses the user's requested pure R4 kernel,
also via its explicitly validated temporal pair adapter. It does not substitute
GL for R. Native spatial R and original GL sources remain unchanged.

## Why the original spatial call was not used

The original GL interface can directly sample bound-wave fields from a known
native spatial first-order spectrum, applying exp(-i*omega*t) to its modes.
The present time-series reconstruction instead constrains the input by a
measured/separated point record and allocates frequency-direction coefficients.
Their dispersion wavevectors generally are not native spatial FFT lattice
points. That input representation requires an explicit, justified adapter.
The assistant chose a pairwise reference and later a triple reference instead
of resolving that native interface. The subsequent numerical low-rank route
was an additional implementation choice, not a requirement of GL, and must
not be presented as the user's original algorithm.

## Completed deletion

Removed the added interpolation/TT-SVD implementation, plan apply/load helpers,
its deployment/controller/test/plot entry points and active method document.
Removed local experiment source archives and deployed experimental source
copies/archives on 93. Removed the 204-file cached operator (634882168 bytes,
about 606 MiB). No active route processes remained at deletion.

Original GL code, scalar GL quadrature helper, HOS CUDA solver and cuFFT code,
R4 implementation, HOS initial conditions and completed simulation records
remain. Numerical reports/figures/logs from the withdrawn experiment remain
only as historical evidence, with a removal manifest. Do not advertise their
timings as original-GL timings or reconstruct this route from scratch unless
the user explicitly requests it.

Removal evidence: local
`artifacts/hos_gpu/fft-timeseries-20260927T155457Z/removal-manifest.json`
and `artifacts/hos_gpu/route-removal-20260927/remote-result.json`; remote
`gpu-fft-timeseries-20260927T155457Z/removal-manifest.json`.

## Earlier OW3D time-series comparison, before HOS

The earlier unidirectional OW3D driver `run_ow3d_eta22_pilot.m` calls
`gl_eta22_time_pairs.m` at GL6/8/12 (6 was the declared primary choice).
The latter transfers the released pure GL ordered-pair kernel, with
Gauss--Laguerre integration, no Stokes/angle correction, and explicit pair
matrices and time exponentials. It is not the native spatial FFT executor.
The tracked implementation history identifies commit a8aa534 for this helper.
Its input rule was kh>=.3 and 2omega below temporal Nyquist; it did NOT apply
the later JONSWAP initial-frequency min/max restriction or energy-ranked cap.

The directional OW3D driver `run_directional_joint_pilot.m` calls
`gl_directional_sum_time.m` at rank 8. This explicitly computes pair kernels,
then aggregates sum-frequency bins and synthesizes with FFT. Its helper's
tracked history identifies a3162d5. Later HOS/JONSWAP calls used rank 16;
do not retroactively label the OW3D run GL16.

The earlier eta33 OW3D trial calls `gl_eta33_time_triples.m` at ranks 4/6/8,
an explicit ordered-triple reference executor. Its historical eta20 trial
used the shared-scale GL difference diagnostic, not the later requested R4.

`test_time_pair_reference.m` checks the temporal pair evaluator against
`gl_spectral_surface` and the frozen ordered GL8 reference on small fixtures.
Thus the original spatial GL WAS used for implementation consistency tests,
but the OW3D time-series driver itself did not directly call the fast spatial
FFT-GL as its reconstruction executor. These historical routines predate the
now-deleted interpolation/TT-SVD route; they do not contain that added method.
