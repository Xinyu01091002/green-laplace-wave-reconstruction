# Energy-ranked free-mode support count

User corrected the targets to kp Hs/2=.02/.12: Hs=1.4336917563/8.6021505376 m.
The existing linear preview has half these amplitudes and is not a simulation.
Uniform amplitude scaling leaves the counts below unchanged.

Actual Cartesian design: gamma3.3, .5--2.5 fp, Gaussian energy directional
sigma17.6777 deg on +/-60 deg, domain50x20 lambda, grid1024x512. Sort each
retained free mode by variance contribution |C_j|^2/2, largest first, and retain
the minimum number reaching the specified fraction. This is relative to the
already frequency-truncated design, not the infinite JONSWAP spectrum.

| Retained variance | Modes | Ordered cross pairs | One double pair array (GiB) |
| --- | ---: | ---: | ---: |
| 100% | 23526 | 553449150 | 4.1235 |
| 99% | 11128 | 123821256 | .9225 |
| 99.5% | 13420 | 180082980 | 1.3417 |
| 99.9% | 17415 | 303264810 | 2.2595 |

99% retains47.30% of modes and22.37% of ordered cross pairs, cutting pair
storage/work proxy by77.63% (about4.47x). This is not a measured MF12 speedup;
the initializer allocates many arrays and intermediates, not one. No MF12
coefficient generation or HOS propagation was performed for this count.

Actual retained energy .99000010236. Normalize selected amplitudes by
1.0050377633 to preserve the target Hs; preserve original phases and apply the
same retained set at both steepnesses. Suggested consistent initialization:
use the same reduced first-order spectrum for MF12, HOS initial eta/psi and
GL initial directional prior. Keeping all linear parents while providing only
their truncated quadratic initialization is a different transient experiment.
Later nonlinear HOS-generated modes must not be confused with this initial
free-parent truncation.

Retained relative frequency bandwidth=.2527422 versus full=.2661706. Direction
RMS=17.3838 deg. The2--2.5 fp band retains73.61% of its own original variance:
global99% does not imply every frequency band retains99%, nor does it bound
nonlinear reconstruction error by1%. Thus99% is a plausible practical pilot,
not an automatic nonlinear-accuracy certificate.99.5/99.9% can be compared in
initialization/postprocessing if needed before adding full campaigns.

99% support second-order axis bounds are9.2 kp in x and3.2 kp in y. Thus the
original1024x256 grid (Nyquists10.24/6.4 kp) can represent this reduced initial
second-order support. That reopens a cheaper grid candidate, but does not
establish HOS propagation convergence or full third-order support (x bound
13.8 kp exceeds10.24). Do not silently halve the published design grid without
an explicit updated numerical design and short broadband check.

Evidence: remote jonswap-energy-support-v1 under the existing random-run root;
local artifacts/hos_ocean/jonswap-energy-support/{counts.csv,report.json}.
MATLAB counting run completed in12.62 s; no nonlinear solver was invoked.
