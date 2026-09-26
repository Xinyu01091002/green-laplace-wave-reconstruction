# Directional JONSWAP random-wave design (not launched)

User request: broader random sea, designed from an energy spectrum, with
**kp Hs / 2 = .02 and .12**, corrected by the user after the initial design preview. These are not potential focusing Akp.
This design supersedes the phase-only amplitude rule for the proposed new
family. The old HOS runs remain unchanged. Only a linear spectral design and
preview were computed; no MF12 coefficients or new HOS run were started.

## Physical specification

- Keep g=9.81, kp=.0279 /m and kph=1, hence h=35.84229391 m,
  lambda_p=225.20377445 m and Tp=13.76199709 s, fp=.07266387238 Hz.
- Target spectral significant height Hm0=4 sqrt(m0): low 1.4336917563 m;
  high 8.6021505376 m. Report empirical time-domain H1/3 separately if used.
- Use a prescribed JONSWAP frequency **variance density**, gamma=3.3,
  sigma=.07 below fp and .09 above. At this finite depth this is a specified
  test spectrum with finite-depth dispersion, not a claim of an equilibrium
  sea or a TMA spectrum. No TMA shape factor is implicitly added.
- Positive-frequency shape proportional to
  f^-5 exp[-1.25(fp/f)^4] gamma^exp[-(f-fp)^2/(2 sigma^2 fp^2)].
  Normalize the retained, discretized spectrum to the target Hm0 rather than
  impose a wind/fetch-dependent alpha and an independent inconsistent Hs.
- Retain .5 fp <= f <= 2.5 fp. Continuous gamma3.3 variance loss relative
  to the broad untruncated shape is 2.0654%; normalize retained variance.
  Declare this explicitly as a truncated JONSWAP realization.
- Keep comparable angular spread: energy D(theta) is a Gaussian with
  sigma_theta=25/sqrt(2)=17.6777 degrees, mean0, normalized on [-60,60] deg.
  The old 25-degree parameter described amplitude, not energy. Choosing
  25 degrees for the new energy Gaussian would broaden directions too.
- Fixed seed20260925. Use identical modal phases for high/low. Four copies
  of one realization have global shifts0/90/180/270; never four new seeds.
  No spatial taper, focusing condition or selection of a favorable seed.

## Cartesian modal construction

For native HOS periodic modes, omega(k)^2=g k tanh(kh), f=omega/(2pi).
Cell variance v_j = S(f_j) D(theta_j) [cg_j/(2pi k_j)] Delta kx Delta ky.
The Jacobian is df dtheta = cg/(2pi k) dkx dky. Normalize sum(v_j)=(Hs/4)^2,
then C_j=sqrt(2 v_j) exp(i phi_j) in eta1=Re sum C_j exp(i k.x-i omega t).
Thus variance is sum |C|^2/2; do not use S itself as modal amplitude.

The earlier linear preview used the superseded kp Hs=.12 normalization and verifies full-domain discrete Hs=4.3010752688 m
to relative error <1e-12. Its elevations must be doubled for the corrected high target; spectral counts and normalized bandwidth do not change. Discrete sigma_f/mean_f=.26617064, versus .09357931
for the old phase-only random field. The new low first-order field is high/6;
nonlinear eta and true surface psi must be recomputed or scaled by order,
never divide the full nonlinear initial field by six.

## Spectrum alternatives measured before scoring

| gamma | upper f/fp | omitted variance (%) | sigma_f/mean_f |
| --- | ---: | ---: | ---: |
| 3.3 | 2 | 4.9285 | .22219 |
| 3.3 | 2.5 | 2.0654 | .26616 |
| 3.3 | 3 | 1.0041 | .29512 |
| 1 | 2.5 | 3.1496 | .28933 |

Recommend gamma3.3 first: it is already substantially broader than the prior
case. Gamma1 is the PM limit and is an optional later bandwidth variation,
not another run silently added to the two-amplitude matrix.
Formula/conventions reference:
https://wavespectra.readthedocs.io/en/latest/construction.html
https://wavespectra.readthedocs.io/en/latest/generated/wavespectra.core.npstats.jonswap.html

## Numerical proposal and limits

- Periodic domain50x20 lambda_p, matching the earlier domain.
- Candidate grid1024x512, rather than1024x256: new first-parent kmax=4.7606 kp.
  The retained parent support has second-order axis bounds |Kx|<=9.52 kp,
  |Ky|<=8.20 kp, whereas old Nyquist_y was6.4 kp. New Nyquists10.24/12.8 kp
  represent these sums. The MATLAB design measures these bounds.
- This is a candidate for the second-harmonic experiment, not a spatial
  convergence claim. Complete third-order x support extends to14.28 kp,
  beyond1024's Nyquist_x. If full third-harmonic support is required, use
 2048x512 (Nyquists20.48/12.8 kp) or redesign support. Do not claim full
  third-harmonic validation on the candidate1024x512 grid.
- Reuse adaptive Cash-Karp5(4), M5, qx=qy3 and abs tolerance1e-12 as initial
  candidate settings, subject to a short broadband high-case check. q3 remains
  partial dealiasing at M5. No claim that a narrowband benchmark certifies
  broadband high-steepness behavior. Ta0 with MF12 order2 eta/true-surface psi;
  no initial31/33, no fitted damping/repair.
- Proposed record80 Tp=1100.959767 s. Score10--70 Tp, plot the entire record.
  Treat these margins as a planned diagnostic, not proof of transient removal.
  Output every.2 s and preserve actual timestamps; 80Tp is not an integer
  multiple of.2, so do not reuse the hard-coded1101-record completion check.
- Retain the same five physical probe coordinates. Add sparse full-domain
  eta/psi snapshots at planned times for spatial checks; use nearest actual
  output times and record them rather than claiming exact Tp multiples.
- Low then high, one family at a time; each family four phases x eight MPI
  ranks. One seed is an initial controlled comparison, not ensemble statistics.
- A periodic random field needs no artificial absorbing strip. Long records
  in a finite periodic domain are not independent infinite-ocean realizations;
  domain/seed dependence remains a separate question.

## Resource and implementation conditions before launching

The unchanged grid-area scaling gives about8.5 GiB sampled aggregate HOS RSS
for four eight-rank cases; use a provisional16 GiB budget and verify in a
short run. If2048x512 is selected, about17 GiB and a32 GiB budget are geometric
estimates, not measurements. Old32 compute ranks plus retained8 OW3D fits
the previously observed48 CPU slots, but recheck live resources at launch.

At1024x512 and80Tp, simple grid-area x duration scaling of the earlier random
run gives13.72 hours per family. This omits FFT scaling, high-frequency time-step
restriction and the higher actual sea-state amplitude. It is a planning scale,
not a benchmark, upper bound or finish-time promise. A short run must establish
runtime before full execution; no duplicated full-length MPI benchmarks.

There are23526 parents and553449150 ordered cross pairs. A single real8-byte
array indexed by all pairs is4.124 GiB; MF12 allocates multiple such arrays.
Do not invoke the current all-pairs initializer blindly. Prepare pairs in
blocks using the established MF12 formulas, accumulating onto the native grid,
or separately design and quantify a support reduction. Never silently prune
the spectrum merely to fit memory. First-order design itself needs no pairs.

GL postprocessing also needs attention:80Tp gives roughly160 temporal parent
bins across the declared band. At1.875 degrees the existing dense pair path
can have order10^4 allocated parents, many more than the prior17-bin case.
Use sequential probes and bounded pair processing; retain directional
conditioning, record-window sensitivity and a small J-quadrature consistency
check appropriate to the wider frequency support. Do not assume J8 narrowband
accuracy transfers without checking. No new GL kernel or empirical correction
is proposed here.

## Design artifacts

Remote `jonswap-design-v1` under the established hos-random-fourphase-20260926-v1
root holds the linear_design.mat and report. Local compact artifacts:
artifacts/hos_ocean/jonswap-design/{report.json,spectrum_options.csv,jonswap_design.png}.
The plot is a linear design preview, not HOS propagation or an MF12 validation.
