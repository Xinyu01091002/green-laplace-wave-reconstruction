# First directional joint-input experiment

This is the first executed test of the user-proposed combination: initial
linear spatial spectral structure plus the observed probe first-order time
record. It is a prescribed conditional-amplitude model, not a unique inverse
recovery of an evolved directional spectrum.

## Data, geometry and information boundary

Read-only OW3D source is
`C:/Research/VWA/VWA Unidirectinal/Directional/test1/`
`kd1.0_spread_25_Akp_0.02_phi_shift_{0,90,180,270}`.
All 364 EP files at steps 0:10:900 were read and hashed. The saved interval is
4 s and time coverage is 0--360 s. h=35.84229 m, g=9.81 m/s^2, kp=0.0279 rad/m.
The probe is the declared domain midpoint (4500,3375) m, selected without
looking at high-order errors. Historical Alpha is not certified here.

The EP file has extra coordinate points. Select physical coordinates
0<=x<Lx and 0<=y<Ly, obtaining the 1024 x 256 FFT grid on 9000 x 6750 m.
This is an explicit periodic extension of the selected physical field, not
a claim that the closed-tank boundary is periodic. At the initial time, the
reported edge/field norm ratio is 6.1625e-5.

At every saved time, form the complex four-phase first class and retain its
positive-kx part to remove the conjugate third-harmonic alias on strict-forward
support. Its real part at the probe is the observed input eta1(t). Current
spatial fields supply no direction weights to the model. Only the t=0 first
class supplies the initial spatial spectrum. Both are labelled first-phase-
sector approximations, not certified exact perturbation-order eta11.
The high-order reference is the raw real second phase combination at the
probe, for eta and psi. It is used only after model predictions are fixed.

## Directional allocation model

1. At the probe and on the actual saved time grid, synthesize each initial
   spatial mode with its known linear frequency. Group modes into fixed
   angular cells, forming a complex prior signal for each angular cell.
2. Project these prior signals onto the same temporal Fourier basis as the
   observed input. Denote the coefficient of angular cell j at frequency m
   by P_mj. The observed analytic coefficient is B_m.
3. Set W_mj=P_mj/sum_j(P_mj), then C_mj=B_m*W_mj. Thus sum_j(C_mj)=B_m exactly,
   while the initial prior's relative complex directional composition is
   retained. This is a first-order input constraint, not a fit to eta22/psi22.

The initial spectrum's linearly propagated *total* record is not required to
match the observed one. It is kept as a separate control. The model still
assumes its conditional directional structure is useful after evolution.
It does not claim to recover nonlinear directional redistribution.

Angular cells of 15 and 7.5 degrees were specified before target scoring.
Their centers discretize the directional distribution; this discretization
is reported, not disguised as exact recovery of the original angles.
Each observed frequency gets its finite-depth radial wave number. It is not
rounded onto the initial spatial FFT grid. The two trials use 444 and 888
direction-frequency parents, respectively, on 37 temporal frequency bins.
The common first-order domain is kh>=0.3 and below the input Nyquist. Initial
spectral energy retained is 99.9999996471%; observed temporal positive-frequency
energy retained is 99.9999999717%. Projection L2 is 0.00168186%.

Directional GL uses the released eta22 rank-8 kernel with full vector dot
products and vector sum magnitudes, and the released true-surface psi22
GL2+2 kernel including its Taylor term. No Stokes/angle correction, external
high-order field, MF12 reference or fitted gain enters either executor.

The initial-only control uses the same 15-degree angular discretization.
The direction-blind control uses the observed eta1 but sets all directions
to zero. Neither control is selected according to its error.

## Propagation convention audit

The historical directional generator contains `exp(i*(kx*x+ky*y+omega*t))`.
The initial positive-kx eta/psi spectra confirm positive time exponent:
relative polarization error is 7.0166e-5 for +i*g/omega and approximately 2
for -i*g/omega. This check uses initial first-order data only.

Convert the stored positive-k/+omega analytic representation to the physical
negative-k/-omega representation by conjugation. Rotate both horizontal axes
by pi to retain the released strict-forward convention; pair dot products
and sum magnitudes are unchanged. No high-order score, probe time shift or
fitted sign chooses this convention. The initial psi field is used only to
audit the stored convention, not as a predictor input.

An initial exploratory execution before this audit used the wrong prior time
sign and produced poor results. Its root-directory artifacts are retained
for debugging and are **not valid physical results**. Use only the
`verified_convention/` result directory below.

## Checks and conditioning

Wolfram passes five checks: allocation sum, identity when observed equals
prior, rotational invariance, and eta/psi temporal-conjugation identities.
MATLAB compares non-collinear complex modes against the existing spatial
GL interface at two times and three positions: maximum absolute errors
2.4243e-19 m for eta22 and 4.4124e-19 m^2/s for psi22. Exact cancellation in
the denominator is explicitly rejected; no arbitrary denominator floor is used.

Measured allocation sum relative errors are 2.02e-16 / 2.06e-16. These test
input consistency, not physical truth of the direction allocation.
The observed-energy-weighted sum(abs(W)) is 1.00000034 / 1.00000117, but
the maximum over all bins is 3575.88 / 8022.49. Extremely weak tail bins can
therefore be ill-conditioned. No clipping or target-selected cutoff was used;
robustness of this division in stronger/broader cases remains unresolved.
The fraction of observed input energy in bins with sum(abs(W))>10 is
1.21094e-9 / 2.56135e-9 for the two angular grids. These are diagnostics only;
no such bins were removed. All new MATLAB functions pass Code Analyzer.

## First-case results

The main window is +/-2Tp around the observed input envelope peak, selected
before high-order scoring. It contains only the original saved sample times.

| Method | eta22 main-group relative L2 (%) | psi22 main-group relative L2 (%) |
|---|---:|---:|
| Joint input, 15-degree cells | 0.645789 | 0.403072 |
| Joint input, 7.5-degree cells | 0.266808 | 0.140161 |
| Initial spectrum only, 15-degree cells | 2.510962 | 3.455475 |
| Observed eta1, all directions zero | 25.703151 | 21.554663 |

Initial-only first-order probe discrepancy is 1.4923% on the full projected
record; it was not a pass/fail gate. Angular refinement changes the main-group
prediction by 0.3958% for eta22 and 0.3437% for psi22. More angular convergence
evidence is needed before claiming a converged directional discretization.

Full-window errors are also retained: the joint 7.5-degree trial is 19.744%
for eta22 and 22.698% for psi22. The reported sub-percent result is therefore
specifically a main-group result, not a whole-record certificate.

This weak-steepness case supports the usefulness of the joint-input model at
the sampled main-group times. It does not establish high-steepness validity,
uniqueness of direction recovery, pure perturbation-order separation, or a
resolved continuous high-harmonic time history. Nonlinear output sums above
the temporal Nyquist are evaluated analytically at the original sample times;
no interpolation of OW3D is used and plot lines only connect those samples.

## Reproduction and artifacts

Add `research/directional_wave_data` to the MATLAB path. Execute:

```matlab
extract_directional_joint_input('C:/Research/VWA/VWA Unidirectinal/Directional/test1');
audit_directional_initial_convention('C:/Research/VWA/VWA Unidirectinal/Directional/test1');
test_directional_joint_model;
run_directional_joint_pilot;
```

Run `wolframscript -file research/directional_wave_data/verify_directional_joint_model.wl`.
Current results: `results/directional_joint_input/verified_convention/`.
Extraction and per-file hashes: `results/directional_joint_input/extracted.mat`
and `source_files.csv`. Checks/logs: `artifacts/directional_joint_input/`.
