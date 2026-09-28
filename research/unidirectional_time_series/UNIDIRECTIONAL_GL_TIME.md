# Single-direction GL time-series operators, 2026-09-27

Active entry point: `gl_unidirectional_time_series.m`.
It takes a declared first-order elevation record and returns the positive
sum components eta22, and optionally eta33. It preserves the original
prescribed Gauss--Laguerre rule, filtered-field products and FFT operators.
It does not enumerate parent pairs/triples, call an interaction reference
executor, fit a numerical kernel, or add Stokes corrections.

The user's new prohibition applies to verification as well as production.
It is recorded at the top of AGENTS.md and the current handoff. Historical
enumerating code and its outputs remain provenance only. None was executed
or used as a numerical comparator in this implementation work.

## What changed mathematically

For positive analytic time harmonics, the GL damping is an analytic time
translation. Consequently it preserves products:

\[
 \mathcal T_\tau(fg)=(\mathcal T_\tau f)(\mathcal T_\tau g),
 \qquad \widehat{\mathcal T_\tau f}(\Omega)
 =e^{-\tau\Omega\sqrt{h/g}}\hat f(\Omega).
\]

The original parent damping can therefore be moved to the output temporal
frequency after the source products are formed. For s=Omega sqrt(h/g) and
a=sqrt(K tanh K), each original second-order GL node is evaluated as

\[
e^{-s\tau}\{-a\sinh(a\tau)S_d+\cosh(a\tau)S_k\}
=\tfrac12\{(-aS_d+S_k)e^{-(s-a)\tau}
+(aS_d+S_k)e^{-(s+a)\tau}\}.
\]

The same node weights and peak-dependent scale are retained. This is the
prescribed finite-node GL graph, **not** the analytic-integral replacement
used in the earlier resolvent experiments. The node loop acts on output
wavenumber/frequency arrays, not on combinations of input waves.

For eta33, all required lower eta2 and flat Phi2 fields are generated from
the input. The original nested graph's distinct scales
2nu_p-a(2q_p) and 3nu_p-a(3q_p) are retained. Its internal eta2 is not
replaced by the separately returned determinant-root-scale eta22 field.
The original cubic forcing and its factors 1/2 and 1/4 are preserved.
Analytic time translation commutes through these generated lower fields,
so no outer-node repetition of the full lower-field calculation is needed.

The six symbolic checks in `verify_gl_time_node_transfer.wl` pass, including
the original node identities, lower-potential and outer-cubic transfers,
balance cancellation and the native test's dispersion relation.

## Computational structure

1. Fourier transform the supplied eta1 record and solve the linear
   finite-depth dispersion relation once per retained frequency.
2. Synthesize exact parent phases on a one-dimensional auxiliary spatial
   grid. No parent wavevector is rounded to that grid.
3. Form a fixed collection of filtered first-order fields using temporal
   FFTs, and construct the original source terms by pointwise products.
4. Transform source fields in time and space and apply the original GL node
   multipliers at each output (K,Omega).
5. For order three, form internally generated lower fields, construct the
   cubic sources, and apply the original outer GL node multipliers.
6. Evaluate the probe output and synthesize at the original time samples.

Support bounds use convex dispersion and endpoint formulas with work linear
in output-bin count. They do not evaluate or accumulate interaction kernels.
Arrays are auxiliary-space by temporal-frequency/time, plus space by parent
frequency for synthesis. There is no parent-count-squared or cubed kernel
array. For fixed auxiliary grid Nx and GL rank J, work has the form
O(Nx M + Nx W log W + Nx B log Nx + J Nx B), with fixed additional fields
for eta33. Nx must still be selected by convergence; this is not a uniform
cost theorem across arbitrary bandwidths and tolerances.

The finite auxiliary domain requires an output-response extension outside
the physically possible sum support. Its error is checked by successive
domain enlargement. `audit.converged` describes that change test, not a
rigorous absolute error bound. Original spatial operators remain unchanged.

## Usage

From the retained project root, with eta1 in metres, t in seconds, h in
metres, g in m/s^2, and kp in rad/m:

```matlab
addpath('research/unidirectional_time_series');
omega_p = sqrt(g*kp*tanh(kp*h));
options = struct( ...
    'order', 3, ...                       % 2 gives eta22 only
    'quadrature_rank', 8, ...
    'omega_max', 4*omega_p, ...            % explicitly declared input band
    'peak_wavenumber', kp, ...
    'domain_lengths', [25 50 100 200 400 800 1600], ... % lengths / h
    'relative_tolerance', 5e-4, ...
    'memory_budget_MiB', 4096);
[result, audit] = gl_unidirectional_time_series(eta1(:),t(:),g,h,options);
assert(audit.converged, 'Auxiliary-domain convergence not established.');
eta22 = result.eta22;
eta33 = result.eta33;
```

The core entry point writes no result files and can be called repeatedly.
Default order is 2 and default quadrature rank is 8. Shared inner/outer
cubic ranks 4/6/8 are supported; order-two ranks 4/6/8/12/16 are accepted.
`omega_max` must be declared. No kh cutoff or amplitude-ranking cutoff is
applied; counts below the previously validated model domains are reported
as extrapolation diagnostics. DC and the real-record Nyquist component are
reported separately. This is not an inversion from total measured eta.

Order two projects outputs onto the native representable positive temporal
band while retaining declared parents. The current order-three implementation
requires all declared cubic sums within that band and errors rather than
silently reducing the input support. Both return the original time vector.
An exact smaller internal complex FFT is used when the sum-degree bound
permits it. eta20, eta31, surface potential, directionality and nonlinear
time evolution are outside this entry point's scope.

## Validation without interaction enumeration

The native fixture has two nonzero complex first-order modes. Choose
tanh(q0)=sqrt(7)/3 and k2=2k1; then omega2/omega1=3/2, so both spatial and
temporal lattices are exactly representable. This permits a direct
comparison against the ORIGINAL spatial FFT-GL graphs at every sampled
time, without changing the input representation or evaluating pair kernels.

Measured relative differences:

- eta22 versus `green_laplace_eta22`, GL8: 2.16e-15.
- eta33 versus `gl_no_stokes_eta33`, GL8/GL8: 2.52e-12.
- First-order input reconstruction: 8.33e-16.
- Native-band check: an input whose nonlinear outputs are out of band has
  no effect on the retained output, to 4.02e-15 relative.
- Quadratic/cubic amplitude homogeneity and time-origin consistency pass.

MATLAB runtime profiling records all called functions and rejects the known
enumerating GL reference entry points. Original native cubic audits also
report zero pair and triple loops. Those checks passed. The original FFT
graphs are the numerical GL comparators; no retired time-pair/triple code
or saved GL prediction array is used.

## Independent raw-record demonstration

The demonstration explicitly loads only t, raw four-phase records and
physical metadata from the saved boundary-corrected single-direction OW3D
cases (kph=1, Alpha=1). First/second/third phase sectors are reconstructed
from those raw records. The input is every nonzero bin through 4omega_p:
641 samples and 64 temporal bins. No high-order reference enters the API.

Final main-group relative L2 against the raw observed phase sectors:

| Akp | eta22 | eta33 | Last domain-change diagnostic |
| ---: | ---: | ---: | ---: |
| .02 | .11329% | 2.08439% | .04940% at L/h=800 |
| .12 | 3.71986% | 10.18406% | .01971% at L/h=1600 |

The main group is fixed from the input envelope, +/-2Tp, with no alignment
or amplitude fitting. These are measured examples, not exact perturbation-
order certification. The first record is a first-harmonic phase sector,
and the references can contain higher-order contributions. Five retained
parents have kh<.3 and eight have kh<=.5; their inclusion is explicit and
does not extend the published physical validation domain generally.

The high case initially missed the .05% grid-change criterion at L/h=800
(.0548%). It was retained and then checked at L/h=1600, which passed.
The preliminary looser-tolerance second-order run also remains saved;
its closer OW3D score was not chosen as the preferred result. The sizeable
high-steepness third-sector discrepancy remains a limitation.

## Measured cost

Fresh local MATLAB R2022b process, low-case input, no profiling during timing.
Each call includes input FFT, dispersion, node preparation and all automatic
domain-refinement levels; MATLAB startup and raw-record extraction are
excluded. No operator plan is cached by the API.

| Call | Measured time |
| --- | ---: |
| First eta22 call | 1.201 s |
| eta22, median of 3 repeats | .782 s |
| Explicit changed-support eta22 call, 48 instead of 64 bins | .257 s |
| First combined eta22+eta33 call in the same process | 5.249 s |
| Combined eta22+eta33, median of 3 repeats | 5.057 s |

The changed-support call is timing-only; it does not alter the demonstrated
4omega_p input. The cubic first call follows the second-order calls in the
same process and is not described as an independent machine-cold launch.
No speed comparison to a prohibited interaction enumerator was made.

## Evidence and reproduction

Under `artifacts/unidirectional_time_series/`:

- `gl_time_node_transfer.json`: six symbolic identities.
- `unidirectional-gl-operator-v3/`: native eta22 comparison and call trace.
- `unidirectional-gl-cubic-operator/`: native eta33 comparison and call trace.
- `unidirectional-gl-demo-refined/`: eta22 raw-data run and unprofiled timing.
- `unidirectional-gl-cubic-demo/`: cubic records, high-case refinement,
  call traces and `final_demo.png/.pdf`.

The native tests are self-contained in this repository and accept a fresh
output directory to preserve previous evidence:

```matlab
test_gl_unidirectional_time_series('artifacts/my-new-native-eta22-check');
test_gl_unidirectional_cubic_operator('artifacts/my-new-native-eta33-check');
```

Observed-record demos require the existing local raw data and refuse to
overwrite their evidence folders. No remote simulation or postprocessing
job was launched in this implementation step. Historical pair-based
experiment commands elsewhere in the repository are superseded by AGENTS.md
and must not be executed under the current user instruction.
