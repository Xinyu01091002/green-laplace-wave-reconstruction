# Joint input: initial linear spatial spectrum and observed eta1(t)

The user permits directional assumptions and proposes knowing the initial
linear spatial wave-number spectrum together with the probe's first-order
time series. The user clarified that linear propagation of the initial
spectrum need not reproduce the probe record at larger steepness. Matching
that linear prediction is therefore not a prerequisite for this research.
No project rename or new OW3D run is implied.

Execution update: the first prescribed conditional-amplitude trial has now
been run. See [JOINT_INPUT_PILOT.md](JOINT_INPUT_PILOT.md) for the explicit
allocation formula, source-convention audit, conditioning and limits. The
assumption below remains a model; one successful weak case is not a general
direction-recovery certification.

## Roles of the two inputs

When available, retain the complex linear spectrum at a declared initial
time: wave vectors (kx,ky), complex
amplitudes/phases, water depth, gravity, spatial origin and time origin.
One-sided analytic and two-sided real-field FFT normalizations must be stated.
The total nonlinear OceanWave3D.init field is not automatically this linear
spectrum; its true first-order component must be sourced independently or
labelled as a phase-sector approximation.

The initial spatial spectrum supplies a prior on directional structure and,
if complex coefficients are available, relative directional phases. The
observed eta1(t) supplies the actual local temporal amplitude and phase after
propagation. It is not replaced by the initial spectrum's linearly propagated
probe signal. The two inputs have distinct roles rather than competing as
alternative predictions of the same first-order record.

Linear propagation of the initial spectrum may be retained as an optional
diagnostic or synthetic implementation fixture. Its disagreement with the
observed probe does not reject the joint-input reconstruction and must not
trigger fitted time shifts or use of high-order references to change the
input. This study is not a nonlinear evolution solver.

## Relation to the eta1(t)-driven objective

The primary goal is now the observed eta1(t) plus initial directional structure.
An explicit model is still needed to distribute each observed temporal
component among its directions. A first candidate is to retain the initial
relative directional structure over the main-group window, while enforcing
that the sum of directional complex amplitudes at the probe reproduces the
observed temporal coefficient. This assumption has not yet been derived,
implemented or validated. It is not an assertion that nonlinear directional
redistribution is absent, nor a unique inversion from the available inputs.

Directional energy weights are not automatically complex amplitude weights.
Spatial phase at the probe and the common time origin must be included.
Frequency-to-wave-number association also needs a declared convention; no
average angle or output-frequency linear-dispersion substitution may replace
true vector interactions. If initial directional components nearly cancel
at the probe, dividing by their summed complex amplitude can be ill-conditioned.
That condition must be audited rather than hidden with an arbitrary floor.
No model coefficient may be selected against an OW3D high-order reference.

If only the initial energy spectrum S(kx,ky) is known, realization phases are
still missing. This is a distinct, weaker assumption and does not uniquely
determine a deterministic bound-wave time record.

## Next bounded comparison

Use the existing small-amplitude directional kh=1, spread=25, Akp=0.02 family.
Identify and verify its initial first-order spectrum and the observed local
first-order record; audit strict-forward support against the released domain.
Define and validate the directional allocation model and its exact recovery
of the supplied eta1(t), then test eta22 and true surface psi22. A known-parent
synthetic case can check vector-kernel implementation independently. The local
4 s snapshots provide spatial fields but are not automatically an adequately
resolved temporal input; probe-record sampling must be checked separately.
Preserve larger-steepness and crossing-wave cases for subsequent checks.
No directional eta20 zero-
frequency projection is inherited from the unidirectional trial: zero temporal
frequency can have a nonzero spatial difference wave vector.
