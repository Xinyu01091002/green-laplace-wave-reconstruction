# Accepted directional information: initial linear spatial spectrum

The user permits directional assumptions and proposes knowing the initial
linear spatial wave-number spectrum. This supplies the directional input for
the next study, without renaming the project or launching new OW3D runs.

## Deterministic baseline interpretation

For the first forward baseline, explicitly interpret this as the complex
linear spectrum at a declared initial time: wave vectors (kx,ky), complex
amplitudes/phases, water depth, gravity, spatial origin and time origin.
One-sided analytic and two-sided real-field FFT normalizations must be stated.
The total nonlinear OceanWave3D.init field is not automatically this linear
spectrum; its true first-order component must be sourced independently or
labelled as a phase-sector approximation.

The linear dispersion relation determines omega from |k|. Linear propagation
of each declared parent then produces a predicted first-order record at the
probe. The directional GL kernels retain pair/triple wave-vector sums and
dot products; temporal output frequencies are the corresponding signed sums
of the parent frequencies. No averaged propagation angle replaces the actual
parent directions. High-order fields are generated only from those parents.

First compare this predicted first-order time record to the OW3D first phase
sector at the same location and saved times, without fitted shifts or scaling.
This distinguishes disagreement in assumed linear parent evolution from a
higher-order GL reconstruction error. The forward baseline is not a nonlinear
evolution solver. The local 4 s snapshots support same-time sample comparisons,
not a claim of resolved high-harmonic temporal FFT data.

## Relation to the eta1(t)-driven objective

Knowing a complete initial complex spectrum permits a deterministic forward
baseline. If the intended input remains the *observed* eta1(t), and the initial
spectrum is used only to supply directional structure, an additional model
is needed to distribute each observed temporal component among its directions.
A fixed directional structure could be a declared assumption, but it has not
been derived, implemented or validated by accepting the initial spectrum.
It must not be silently replaced by a unique angle for each frequency.

If only the initial energy spectrum S(kx,ky) is known, realization phases are
still missing. This is a distinct, weaker assumption and does not uniquely
determine a deterministic bound-wave time record.

## Next bounded comparison

Use the existing small-amplitude directional kh=1, spread=25, Akp=0.02 family.
Identify and verify its actual first-order initial spectrum and conventions;
audit strict-forward support against the released executor's domain. Establish
linear probe-record consistency, then eta22 and true surface psi22 spatial/
temporal parity, and compare at OW3D's saved times. Preserve the larger-steepness
and crossing-wave cases for subsequent checks. No directional eta20 zero-
frequency projection is inherited from the unidirectional trial: zero temporal
frequency can have a nonzero spatial difference wave vector.
