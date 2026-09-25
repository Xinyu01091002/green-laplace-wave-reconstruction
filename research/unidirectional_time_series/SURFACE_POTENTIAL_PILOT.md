# First surface-potential time-series pilot

Output is the true surface potential `psi(t)=phi(x,eta(x,t),t)`, in m^2/s,
not the flat bulk-potential trace or velocity. Inputs are precisely the saved
eta22 parent amplitudes, frequencies and wave numbers from the verified
boundary-corrected kh=1, Alpha=1 cases, Akp=0.02/0.12. The OW3D psi records
are held-out outputs and never supplied to the GL executor.

Linear prediction: `Re sum (-i*g/omega_j)*A_j*exp(-i*omega_j*(t-t0))`.
For psi22, `gl_psi22_time_pairs` transfers the released
`finite_depth_directional_dual_branch_green_laplace_psi22` graph. It uses
two prescribed Gauss-Laguerre nodes on each slow/fast branch, with scales
`2*nu_p - nu(2qp)` and `2*nu_p + nu(2qp)`. It retains the surface evaluation
term `0.5*eta1*phi1_z`; in the symmetric ordered-pair kernel this is
`-i*(nu_i+nu_j)/4`. Physical conversion multiplies the dimensionless kernel
and A_i*A_j by sqrt(g/h). No angular or Stokes repair is applied.

The spectral MF12 comparison requests **order 2 only**, using the same
surface coefficients `i*mu_2*(A_2+i*B_2)` and `i*mu_npm*(A_npm+i*B_npm)`
as `src/gl_mf12_order2_comparison.m`. No third-order MF12 call is made.

OW3D first and second psi phase sectors use the same temporal Hilbert/four-phase
operators as eta. No reference filtering, fitted gauge constant, sign, scale
or time shift is used. The common eta1-envelope maximum fixes the +/-2Tp
main-group window. Higher-order first/second-harmonic content remains present
in OW3D and is not called exact perturbation-order truth.

## Checks

- Wolfram: dual-branch exponential identity and symmetric surface Taylor term,
  both true (`verify_psi22_transfer.wl`).
- MATLAB: three times, three spatial points, nonzero complex phases and h=0.7;
  maximum GL temporal/spatial absolute difference 6.5278e-19 m^2/s.
- MF12 temporal/spatial maximum absolute difference 4.3368e-19 m^2/s.
- New MATLAB files pass Code Analyzer; input/output fields are finite.

## Main-group relative L2 against OW3D

| Prediction | Akp=0.02 (%) | Akp=0.12 (%) |
|---|---:|---:|
| Linear psi11 from eta1 | 0.026880 | 0.900116 |
| GL2+2 psi22 | 0.259661 | 2.085990 |
| Spectral MF12 psi22 | 0.054544 | 1.976162 |

The psi11 check is a linear polarization test of the declared eta1 input;
it is not a nonlinear GL correction. GL2+2 denotes two nodes in each branch,
not fourth perturbation order. The evidence supports surface-potential
reconstruction on these cases without claiming psi20/psi33 or arbitrary-depth
bulk-potential validation.

```powershell
wolframscript -file research/unidirectional_time_series/verify_psi22_transfer.wl
matlab -batch "restoredefaultpath; addpath('research/unidirectional_time_series'); test_surface_potential_time('C:/Research/spectral domain implementation of wave interaction theory'); run_surface_potential_pilot('C:/Research/spectral domain implementation of wave interaction theory');"
```

Results: `results/unidirectional_time_series/surface_potential_boundary_alpha1_akp002/`
and `..._akp012/`, with fields, raw metrics, source metadata, PNG/PDF main-group
figures and execution hashes. All earlier elevation results remain intact.
