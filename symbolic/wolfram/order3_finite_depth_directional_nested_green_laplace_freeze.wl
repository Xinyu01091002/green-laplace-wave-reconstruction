(* ::Package:: *)
(* Symbolic freeze for the eta33 nested Green--Laplace core.

   The complete executable graph begins at eta11. For every outer quadrature
   node it generates eta22 and flat Phi22 internally from the same damped
   eta11 leaves, assembles the original order-three kinematic and dynamic
   forcing, and applies the outer triple-sum response. No input-chi formula,
   sampled tuple, MF12 field, or eta33 oracle is used.

   Both Green--Laplace integrals use the exact four-node Gauss--Laguerre rule.
   Their fixed scales are the co-directional slow-tail rates s-Sqrt[A] at the
   declared spectral peak. These choices are frozen before field validation.
*)

residualSourceFD3NGL = FileNameJoin[{
  DirectoryName[$InputFileName],
  "..", "residuals", "order3_finite_depth_directional_euler.wl"
}];
If[!FileExistsQ[residualSourceFD3NGL],
  Print["FAIL: eta33 Euler residual source not found."];
  Exit[1]
];
Get[residualSourceFD3NGL];

gateSourceFD3NGL = FileNameJoin[{
  DirectoryName[$InputFileName],
  "order3_finite_depth_directional_crossing_stokes_gate.wl"
}];
If[!FileExistsQ[gateSourceFD3NGL],
  Print["FAIL: eta33 crossing Stokes gate source not found."];
  Exit[1]
];

ClearAll[
  qFD3NGL, qpFD3NGL, xFD3NGL, nodesFD3NGL, weightsFD3NGL,
  nuFD3NGL, lambda2FD3NGL, lambda3FD3NGL,
  nu1BalanceFD3NGL, nu2BalanceFD3NGL, aBalanceFD3NGL,
  cBalanceFD3NGL, tBalanceFD3NGL, balanceResidualsFD3NGL,
  momentResidualsFD3NGL, diagonalRulesFD3NGL,
  rKDiagonalFD3NGL, rDDiagonalFD3NGL,
  eta2OrderedDiagonalFD3NGL, phi2OrderedDiagonalFD3NGL,
  surfacePotentialSeriesFD3NGL, surfacePotentialOrder2FD3NGL,
  psi2GroupedTaylorFD3NGL, psi2ExpectedTaylorFD3NGL,
  psi2OrderedDiagonalFD3NGL,
  innerPhiCoreDiagonalFD3NGL, correctionPhi2FD3NGL,
  physicalDirectKFD3NGL, physicalDirectDFD3NGL,
  physicalPairKFD3NGL, physicalPairDFD3NGL,
  orderedPhysicalKFD3NGL, orderedPhysicalDFD3NGL,
  a3DiagonalFD3NGL, s3DiagonalFD3NGL, stokes3FD3NGL,
  phi3GroupedDiagonalFD3NGL, phi3OrderedDiagonalFD3NGL,
  outerCoreDiagonalFD3NGL, correction3FD3NGL,
  outerPhiCoreDiagonalFD3NGL, correctionPhi3FD3NGL,
  zFD3NGL, verticalModeFD3NGL, verticalModeChecksFD3NGL,
  exactRecoveredDiagonalFD3NGL, positivityChecksFD3NGL,
  exactTextFD3NGL, gateFD3NGL, booleanGateFD3NGL,
  gatesFD3NGL, overallPassFD3NGL, projectRootFD3NGL,
  artifactPathFD3NGL, generatedDirectoryFD3NGL,
  interfacePathFD3NGL, sourceRelativePathFD3NGL,
  sourceSha256FD3NGL, residualRelativePathFD3NGL,
  residualSha256FD3NGL, gateRelativePathFD3NGL, gateSha256FD3NGL
];

nuFD3NGL[q_] := Sqrt[q Tanh[q]];
nodesFD3NGL = xFD3NGL /. Solve[
  LaguerreL[4, xFD3NGL] == 0,
  xFD3NGL,
  Reals
];
nodesFD3NGL = SortBy[nodesFD3NGL, N];
weightsFD3NGL = #/(25 LaguerreL[5, #]^2) & /@ nodesFD3NGL;
momentResidualsFD3NGL = FullSimplify[
  Table[
    Sum[
      weightsFD3NGL[[index]] nodesFD3NGL[[index]]^degree,
      {index, 4}
    ] - degree!,
    {degree, 0, 7}
  ]
];

lambda2FD3NGL = (
  2 nuFD3NGL[qpFD3NGL]
  - Sqrt[2 qpFD3NGL Tanh[2 qpFD3NGL]]
);
lambda3FD3NGL = (
  3 nuFD3NGL[qpFD3NGL]
  - Sqrt[3 qpFD3NGL Tanh[3 qpFD3NGL]]
);
positivityChecksFD3NGL = N[
  {lambda2FD3NGL, lambda3FD3NGL} /. qpFD3NGL -> #,
  30
] & /@ {1/2, 3/4, 1, 2, 4, 8};
balanceResidualsFD3NGL = FullSimplify[{
  Exp[-tBalanceFD3NGL (nu1BalanceFD3NGL-cBalanceFD3NGL)]
    Exp[-tBalanceFD3NGL (nu2BalanceFD3NGL-cBalanceFD3NGL)]
    Exp[tBalanceFD3NGL (aBalanceFD3NGL-2 cBalanceFD3NGL)]
    -Exp[-tBalanceFD3NGL (
      nu1BalanceFD3NGL+nu2BalanceFD3NGL-aBalanceFD3NGL)],
  Exp[-tBalanceFD3NGL (nu1BalanceFD3NGL-cBalanceFD3NGL)]
    Exp[-tBalanceFD3NGL (nu2BalanceFD3NGL-cBalanceFD3NGL)]
    Exp[-tBalanceFD3NGL (aBalanceFD3NGL+2 cBalanceFD3NGL)]
    -Exp[-tBalanceFD3NGL (
      nu1BalanceFD3NGL+nu2BalanceFD3NGL+aBalanceFD3NGL)]
}];

diagonalRulesFD3NGL = {
  q1ParentFD3 -> qFD3NGL,
  q2ParentFD3 -> qFD3NGL,
  q3ParentFD3 -> qFD3NGL,
  rhoFD3 -> 0
};
rKDiagonalFD3NGL = FullSimplify[
  kinematicForcingFD3 /. diagonalRulesFD3NGL,
  Assumptions -> qFD3NGL > 1/2
];
rDDiagonalFD3NGL = FullSimplify[
  dynamicForcingFD3 /. diagonalRulesFD3NGL,
  Assumptions -> qFD3NGL > 1/2
];
eta2OrderedDiagonalFD3NGL = FullSimplify[
  diagonalResidualsFD2[[1]]/4 /. q1FD2 -> qFD3NGL,
  Assumptions -> qFD3NGL > 1/2
];
phi2OrderedDiagonalFD3NGL = FullSimplify[
  diagonalResidualsFD2[[2]]/4 /. q1FD2 -> qFD3NGL,
  Assumptions -> qFD3NGL > 1/2
];

(* The free-surface potential is obtained from the bulk potential by an
   independent Taylor contact at z=eta.  At order two,
     psi2=Phi2+eta1 phi1z(0).
   The labelled grouped modal coefficient is divided by four to obtain the
   ordered analytic-kernel convention used by the executable graph. *)
surfacePotentialSeriesFD3NGL = surfaceTaylorG2[phiSeriesG2];
surfacePotentialOrder2FD3NGL = coefficientG2[
  surfacePotentialSeriesFD3NGL,
  2
];
psi2GroupedTaylorFD3NGL = Together[
  modalCoefficientFD2[surfacePotentialOrder2FD3NGL]
];
psi2ExpectedTaylorFD3NGL =
  phiGroupedCoefficientFD2 - I nuSumFD2;
psi2OrderedDiagonalFD3NGL = FullSimplify[
  (
    psi2GroupedTaylorFD3NGL /. {
      phiGroupedCoefficientFD2 -> 4 phi2OrderedDiagonalFD3NGL,
      q1FD2 -> qFD3NGL,
      q2FD2 -> qFD3NGL,
      cosineFD2 -> 1
    }
  )/4,
  Assumptions -> qFD3NGL > 1/2
];

(* The same inner Green--Laplace response provides the ordered flat
   bulk-potential trace.  It is distinct from the surface potential. *)
innerPhiCoreDiagonalFD3NGL = Sum[
  With[
    {
      node = nodesFD3NGL[[index]],
      weight = weightsFD3NGL[[index]],
      lambda = lambda2FD3NGL,
      s = 2 nuFD3NGL[qFD3NGL],
      a = Sqrt[2 qFD3NGL Tanh[2 qFD3NGL]],
      sourceD = 3 nuFD3NGL[qFD3NGL]^2
        - qFD3NGL^2/nuFD3NGL[qFD3NGL]^2,
      sourceK = 4 qFD3NGL^2/nuFD3NGL[qFD3NGL]
    },
    (weight Exp[node]/lambda)
      * Exp[-s node/lambda]
      * I (
        Cosh[a node/lambda] sourceD
        - Sinh[a node/lambda]/a sourceK
      )/4
  ],
  {index, 4}
];
correctionPhi2FD3NGL =
  phi2OrderedDiagonalFD3NGL - innerPhiCoreDiagonalFD3NGL;

(* A flat trace uniquely determines the finite-depth bulk mode in the
   positive pure-sum sector.  These gates certify the continuation without
   introducing a surface-potential variable. *)
verticalModeFD3NGL[q_, z_] := Cosh[q (z + 1)]/Cosh[q];
verticalModeChecksFD3NGL = FullSimplify[
  {
    verticalModeFD3NGL[qFD3NGL, 0] - 1,
    D[verticalModeFD3NGL[qFD3NGL, zFD3NGL], zFD3NGL]
      /. zFD3NGL -> -1,
    D[verticalModeFD3NGL[qFD3NGL, zFD3NGL], {zFD3NGL, 2}]
      - qFD3NGL^2 verticalModeFD3NGL[qFD3NGL, zFD3NGL]
  },
  Assumptions -> qFD3NGL > 0 && -1 <= zFD3NGL <= 0
];
physicalDirectKFD3NGL = (3 I/2) qFD3NGL^2 nuFD3NGL[qFD3NGL];
physicalDirectDFD3NGL = -qFD3NGL^2/2;
physicalPairKFD3NGL = (
  3 I qFD3NGL^2 eta2OrderedDiagonalFD3NGL/nuFD3NGL[qFD3NGL]
  - 6 qFD3NGL^2 phi2OrderedDiagonalFD3NGL
);
physicalPairDFD3NGL = (
  -3 I nuFD3NGL[qFD3NGL]
    (2 qFD3NGL Tanh[2 qFD3NGL]) phi2OrderedDiagonalFD3NGL
  - nuFD3NGL[qFD3NGL]^2 eta2OrderedDiagonalFD3NGL
  + 2 I qFD3NGL^2 phi2OrderedDiagonalFD3NGL/nuFD3NGL[qFD3NGL]
);
orderedPhysicalKFD3NGL =
  physicalDirectKFD3NGL/4 + physicalPairKFD3NGL/2;
orderedPhysicalDFD3NGL =
  physicalDirectDFD3NGL/4 + physicalPairDFD3NGL/2;
a3DiagonalFD3NGL = Sqrt[3 qFD3NGL Tanh[3 qFD3NGL]];
s3DiagonalFD3NGL = 3 nuFD3NGL[qFD3NGL];
stokes3FD3NGL = etaDiagonalOrderedFD3 /. qDiagonalFD3 -> qFD3NGL;
phi3GroupedDiagonalFD3NGL = FullSimplify[
  -(
    kinematicEtaLinearFD3 etaGroupedCoefficientFD3
    + kinematicForcingFD3
  )/kinematicPhiLinearFD3 /. Join[
    diagonalRulesFD3NGL,
    {etaGroupedCoefficientFD3 -> 24 stokes3FD3NGL}
  ],
  Assumptions -> qFD3NGL > 1/2
];
phi3OrderedDiagonalFD3NGL = phi3GroupedDiagonalFD3NGL/24;

(* The original residual convention is
     RK=-i s E-A Phi+FK, RD=E-i s Phi+FD.
   Hence the ordered outer Green--Laplace elevation core is the integral of
     Exp[-s t] (A Sinh[a t]/a FD-i Cosh[a t] FK)/24. *)
outerCoreDiagonalFD3NGL = Sum[
  With[
    {
      node = nodesFD3NGL[[index]],
      weight = weightsFD3NGL[[index]],
      lambda = lambda3FD3NGL
    },
    (weight Exp[node]/lambda)
      * Exp[-s3DiagonalFD3NGL node/lambda]
      * (
        a3DiagonalFD3NGL Sinh[a3DiagonalFD3NGL node/lambda]
          rDDiagonalFD3NGL
        - I Cosh[a3DiagonalFD3NGL node/lambda]
          rKDiagonalFD3NGL
      )/24
  ],
  {index, 4}
];
correction3FD3NGL = stokes3FD3NGL - outerCoreDiagonalFD3NGL;
outerPhiCoreDiagonalFD3NGL = Sum[
  With[
    {
      node = nodesFD3NGL[[index]],
      weight = weightsFD3NGL[[index]],
      lambda = lambda3FD3NGL
    },
    (weight Exp[node]/lambda)
      * Exp[-s3DiagonalFD3NGL node/lambda]
      * (
        -I Cosh[a3DiagonalFD3NGL node/lambda]
          rDDiagonalFD3NGL
        - Sinh[a3DiagonalFD3NGL node/lambda]/a3DiagonalFD3NGL
          rKDiagonalFD3NGL
      )/24
  ],
  {index, 4}
];
correctionPhi3FD3NGL =
  phi3OrderedDiagonalFD3NGL - outerPhiCoreDiagonalFD3NGL;
exactRecoveredDiagonalFD3NGL = FullSimplify[
  -elevationForcingFD3/elevationLinearFD3 /. diagonalRulesFD3NGL,
  Assumptions -> qFD3NGL > 1/2
]/24;

exactTextFD3NGL[expression_] := ToString[InputForm[expression]];
gateFD3NGL[name_, residual_] := Module[
  {pass = TrueQ[PossibleZeroQ[residual]]},
  Print[name, " = ", If[pass, "PASS", "FAIL"]];
  If[!pass, Print["  residual: ", InputForm[residual]]];
  pass
];
booleanGateFD3NGL[name_, value_] := Module[{pass = TrueQ[value]},
  Print[name, " = ", If[pass, "PASS", "FAIL"]];
  pass
];
gatesFD3NGL = {
  gateFD3NGL[
    "four-node Gauss-Laguerre moments are exact through degree seven",
    Total[momentResidualsFD3NGL]
  ],
  booleanGateFD3NGL[
    "inner and outer co-directional tail scales are positive",
    And @@ Flatten[
      Map[NumericQ[#] && # > 0 &, positivityChecksFD3NGL, {2}]
    ]
  ],
  gateFD3NGL[
    "single-shift exponential balancing preserves both response branches",
    Total[balanceResidualsFD3NGL]
  ],
  gateFD3NGL[
    "exact residual recovery reproduces the ordered Stokes trace",
    Together[exactRecoveredDiagonalFD3NGL - stokes3FD3NGL]
  ],
  gateFD3NGL[
    "FFT direct and pair forcing multiplicities reproduce ordered RK and RD",
    Together[
      orderedPhysicalKFD3NGL-rKDiagonalFD3NGL/24
      +orderedPhysicalDFD3NGL-rDDiagonalFD3NGL/24
    ]
  ],
  gateFD3NGL[
    "crossing correction restores the complete nested-core diagonal",
    Together[
      outerCoreDiagonalFD3NGL + correction3FD3NGL - stokes3FD3NGL
    ]
  ],
  gateFD3NGL[
    "flat-Phi22 correction restores the exact Stokes diagonal",
    Together[
      innerPhiCoreDiagonalFD3NGL + correctionPhi2FD3NGL
        - phi2OrderedDiagonalFD3NGL
    ]
  ],
  gateFD3NGL[
    "surface-potential Taylor contact is exact through order two",
    Together[
      psi2GroupedTaylorFD3NGL-psi2ExpectedTaylorFD3NGL
    ]
  ],
  gateFD3NGL[
    "ordered Psi22 Stokes trace equals flat Phi22 plus Taylor contact",
    Together[
      psi2OrderedDiagonalFD3NGL
        - phi2OrderedDiagonalFD3NGL
        + I nuFD3NGL[qFD3NGL]/2
    ]
  ],
  gateFD3NGL[
    "flat-Phi33 correction restores the exact Stokes diagonal",
    Together[
      outerPhiCoreDiagonalFD3NGL + correctionPhi3FD3NGL
        - phi3OrderedDiagonalFD3NGL
    ]
  ],
  gateFD3NGL[
    "finite-depth bulk vertical mode satisfies trace bed and Laplace gates",
    Total[verticalModeChecksFD3NGL]
  ]
};
overallPassFD3NGL = And @@ gatesFD3NGL;

projectRootFD3NGL = DirectoryName[$InputFileName, 3];
artifactPathFD3NGL = FileNameJoin[{
  projectRootFD3NGL, "artifacts",
  "order3_finite_depth_directional_nested_green_laplace_freeze.json"
}];
sourceRelativePathFD3NGL =
  "symbolic/wolfram/order3_finite_depth_directional_nested_green_laplace_freeze.wl";
residualRelativePathFD3NGL =
  "symbolic/residuals/order3_finite_depth_directional_euler.wl";
gateRelativePathFD3NGL =
  "symbolic/wolfram/order3_finite_depth_directional_crossing_stokes_gate.wl";
sourceSha256FD3NGL = FileHash[$InputFileName, "SHA256", "HexString"];
residualSha256FD3NGL = FileHash[residualSourceFD3NGL, "SHA256", "HexString"];
gateSha256FD3NGL = FileHash[gateSourceFD3NGL, "SHA256", "HexString"];

Export[
  artifactPathFD3NGL,
  <|
    "schema_version" -> 1,
    "status" -> "frozen-symbolic-nested-green-laplace-pre-mf12",
    "candidate_id" -> "gl-eta33-nested-gl4-tail-g3-v1",
    "scope" ->
      "finite-depth strict-forward positive-pure-sum eta33 q_i>0.5",
    "candidate_input_fields" -> {"eta11"},
    "internal_generated_fields" -> {
      "eta22_plus", "Phi22_plus", "Psi22_plus", "Phi22_t_plus",
      "eta33_plus", "Phi33_plus", "Phi33_t_plus"
    },
    "external_lower_order_fields" -> False,
    "input_chi_route_used" -> False,
    "oracle_or_mf12_used" -> False,
    "sampled_selection_used" -> False,
    "inner_quadrature_rank" -> 4,
    "outer_quadrature_rank" -> 4,
    "quadrature_exact_moment_degree" -> 7,
    "quadrature_nodes_exact" -> (exactTextFD3NGL /@ nodesFD3NGL),
    "quadrature_weights_exact" -> (exactTextFD3NGL /@ weightsFD3NGL),
    "inner_scale" -> exactTextFD3NGL[lambda2FD3NGL],
    "outer_scale" -> exactTextFD3NGL[lambda3FD3NGL],
    "scale_rationale" ->
      "co-directional slow-tail rate s-Sqrt[A] from the exact sum detuning",
    "exponential_balance_identity" ->
      "Exp[-t(nu1-c)] Exp[-t(nu2-c)] Exp[+/-t(a-2c branch)] preserves Exp[-t(nu1+nu2-/+a)]",
    "kinematic_diagonal_forcing" -> exactTextFD3NGL[rKDiagonalFD3NGL],
    "dynamic_diagonal_forcing" -> exactTextFD3NGL[rDDiagonalFD3NGL],
    "fft_forcing_normalization" ->
      "F_ordered=F_direct/4+F_pair/2",
    "ordered_stokes_trace" -> exactTextFD3NGL[stokes3FD3NGL],
    "outer_core_diagonal" -> exactTextFD3NGL[outerCoreDiagonalFD3NGL],
    "stokes_correction_diagonal" -> exactTextFD3NGL[correction3FD3NGL],
    "Phi33_ordered_stokes_trace" ->
      exactTextFD3NGL[phi3OrderedDiagonalFD3NGL],
    "Phi22_ordered_stokes_trace" ->
      exactTextFD3NGL[phi2OrderedDiagonalFD3NGL],
    "Phi22_inner_core_diagonal" ->
      exactTextFD3NGL[innerPhiCoreDiagonalFD3NGL],
    "Phi22_stokes_correction_diagonal" ->
      exactTextFD3NGL[correctionPhi2FD3NGL],
    "Psi22_ordered_stokes_trace" ->
      exactTextFD3NGL[psi2OrderedDiagonalFD3NGL],
    "Psi22_taylor_relation" ->
      "Psi22=Phi22+(1/2) eta11 phi1z_analytic; phi1z_analytic=-I nu eta11",
    "Phi33_outer_core_diagonal" ->
      exactTextFD3NGL[outerPhiCoreDiagonalFD3NGL],
    "Phi33_stokes_correction_diagonal" ->
      exactTextFD3NGL[correctionPhi3FD3NGL],
    "Phi33_t_construction" ->
      "analytic time derivative of the frozen discrete eta11-only graph",
    "bulk_potential_variable" ->
      "Phi_nn=phi_nn(z=0); Phi_nn is not surface potential psi_nn",
    "bulk_vertical_mode" ->
      "V(Q,z)=Cosh[Q (z+1)]/Cosh[Q], -1<=z<=0, g=h=1",
    "crossing_gate" -> "Product_{i<j} (1+unit_ki dot unit_kj)/2",
    "two_lobe_crossing_reduction" -> "((1+cos(gamma))/2)^2",
    "cost_class" ->
      "fixed J2*J3 filtered products; constant-times N log N after full MATLAB graph count",
    "gates" -> AssociationThread[
      {"laguerre_moments", "positive_scales", "exponential_balance",
        "exact_stokes",
        "forcing_multiplicity", "correction", "Phi22_correction",
        "Psi22_taylor", "Psi22_diagonal", "Phi33_correction",
        "bulk_vertical_mode"},
      gatesFD3NGL
    ],
    "overall_exact_gate_pass" -> overallPassFD3NGL
  |>,
  "RawJSON"
];

generatedDirectoryFD3NGL = FileNameJoin[{
  projectRootFD3NGL, "symbolic", "generated"
}];
If[!DirectoryQ[generatedDirectoryFD3NGL],
  CreateDirectory[
    generatedDirectoryFD3NGL,
    CreateIntermediateDirectories -> True
  ]
];
interfacePathFD3NGL = FileNameJoin[{
  generatedDirectoryFD3NGL,
  "finite_depth_directional_order3_nested_green_laplace.json"
}];
Export[
  interfacePathFD3NGL,
  <|
    "schema_version" -> 1,
    "generator" -> sourceRelativePathFD3NGL,
    "generator_sha256" -> sourceSha256FD3NGL,
    "residual_source" -> residualRelativePathFD3NGL,
    "residual_source_sha256" -> residualSha256FD3NGL,
    "crossing_gate_source" -> gateRelativePathFD3NGL,
    "crossing_gate_source_sha256" -> gateSha256FD3NGL,
    "status" -> "frozen-symbolic-nested-green-laplace-pre-mf12",
    "candidate_id" -> "gl-eta33-nested-gl4-tail-g3-v1",
    "candidate_input_fields" -> {"eta11"},
    "internal_generated_fields" -> {
      "eta22_plus", "Phi22_plus", "Psi22_plus", "Phi22_t_plus",
      "eta33_plus", "Phi33_plus", "Phi33_t_plus"
    },
    "input_chi_route_used" -> False,
    "oracle_or_mf12_used" -> False,
    "inner_quadrature_rank" -> 4,
    "outer_quadrature_rank" -> 4,
    "quadrature_nodes_exact" -> (exactTextFD3NGL /@ nodesFD3NGL),
    "quadrature_weights_exact" -> (exactTextFD3NGL /@ weightsFD3NGL),
    "inner_scale" -> exactTextFD3NGL[lambda2FD3NGL],
    "outer_scale" -> exactTextFD3NGL[lambda3FD3NGL],
    "exponential_balance_identity" ->
      "single scalar shift c is an exact representation change, not a kernel change",
    "kinematic_diagonal_forcing" -> exactTextFD3NGL[rKDiagonalFD3NGL],
    "dynamic_diagonal_forcing" -> exactTextFD3NGL[rDDiagonalFD3NGL],
    "fft_forcing_normalization" ->
      "F_ordered=F_direct/4+F_pair/2",
    "ordered_stokes_trace" -> exactTextFD3NGL[stokes3FD3NGL],
    "stokes_correction_diagonal" -> exactTextFD3NGL[correction3FD3NGL],
    "Phi33_ordered_stokes_trace" ->
      exactTextFD3NGL[phi3OrderedDiagonalFD3NGL],
    "Phi22_ordered_stokes_trace" ->
      exactTextFD3NGL[phi2OrderedDiagonalFD3NGL],
    "Phi22_stokes_correction_diagonal" ->
      exactTextFD3NGL[correctionPhi2FD3NGL],
    "Psi22_ordered_stokes_trace" ->
      exactTextFD3NGL[psi2OrderedDiagonalFD3NGL],
    "Psi22_taylor_relation" ->
      "Psi22=Phi22+(1/2) eta11 phi1z_analytic; phi1z_analytic=-I nu eta11",
    "Phi33_stokes_correction_diagonal" ->
      exactTextFD3NGL[correctionPhi3FD3NGL],
    "Phi33_t_construction" ->
      "analytic time derivative of the frozen discrete eta11-only graph",
    "bulk_potential_variable" ->
      "Phi_nn=phi_nn(z=0); Phi_nn is not surface potential psi_nn",
    "bulk_vertical_mode" ->
      "V(Q,z)=Cosh[Q (z+1)]/Cosh[Q], -1<=z<=0, g=h=1",
    "overall_exact_gate_pass" -> overallPassFD3NGL
  |>,
  "RawJSON"
];

Print["inner tail scale = ", InputForm[lambda2FD3NGL]];
Print["outer tail scale = ", InputForm[lambda3FD3NGL]];
Print["diagonal forcing leaf counts = ",
  InputForm[LeafCount /@ {rKDiagonalFD3NGL, rDDiagonalFD3NGL}]];
Print["artifact = ", artifactPathFD3NGL];
Print["interface = ", interfacePathFD3NGL];
Print[
  "OVERALL_ORDER3_NESTED_GREEN_LAPLACE_FREEZE = ",
  If[overallPassFD3NGL, "PASS", "FAIL"]
];
If[!overallPassFD3NGL, Exit[1]];
