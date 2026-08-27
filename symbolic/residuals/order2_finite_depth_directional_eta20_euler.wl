(* ::Package:: *)
(* Exact finite-depth directional order-two difference-frequency Euler
   residual. This is a residual and reference source, not a candidate. *)

unidirectionalSourceFD2D20 = FileNameJoin[{
  DirectoryName[$InputFileName],
  "order2_finite_depth_eta20_euler.wl"
}];
directionalExpansionSourceFD2D20 = FileNameJoin[{
  DirectoryName[$InputFileName],
  "order1_to_order3_symbolic_expansion_2d.wl"
}];
If[
  !FileExistsQ[unidirectionalSourceFD2D20]
  || !FileExistsQ[directionalExpansionSourceFD2D20],
  Print["FAIL: required finite-depth Euler source is missing."];
  Exit[1]
];

(* Preserve the independently derived one-dimensional ordered kernel before
   loading the generic two-dimensional expansion. *)
Get[unidirectionalSourceFD2D20];
unidirectionalOrderedFD2D20 = etaOrdereddFD20;
Get[directionalExpansionSourceFD2D20];

ClearAll[
  amplitude1FD2D20, amplitude2BarFD2D20,
  q1FD2D20, q2FD2D20, cosineFD2D20, sineFD2D20,
  nuFD2D20, dnoFD2D20, nu1FD2D20, nu2FD2D20,
  sigmaDeltaFD2D20, qxDeltaFD2D20, qyDeltaFD2D20,
  qDeltaFD2D20, dot12FD2D20,
  etaGroupedFD2D20, phiGroupedFD2D20,
  verticalJetFD2D20, eta1ModalFD2D20, phi1ModalFD2D20,
  eta2ModalFD2D20, phi2ModalFD2D20,
  zeroEtaFD2D20, zeroPhiFD2D20, modalRulesFD2D20,
  modalCoefficientFD2D20, assumptionsFD2D20,
  rKFD2D20, rDFD2D20,
  kinematicEtaLinearFD2D20, kinematicPhiLinearFD2D20,
  dynamicEtaLinearFD2D20, dynamicPhiLinearFD2D20,
  forcingKFD2D20, forcingDFD2D20, determinantFD2D20,
  etaNumeratorFD2D20, phiNumeratorFD2D20,
  etaGroupedExactFD2D20, phiGroupedExactFD2D20,
  etaOrderedExactFD2D20, phiOrderedExactFD2D20,
  rKRecoveredFD2D20, rDRecoveredFD2D20,
  etaResidualFD2D20,
  unidirectionalForcingKFD2D20, unidirectionalForcingDFD2D20,
  unidirectionalRegressionFD2D20,
  carrierQFD2D20, equalRadialKernelFD2D20,
  equalRadialAngleLimitFD2D20, radialOneSidedLimitFD2D20,
  strictMatrixFD2D20, strictRankFD2D20,
  exactTextFD2D20, gateFD2D20, gateResultsFD2D20,
  overallPassFD2D20, projectRootFD2D20,
  artifactDirectoryFD2D20, artifactPathFD2D20
];

assumptionsFD2D20 = (
  q1FD2D20 > q2FD2D20 > 0
  && -1 < cosineFD2D20 <= 1
);
nuFD2D20[q_] := Sqrt[q Tanh[q]];
dnoFD2D20[q_] := q Tanh[q];
sineFD2D20 = Sqrt[1-cosineFD2D20^2];
nu1FD2D20 = nuFD2D20[q1FD2D20];
nu2FD2D20 = nuFD2D20[q2FD2D20];
sigmaDeltaFD2D20 = nu1FD2D20-nu2FD2D20;
dot12FD2D20 = q1FD2D20 q2FD2D20 cosineFD2D20;
qxDeltaFD2D20 = q1FD2D20-q2FD2D20 cosineFD2D20;
qyDeltaFD2D20 = -q2FD2D20 sineFD2D20;
qDeltaFD2D20 = Sqrt[
  q1FD2D20^2+q2FD2D20^2-2 dot12FD2D20
];

verticalJetFD2D20[q_, zz_] := (
  1+dnoFD2D20[q] zz+q^2 zz^2/2
);
eta1ModalFD2D20[xx_, yy_, tt_] := (
  amplitude1FD2D20 Exp[
    I (q1FD2D20 xx-nu1FD2D20 tt)
  ]
  +amplitude2BarFD2D20 Exp[
    -I (
      q2FD2D20 cosineFD2D20 xx
      +q2FD2D20 sineFD2D20 yy
      -nu2FD2D20 tt
    )
  ]
);
phi1ModalFD2D20[xx_, yy_, zz_, tt_] := (
  -I amplitude1FD2D20/nu1FD2D20
    verticalJetFD2D20[q1FD2D20, zz]
    Exp[I (q1FD2D20 xx-nu1FD2D20 tt)]
  +I amplitude2BarFD2D20/nu2FD2D20
    verticalJetFD2D20[q2FD2D20, zz]
    Exp[
      -I (
        q2FD2D20 cosineFD2D20 xx
        +q2FD2D20 sineFD2D20 yy
        -nu2FD2D20 tt
      )
    ]
);
eta2ModalFD2D20[xx_, yy_, tt_] := (
  etaGroupedFD2D20 amplitude1FD2D20 amplitude2BarFD2D20
    Exp[
      I (
        qxDeltaFD2D20 xx+qyDeltaFD2D20 yy
        -sigmaDeltaFD2D20 tt
      )
    ]
);
phi2ModalFD2D20[xx_, yy_, zz_, tt_] := (
  phiGroupedFD2D20 amplitude1FD2D20 amplitude2BarFD2D20
    verticalJetFD2D20[qDeltaFD2D20, zz]
    Exp[
      I (
        qxDeltaFD2D20 xx+qyDeltaFD2D20 yy
        -sigmaDeltaFD2D20 tt
      )
    ]
);
zeroEtaFD2D20[xx_, yy_, tt_] := 0;
zeroPhiFD2D20[xx_, yy_, zz_, tt_] := 0;
modalRulesFD2D20 = {
  etaG2[1] -> eta1ModalFD2D20,
  etaG2[2] -> eta2ModalFD2D20,
  etaG2[3] -> zeroEtaFD2D20,
  phiG2[1] -> phi1ModalFD2D20,
  phiG2[2] -> phi2ModalFD2D20,
  phiG2[3] -> zeroPhiFD2D20
};
modalCoefficientFD2D20[expression_] :=
  Coefficient[
    Coefficient[
      Expand[
        expression /. modalRulesFD2D20 /. {
          gravityG2 -> 1,
          xG2 -> 0,
          yG2 -> 0,
          tG2 -> 0
        }
      ],
      amplitude1FD2D20,
      1
    ],
    amplitude2BarFD2D20,
    1
  ];

Print["Extracting finite-depth directional eta20 Euler residuals."];
rKFD2D20 = modalCoefficientFD2D20[kinematicActualG2[[2]]];
Print[
  "directional eta20 raw RK leaf count = ",
  LeafCount[rKFD2D20]
];
rDFD2D20 = modalCoefficientFD2D20[dynamicActualG2[[2]]];
Print[
  "directional eta20 raw RD leaf count = ",
  LeafCount[rDFD2D20]
];

kinematicEtaLinearFD2D20 = Coefficient[rKFD2D20, etaGroupedFD2D20];
kinematicPhiLinearFD2D20 = Coefficient[rKFD2D20, phiGroupedFD2D20];
dynamicEtaLinearFD2D20 = Coefficient[rDFD2D20, etaGroupedFD2D20];
dynamicPhiLinearFD2D20 = Coefficient[rDFD2D20, phiGroupedFD2D20];
forcingKFD2D20 = rKFD2D20 /. {
  etaGroupedFD2D20 -> 0,
  phiGroupedFD2D20 -> 0
};
forcingDFD2D20 = rDFD2D20 /. {
  etaGroupedFD2D20 -> 0,
  phiGroupedFD2D20 -> 0
};
Print["Assembling directional eta20 affine system."];
determinantFD2D20 = (
  kinematicEtaLinearFD2D20 dynamicPhiLinearFD2D20
  -dynamicEtaLinearFD2D20 kinematicPhiLinearFD2D20
);
etaNumeratorFD2D20 = (
  kinematicPhiLinearFD2D20 forcingDFD2D20
  -dynamicPhiLinearFD2D20 forcingKFD2D20
);
phiNumeratorFD2D20 = (
  dynamicEtaLinearFD2D20 forcingKFD2D20
  -kinematicEtaLinearFD2D20 forcingDFD2D20
);
etaGroupedExactFD2D20 = etaNumeratorFD2D20/determinantFD2D20;
phiGroupedExactFD2D20 = phiNumeratorFD2D20/determinantFD2D20;
etaOrderedExactFD2D20 = etaGroupedExactFD2D20/2;
phiOrderedExactFD2D20 = phiGroupedExactFD2D20/2;

(* Cross-multiplied original-equation recovery avoids expanding the
   nonzero resolvent denominator. *)
rKRecoveredFD2D20 = Expand[
  kinematicEtaLinearFD2D20 etaNumeratorFD2D20
  +kinematicPhiLinearFD2D20 phiNumeratorFD2D20
  +forcingKFD2D20 determinantFD2D20
];
rDRecoveredFD2D20 = Expand[
  dynamicEtaLinearFD2D20 etaNumeratorFD2D20
  +dynamicPhiLinearFD2D20 phiNumeratorFD2D20
  +forcingDFD2D20 determinantFD2D20
];
etaResidualFD2D20 = (
  dnoFD2D20[qDeltaFD2D20] rDFD2D20
  -I sigmaDeltaFD2D20 rKFD2D20
);
Print["Directional eta20 exact nonzero recovery assembled."];

unidirectionalForcingKFD2D20 = forcingKFD20 /. {
  gravity -> 1,
  hFD20 -> 1,
  k1FD20 -> q1FD2D20,
  k2FD20 -> q2FD2D20,
  omega1FD20 -> nu1FD2D20,
  omega2FD20 -> nu2FD2D20
};
unidirectionalForcingDFD2D20 = forcingDFD20 /. {
  gravity -> 1,
  hFD20 -> 1,
  k1FD20 -> q1FD2D20,
  k2FD20 -> q2FD2D20,
  omega1FD20 -> nu1FD2D20,
  omega2FD20 -> nu2FD2D20
};
unidirectionalRegressionFD2D20 = {
  Together[
    PowerExpand[
      (forcingKFD2D20 /. cosineFD2D20 -> 1)
      -unidirectionalForcingKFD2D20
    ]
  ],
  Together[
    PowerExpand[
      (forcingDFD2D20 /. cosineFD2D20 -> 1)
      -unidirectionalForcingDFD2D20
    ]
  ]
};
Print["Directional eta20 unidirectional regression evaluated."];

equalRadialKernelFD2D20 = FullSimplify[
  etaOrderedExactFD2D20 /. {
    q1FD2D20 -> carrierQFD2D20,
    q2FD2D20 -> carrierQFD2D20
  },
  Assumptions -> (
    carrierQFD2D20 > 0 && -1 < cosineFD2D20 < 1
  )
];
Print["Directional eta20 equal-radial kernel evaluated."];
equalRadialAngleLimitFD2D20 = FullSimplify[
  Limit[
    equalRadialKernelFD2D20,
    cosineFD2D20 -> 1,
    Direction -> "FromBelow"
  ],
  Assumptions -> carrierQFD2D20 > 0
];
radialOneSidedLimitFD2D20 = oneSidedEtaOrderedFD20 /. {
  carrierQFD20 -> carrierQFD2D20
};

strictMatrixFD2D20 = {
  {
    kinematicEtaLinearFD2D20,
    kinematicPhiLinearFD2D20
  },
  {
    dynamicEtaLinearFD2D20,
    dynamicPhiLinearFD2D20
  }
} /. {
  q1FD2D20 -> carrierQFD2D20,
  q2FD2D20 -> carrierQFD2D20,
  cosineFD2D20 -> 1
};
strictMatrixFD2D20 = FullSimplify[
  strictMatrixFD2D20,
  Assumptions -> carrierQFD2D20 > 0
];
strictRankFD2D20 = MatrixRank[strictMatrixFD2D20];

gateFD2D20[name_, residual_] := Module[{
  pass = And @@ (TrueQ[# === 0] & /@ Flatten[{residual}])
},
  Print[name <> " = " <> If[pass, "PASS", "FAIL"]];
  If[!pass, Print["  residual: ", InputForm[residual]]];
  pass
];
gateResultsFD2D20 = {
  gateFD2D20[
    "kinematic current-order linear block",
    FullSimplify[
      {
        kinematicEtaLinearFD2D20+I sigmaDeltaFD2D20,
        kinematicPhiLinearFD2D20+dnoFD2D20[qDeltaFD2D20]
      },
      Assumptions -> assumptionsFD2D20
    ]
  ],
  gateFD2D20[
    "dynamic current-order linear block",
    FullSimplify[
      {
        dynamicEtaLinearFD2D20-1,
        dynamicPhiLinearFD2D20+I sigmaDeltaFD2D20
      },
      Assumptions -> assumptionsFD2D20
    ]
  ],
  gateFD2D20[
    "determinant is finite-depth vector detuning",
    FullSimplify[
      determinantFD2D20
      -(dnoFD2D20[qDeltaFD2D20]-sigmaDeltaFD2D20^2),
      Assumptions -> assumptionsFD2D20
    ]
  ],
  gateFD2D20[
    "original kinematic residual after recovery",
    FullSimplify[rKRecoveredFD2D20, Assumptions -> assumptionsFD2D20]
  ],
  gateFD2D20[
    "original dynamic residual after recovery",
    FullSimplify[rDRecoveredFD2D20, Assumptions -> assumptionsFD2D20]
  ],
  gateFD2D20[
    "unidirectional exact-kernel regression",
    unidirectionalRegressionFD2D20
  ],
  gateFD2D20[
    "strict vector zero current-order rank is one",
    strictRankFD2D20-1
  ]
};
overallPassFD2D20 = And @@ gateResultsFD2D20;

exactTextFD2D20[expression_] := ToString[InputForm[expression]];
projectRootFD2D20 = DirectoryName[$InputFileName, 3];
artifactDirectoryFD2D20 = FileNameJoin[{projectRootFD2D20, "artifacts"}];
If[!DirectoryQ[artifactDirectoryFD2D20],
  CreateDirectory[
    artifactDirectoryFD2D20,
    CreateIntermediateDirectories -> True
  ]
];
artifactPathFD2D20 = FileNameJoin[{
  artifactDirectoryFD2D20,
  "order2_finite_depth_directional_eta20_exact.json"
}];
Export[
  artifactPathFD2D20,
  <|
    "schema_version" -> 1,
    "status" -> "directional-residual-only-before-candidate-grammar",
    "scope" ->
      "finite-depth forward directional strictly-nonzero eta20",
    "candidate_input_fields" -> {"eta11"},
    "root_result_used" -> False,
    "oracle_or_mf12_used" -> False,
    "candidate_selection_or_fitting" -> False,
    "domain" -> <|
      "radial" -> "q1>q2>0",
      "angular" -> "-1<cosine<=1",
      "strict_vector_zero" -> "separate rank-deficient sector"
    |>,
    "linear_system" -> <|
      "determinant" -> exactTextFD2D20[determinantFD2D20],
      "strict_zero_matrix" -> exactTextFD2D20[strictMatrixFD2D20],
      "strict_zero_rank" -> strictRankFD2D20
    |>,
    "exact_nonzero" -> <|
      "eta20_ordered" -> exactTextFD2D20[etaOrderedExactFD2D20],
      "Phi20_ordered" -> exactTextFD2D20[phiOrderedExactFD2D20]
    |>,
    "path_limits" -> <|
      "equal_radial_nonzero_angle" ->
        exactTextFD2D20[equalRadialKernelFD2D20],
      "equal_radial_angle_to_zero" ->
        exactTextFD2D20[equalRadialAngleLimitFD2D20],
      "unidirectional_radial_one_sided" ->
        exactTextFD2D20[radialOneSidedLimitFD2D20],
      "paths_identified" -> False
    |>,
    "gates" -> AssociationThread[
      {
        "kinematic current-order linear block",
        "dynamic current-order linear block",
        "determinant is finite-depth vector detuning",
        "original kinematic residual after recovery",
        "original dynamic residual after recovery",
        "unidirectional exact-kernel regression",
        "strict vector zero current-order rank is one"
      },
      gateResultsFD2D20
    ],
    "overall_pass" -> overallPassFD2D20
  |>,
  "RawJSON"
];

Print["=== Finite-depth directional eta20 Euler gate ==="];
Print["strict zero matrix = ", InputForm[strictMatrixFD2D20]];
Print[
  "equal-radial nonzero-angle kernel = ",
  InputForm[equalRadialKernelFD2D20]
];
Print[
  "equal-radial angle-to-zero limit = ",
  InputForm[equalRadialAngleLimitFD2D20]
];
Print[
  "unidirectional radial one-sided limit = ",
  InputForm[radialOneSidedLimitFD2D20]
];
Print["artifact = ", artifactPathFD2D20];
Print[
  "OVERALL_ORDER2_FINITE_DEPTH_DIRECTIONAL_ETA20_EULER = ",
  If[overallPassFD2D20, "PASS", "FAIL"]
];
If[!overallPassFD2D20, Exit[1]];
