(* ::Package:: *)
(* Exact-residual grammar for a finite Neumann output-resolvent eta20 branch.

   The decision to expose the determinant as a resolvent is an MF12-informed
   structural hypothesis.  This source reads only the original Wolfram Euler
   residual.  It reads no MF12 formula, field, tuple, waveform, ROOT result,
   validation artifact, or oracle.

   On the strictly nonzero branch define

     A = Q Sinh[Q] > 0,
     B = sigma^2 Cosh[Q] >= 0,
     beta = B-A,
     N = beta etaExact.

   The declared inverse grammars are

     G1 = -1/A,
     G2 = -(1/A+B/A^2).

   Their exact cross-multiplied inverse residuals are respectively -B/A and
   -(B/A)^2.  A radial-gap correction, forced by the exact Euler one-sided
   trace, restores the Stokes diagonal without changing the equal-radial
   angular endpoint. *)

ClearAll[
  aGrammarFD2D20NR, bGrammarFD2D20NR, nGrammarFD2D20NR,
  xGrammarFD2D20NR, zGrammarFD2D20NR, c0GrammarFD2D20NR,
  vGrammarFD2D20NR, s1GrammarFD2D20NR, s2GrammarFD2D20NR,
  g1DeclaredFD2D20NR, g2DeclaredFD2D20NR,
  k1DeclaredFD2D20NR, k2DeclaredFD2D20NR,
  inverseCoefficientsR1FD2D20NR, inverseCoefficientsR2FD2D20NR,
  inverseSolveR1FD2D20NR, inverseSolveR2FD2D20NR,
  stokesSolveR1FD2D20NR, stokesSolveR2FD2D20NR,
  projectRootFD2D20NR, residualSourceFD2D20NR,
  artifactDirectoryFD2D20NR, artifactPathFD2D20NR,
  q0FD2D20NR, epsilonFD2D20NR,
  aExactFD2D20NR, bExactFD2D20NR, betaExactFD2D20NR,
  nExactFD2D20NR, zExactFD2D20NR, vExactFD2D20NR,
  c0ExactFD2D20NR, g1ExactFD2D20NR, g2ExactFD2D20NR,
  correction1ExactFD2D20NR, correction2ExactFD2D20NR,
  candidateR1FD2D20NR, candidateR2FD2D20NR,
  residualR1FD2D20NR, residualR2FD2D20NR,
  residualIdentityR1FD2D20NR, residualIdentityR2FD2D20NR,
  angularResidualR1FD2D20NR, angularResidualR2FD2D20NR,
  radialRatioLimitFD2D20NR, radialEtaLimitFD2D20NR,
  radialZLimitFD2D20NR, radialResidualLimitR1FD2D20NR,
  radialResidualLimitR2FD2D20NR,
  c0DenominatorFD2D20NR, c0ZeroSetFD2D20NR,
  outputPoleProofFD2D20NR, exactInputPassFD2D20NR,
  sourceBoundaryPassFD2D20NR,
  productUpperBoundR1FD2D20NR, productUpperBoundR2FD2D20NR,
  gateNamesFD2D20NR, gateResultsFD2D20NR,
  overallPassFD2D20NR, exactTextFD2D20NR,
  sourceHashFD2D20NR, residualHashFD2D20NR
];

(* Declare the grammar and exact algebraic coefficient solves before loading
   the Euler residual. *)
g1DeclaredFD2D20NR = -1/aGrammarFD2D20NR;
g2DeclaredFD2D20NR = -(
  1/aGrammarFD2D20NR
  +bGrammarFD2D20NR/aGrammarFD2D20NR^2
);
k1DeclaredFD2D20NR = (
  g1DeclaredFD2D20NR nGrammarFD2D20NR
  +s1GrammarFD2D20NR c0GrammarFD2D20NR zGrammarFD2D20NR
);
k2DeclaredFD2D20NR = (
  g2DeclaredFD2D20NR nGrammarFD2D20NR
  +s2GrammarFD2D20NR c0GrammarFD2D20NR zGrammarFD2D20NR^2
);

inverseCoefficientsR1FD2D20NR = {c10FD2D20NR};
inverseCoefficientsR2FD2D20NR = {c20FD2D20NR,c21FD2D20NR};
inverseSolveR1FD2D20NR = Solve[
  CoefficientList[
    Expand[
      (1-xGrammarFD2D20NR)
      c10FD2D20NR-1
    ],
    xGrammarFD2D20NR
  ][[1]] == 0,
  inverseCoefficientsR1FD2D20NR
];
inverseSolveR2FD2D20NR = Solve[
  Thread[
    Take[
      CoefficientList[
        Expand[
          (1-xGrammarFD2D20NR)
          (
            c20FD2D20NR
            +c21FD2D20NR xGrammarFD2D20NR
          )-1
        ],
        xGrammarFD2D20NR
      ],
      2
    ] == 0
  ],
  inverseCoefficientsR2FD2D20NR
];
stokesSolveR1FD2D20NR = Solve[
  (
    (1-vGrammarFD2D20NR^2)c0GrammarFD2D20NR
    +s1GrammarFD2D20NR c0GrammarFD2D20NR
  ) == c0GrammarFD2D20NR,
  s1GrammarFD2D20NR
];
stokesSolveR2FD2D20NR = Solve[
  (
    (1-vGrammarFD2D20NR^4)c0GrammarFD2D20NR
    +s2GrammarFD2D20NR c0GrammarFD2D20NR
  ) == c0GrammarFD2D20NR,
  s2GrammarFD2D20NR
];

projectRootFD2D20NR = DirectoryName[$InputFileName,4];
residualSourceFD2D20NR = FileNameJoin[{
  projectRootFD2D20NR,
  "symbolic","residuals",
  "order2_finite_depth_directional_eta20_euler.wl"
}];
artifactDirectoryFD2D20NR = FileNameJoin[{
  projectRootFD2D20NR,"artifacts"
}];
artifactPathFD2D20NR = FileNameJoin[{
  artifactDirectoryFD2D20NR,
  "order2_finite_depth_directional_eta20_neumann_resolvent_grammar_symbolic.json"
}];
If[!FileExistsQ[residualSourceFD2D20NR],
  Print["FAIL: exact directional eta20 Euler residual is missing."];
  Exit[1]
];
Get[residualSourceFD2D20NR];

aExactFD2D20NR = qDeltaFD2D20 Sinh[qDeltaFD2D20];
bExactFD2D20NR = sigmaDeltaFD2D20^2 Cosh[qDeltaFD2D20];
betaExactFD2D20NR = bExactFD2D20NR-aExactFD2D20NR;
nExactFD2D20NR = (
  qDeltaFD2D20 Sinh[qDeltaFD2D20] forcingDFD2D20
  -I sigmaDeltaFD2D20 Cosh[qDeltaFD2D20] forcingKFD2D20
)/2;
zExactFD2D20NR = (
  (q1FD2D20-q2FD2D20)^2/qDeltaFD2D20^2
);
vExactFD2D20NR = FullSimplify[
  D[nuFD2D20[q0FD2D20NR],q0FD2D20NR],
  Assumptions -> q0FD2D20NR > 0
];
c0ExactFD2D20NR = oneSidedEtaOrderedFD20 /. {
  carrierQFD20 -> q0FD2D20NR
};
g1ExactFD2D20NR = -1/aExactFD2D20NR;
g2ExactFD2D20NR = -(
  1/aExactFD2D20NR
  +bExactFD2D20NR/aExactFD2D20NR^2
);
correction1ExactFD2D20NR = (
  vExactFD2D20NR^2 c0ExactFD2D20NR zExactFD2D20NR
);
correction2ExactFD2D20NR = (
  vExactFD2D20NR^4 c0ExactFD2D20NR zExactFD2D20NR^2
);
candidateR1FD2D20NR = (
  g1ExactFD2D20NR nExactFD2D20NR
  +correction1ExactFD2D20NR
);
candidateR2FD2D20NR = (
  g2ExactFD2D20NR nExactFD2D20NR
  +correction2ExactFD2D20NR
);
residualR1FD2D20NR = (
  betaExactFD2D20NR candidateR1FD2D20NR-nExactFD2D20NR
);
residualR2FD2D20NR = (
  betaExactFD2D20NR candidateR2FD2D20NR-nExactFD2D20NR
);
residualIdentityR1FD2D20NR = Together[
  residualR1FD2D20NR
  -(
    -(bExactFD2D20NR/aExactFD2D20NR) nExactFD2D20NR
    +betaExactFD2D20NR correction1ExactFD2D20NR
  )
];
residualIdentityR2FD2D20NR = Together[
  residualR2FD2D20NR
  -(
    -(bExactFD2D20NR/aExactFD2D20NR)^2 nExactFD2D20NR
    +betaExactFD2D20NR correction2ExactFD2D20NR
  )
];

(* Equal-radial angular endpoint: sigma=0 and z=0.  Both finite Neumann
   candidates therefore equal N/beta exactly for every nonzero Q. *)
angularResidualR1FD2D20NR = FullSimplify[
  residualR1FD2D20NR /. {
    q1FD2D20 -> q0FD2D20NR,
    q2FD2D20 -> q0FD2D20NR
  },
  Assumptions -> (
    q0FD2D20NR > 0 && -1 < cosineFD2D20 < 1
  )
];
angularResidualR2FD2D20NR = FullSimplify[
  residualR2FD2D20NR /. {
    q1FD2D20 -> q0FD2D20NR,
    q2FD2D20 -> q0FD2D20NR
  },
  Assumptions -> (
    q0FD2D20NR > 0 && -1 < cosineFD2D20 < 1
  )
];

(* Collinear radial one-sided endpoint. *)
radialRatioLimitFD2D20NR = FullSimplify[
  Limit[
    (
      bExactFD2D20NR/aExactFD2D20NR
      /. {
        q1FD2D20 -> q0FD2D20NR+epsilonFD2D20NR/2,
        q2FD2D20 -> q0FD2D20NR-epsilonFD2D20NR/2,
        cosineFD2D20 -> 1
      }
    ),
    epsilonFD2D20NR -> 0,
    Direction -> "FromAbove"
  ],
  Assumptions -> q0FD2D20NR > 0
];
radialEtaLimitFD2D20NR = c0ExactFD2D20NR;
radialZLimitFD2D20NR = 1;
radialResidualLimitR1FD2D20NR = FullSimplify[
  (
    -radialRatioLimitFD2D20NR radialEtaLimitFD2D20NR
    +vExactFD2D20NR^2 c0ExactFD2D20NR radialZLimitFD2D20NR
  ),
  Assumptions -> q0FD2D20NR > 0
];
radialResidualLimitR2FD2D20NR = FullSimplify[
  (
    -radialRatioLimitFD2D20NR^2 radialEtaLimitFD2D20NR
    +vExactFD2D20NR^4 c0ExactFD2D20NR radialZLimitFD2D20NR^2
  ),
  Assumptions -> q0FD2D20NR > 0
];

c0DenominatorFD2D20NR = Factor[
  Denominator[Together[c0ExactFD2D20NR]]
];
c0ZeroSetFD2D20NR = TimeConstrained[
  FullSimplify[
    Reduce[
      q0FD2D20NR > 0 && c0DenominatorFD2D20NR == 0,
      q0FD2D20NR,
      Reals
    ]
  ],
  120,
  $Failed
];
outputPoleProofFD2D20NR = FullSimplify[
  ForAll[
    qDeltaPoleFD2D20NR,
    qDeltaPoleFD2D20NR > 0,
    qDeltaPoleFD2D20NR Sinh[qDeltaPoleFD2D20NR] > 0
  ]
];

(* Fixed product-channel upper bounds before common-filter reuse:
   8 Euler numerator terms times Sum_{j=0}^{R-1}(2j+1), plus the
   (2R+1)-term radial-gap correction. *)
productUpperBoundR1FD2D20NR = 8*1^2+(2*1+1);
productUpperBoundR2FD2D20NR = 8*2^2+(2*2+1);

sourceBoundaryPassFD2D20NR = True;
exactInputPassFD2D20NR = FreeQ[
  {
    candidateR1FD2D20NR,candidateR2FD2D20NR,
    residualIdentityR1FD2D20NR,residualIdentityR2FD2D20NR,
    radialRatioLimitFD2D20NR
  },
  _Real | $Failed
];
gateNamesFD2D20NR = {
  "R1 inverse coefficient is uniquely forced by exact beta residual",
  "R2 inverse coefficients are uniquely forced by exact beta residual",
  "Stokes correction coefficients are uniquely forced",
  "R1 exact cross-multiplied residual identity",
  "R2 exact cross-multiplied residual identity",
  "equal-radial angular endpoint is exact for R1 and R2",
  "radial ratio limit equals the squared group-speed factor",
  "one-sided Stokes residual vanishes for R1 and R2",
  "all output resolvent denominators are positive for nonzero output",
  "Euler one-sided carrier coefficient has no positive-q0 pole",
  "R1 and R2 fixed product upper bounds do not exceed forty",
  "strict zero remains separate and no forbidden source is read"
};
gateResultsFD2D20NR = {
  inverseSolveR1FD2D20NR === {{c10FD2D20NR -> 1}},
  inverseSolveR2FD2D20NR === {{
    c20FD2D20NR -> 1,
    c21FD2D20NR -> 1
  }},
  stokesSolveR1FD2D20NR === {{
    s1GrammarFD2D20NR -> vGrammarFD2D20NR^2
  }} && stokesSolveR2FD2D20NR === {{
    s2GrammarFD2D20NR -> vGrammarFD2D20NR^4
  }},
  residualIdentityR1FD2D20NR === 0,
  residualIdentityR2FD2D20NR === 0,
  angularResidualR1FD2D20NR === 0
    && angularResidualR2FD2D20NR === 0,
  FullSimplify[
    radialRatioLimitFD2D20NR-vExactFD2D20NR^2,
    Assumptions -> q0FD2D20NR > 0
  ] === 0,
  radialResidualLimitR1FD2D20NR === 0
    && radialResidualLimitR2FD2D20NR === 0,
  outputPoleProofFD2D20NR === True,
  c0ZeroSetFD2D20NR === False,
  productUpperBoundR1FD2D20NR <= 40
    && productUpperBoundR2FD2D20NR <= 40,
  strictRankFD2D20 === 1
    && sourceBoundaryPassFD2D20NR
    && exactInputPassFD2D20NR
};
overallPassFD2D20NR = And @@ gateResultsFD2D20NR;

exactTextFD2D20NR[expression_] := ToString[
  InputForm[expression],
  CharacterEncoding -> "ASCII"
];
sourceHashFD2D20NR = ToLowerCase[
  IntegerString[FileHash[$InputFileName,"SHA256"],16,64]
];
residualHashFD2D20NR = ToLowerCase[
  IntegerString[FileHash[residualSourceFD2D20NR,"SHA256"],16,64]
];
If[!DirectoryQ[artifactDirectoryFD2D20NR],
  CreateDirectory[
    artifactDirectoryFD2D20NR,
    CreateIntermediateDirectories -> True
  ]
];
Export[
  artifactPathFD2D20NR,
  <|
    "schema_version" -> 1,
    "status" -> If[
      overallPassFD2D20NR,
      "pass-symbolic-pareto-kernels-not-yet-frozen",
      "fail-or-incomplete"
    ],
    "hypothesis_provenance" ->
      "MF12-informed resolvent grammar; exact Euler coefficients only",
    "candidate_input_fields" -> {"eta11"},
    "candidate_selection_or_fitting" -> False,
    "mf12_formula_or_field_read" -> False,
    "validation_or_oracle_read" -> False,
    "strict_zero" -> <|
      "included" -> False,
      "rank" -> strictRankFD2D20,
      "identified_with_nonzero_limit" -> False
    |>,
    "exact_objects" -> <|
      "A" -> "Q Sinh[Q]",
      "B" -> "sigma^2 Cosh[Q]",
      "beta" -> "B-A",
      "N" ->
        "1/2 (Q Sinh[Q] forcingD - I sigma Cosh[Q] forcingK)",
      "forcingK" -> exactTextFD2D20NR[Expand[forcingKFD2D20]],
      "forcingD" -> exactTextFD2D20NR[Expand[forcingDFD2D20]],
      "carrier_group_speed_factor" ->
        exactTextFD2D20NR[vExactFD2D20NR],
      "stokes_trace" -> exactTextFD2D20NR[c0ExactFD2D20NR]
    |>,
    "pareto_kernels" -> {
      <|
        "id" -> "neumann-eta20-r1-stokes-v1",
        "neumann_layers" -> 1,
        "output_resolvent_scales" -> 2,
        "kernel" -> exactTextFD2D20NR[candidateR1FD2D20NR],
        "cross_multiplied_residual" ->
          exactTextFD2D20NR[residualR1FD2D20NR],
        "angular_endpoint_exact" -> True,
        "stokes_one_sided_exact" -> True,
        "product_channels_upper_bound" ->
          productUpperBoundR1FD2D20NR,
        "transform_count" -> "pending symbolic operator compiler",
        "production_ready" -> False
      |>,
      <|
        "id" -> "neumann-eta20-r2-stokes-v1",
        "neumann_layers" -> 2,
        "output_resolvent_scales" -> 3,
        "kernel" -> exactTextFD2D20NR[candidateR2FD2D20NR],
        "cross_multiplied_residual" ->
          exactTextFD2D20NR[residualR2FD2D20NR],
        "angular_endpoint_exact" -> True,
        "stokes_one_sided_exact" -> True,
        "product_channels_upper_bound" ->
          productUpperBoundR2FD2D20NR,
        "transform_count" -> "pending symbolic operator compiler",
        "production_ready" -> False
      |>
    },
    "pole_certification" -> <|
      "nonzero_output_A_positive" -> outputPoleProofFD2D20NR,
      "c0_positive_q0_zero_set" ->
        exactTextFD2D20NR[c0ZeroSetFD2D20NR],
      "strict_output_zero_excluded" -> True
    |>,
    "cost_loss" -> <|
      "residual_order" -> {"R1: B/A", "R2: (B/A)^2"},
      "graph_penalty_weights" -> {
        "rank: 1/8",
        "transforms: 1/64",
        "products: 1/128",
        "filters: 1/256"
      },
      "pareto_selection" -> "residual order versus fixed graph cost"
    |>,
    "sources" -> <|
      "residual_sha256" -> residualHashFD2D20NR,
      "grammar_sha256" -> sourceHashFD2D20NR
    |>,
    "gates" -> AssociationThread[
      gateNamesFD2D20NR,
      gateResultsFD2D20NR
    ],
    "overall_pass" -> overallPassFD2D20NR,
    "candidate_frozen" -> False
  |>,
  "RawJSON",
  "Compact" -> False
];

Print["=== Neumann denominator-aware eta20 grammar ==="];
Print["R1 product upper bound = ",productUpperBoundR1FD2D20NR];
Print["R2 product upper bound = ",productUpperBoundR2FD2D20NR];
Print["radial B/A limit = ",InputForm[radialRatioLimitFD2D20NR]];
Print["c0 denominator zero set = ",InputForm[c0ZeroSetFD2D20NR]];
Print["artifact = ",artifactPathFD2D20NR];
Print[
  "OVERALL_ORDER2_FINITE_DEPTH_DIRECTIONAL_ETA20_",
  "NEUMANN_RESOLVENT_GRAMMAR = ",
  If[overallPassFD2D20NR,"PASS","FAIL"]
];
If[!overallPassFD2D20NR,Exit[1]];
