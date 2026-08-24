(* ::Package:: *)
(* Exact leading-directional finite-depth order-three Euler residual.

   The original two-dimensional Euler free-surface residual is the governing
   source. Exact finite-depth directional pair fields are derived by
   order2_finite_depth_directional_euler.wl and remain inside this Wolfram
   derivation.

   Dimensionless convention: g = h = 1, q_i > 1/2,
     nu_i^2 = q_i Tanh[q_i].
   The radial magnitudes q1,q2,q3 remain exact. The directional jet is
   retained through rho^2, which is chi1.
*)

pairSourceFD3 = FileNameJoin[{
  DirectoryName[$InputFileName],
  "order2_finite_depth_directional_euler.wl"
}];
If[!FileExistsQ[pairSourceFD3],
  Print["FAIL: finite-depth directional order-two source not found."];
  Exit[1]
];
Get[pairSourceFD3];

ClearAll[
  amplitude1FD3, amplitude2FD3, amplitude3FD3,
  q1ParentFD3, q2ParentFD3, q3ParentFD3, qDiagonalFD3, rhoFD3,
  u1FD3, u2FD3, u3FD3, amplitudesFD3, magnitudesFD3,
  directionParametersFD3, pairIndicesFD3, nuFD3, dnoFD3,
  jetFD3, directionXFD3, directionYFD3, cosineFD3,
  kxFD3, kyFD3, pairMagnitudeFD3, outputKxFD3,
  outputKyFD3, outputMagnitudeFD3, pairInvariantFD3,
  pairChi1FD3, pairEtaJetFD3, pairPhiJetFD3, verticalFactorFD3,
  etaGroupedCoefficientFD3, phiGroupedCoefficientFD3,
  etaModesFD3, phiModesFD3, etaDerivativeFD3, phiDerivativeFD3,
  modeMasksFD3, etaPolyDerivativeFD3, phiPolyDerivativeFD3,
  chiFromExpressionFD3, chiToExpressionFD3, chiAddFD3, chiMultiplyFD3,
  polyAddFD3, polyMultiplyFD3, polyPowerFD3, polyEvaluateFD3,
  modalCoefficientFD3,
  rKDirectionalFD3, rDDirectionalFD3,
  kinematicEtaLinearFD3, kinematicPhiLinearFD3,
  kinematicForcingFD3, dynamicEtaLinearFD3,
  dynamicPhiLinearFD3, dynamicForcingFD3,
  elevationLinearFD3, elevationForcingFD3,
  etaDiagonalGroupedFD3, etaDiagonalOrderedFD3,
  diagonalRulesFD3, diagonalSolveFD3,
  pairCollinearCheckFD3, affinePotentialCheckFD3,
  certificationRulesFD3, certificationValuesFD3,
  projectRootFD3, artifactPathFD3, inputFormStringFD3,
  gateFD3, gateResultsFD3, overallPassFD3
];

amplitudesFD3 = {amplitude1FD3, amplitude2FD3, amplitude3FD3};
magnitudesFD3 = {
  q1ParentFD3,
  q2ParentFD3,
  q3ParentFD3
};
directionParametersFD3 = {u1FD3, u2FD3, u3FD3};
pairIndicesFD3 = {{1, 2}, {2, 3}, {3, 1}};

$Assumptions = (
  q1ParentFD3 > 1/2
  && q2ParentFD3 > 1/2
  && q3ParentFD3 > 1/2
);
nuFD3[q_] := Sqrt[q Tanh[q]];
dnoFD3[q_] := q Tanh[q];
jetFD3[expression_] := Normal[
  Series[expression, {rhoFD3, 0, 2}]
];

directionXFD3[index_] := (
  1 - 2 rhoFD3^2 directionParametersFD3[[index]]^2
);
directionYFD3[index_] := (
  2 rhoFD3 directionParametersFD3[[index]]
);
cosineFD3[first_, second_] := (
  1 - 2 rhoFD3^2 (
    directionParametersFD3[[first]]
    - directionParametersFD3[[second]]
  )^2
);
kxFD3[index_] := (
  magnitudesFD3[[index]] directionXFD3[index]
);
kyFD3[index_] := (
  magnitudesFD3[[index]] directionYFD3[index]
);
pairMagnitudeFD3[first_, second_] := pairMagnitudeFD3[first, second] =
  With[
    {
      qa = magnitudesFD3[[first]],
      qb = magnitudesFD3[[second]],
      du = directionParametersFD3[[first]]
        - directionParametersFD3[[second]]
    },
    qa + qb - 2 rhoFD3^2 qa qb du^2/(qa + qb)
  ];
outputKxFD3 = Total[kxFD3 /@ Range[3]];
outputKyFD3 = Total[kyFD3 /@ Range[3]];
outputMagnitudeFD3 = With[
  {
    qsum = Total[magnitudesFD3],
    weightedU = Total[
      magnitudesFD3 directionParametersFD3
    ],
    weightedU2 = Total[
      magnitudesFD3 directionParametersFD3^2
    ]
  },
  qsum + 2 rhoFD3^2 (
    weightedU^2/qsum - weightedU2
  )
];

pairInvariantFD3[type_, qa_, qb_, dot_, qout_] := Module[
  {nua, nub, nusum, sourceK, sourceD, detuning},
  nua = nuFD3[qa];
  nub = nuFD3[qb];
  nusum = nua + nub;
  sourceK = (qa^2 + dot)/nua + (qb^2 + dot)/nub;
  sourceD = nua^2 + nua nub + nub^2 - dot/(nua nub);
  detuning = dnoFD3[qout] - nusum^2;
  If[
    type === "eta",
    (dnoFD3[qout] sourceD - nusum sourceK)/detuning,
    -I (nusum sourceD - sourceK)/detuning
  ]
];
pairChi1FD3[type_, first_, second_] := Module[
  {
    qa, qb, du, dotDummy, outputDummy, expression,
    dotBase, outputBase, dotSecond, outputSecond
  },
  qa = magnitudesFD3[[first]];
  qb = magnitudesFD3[[second]];
  du = directionParametersFD3[[first]]
    - directionParametersFD3[[second]];
  dotBase = qa qb;
  outputBase = qa + qb;
  dotSecond = -2 qa qb du^2;
  outputSecond = -2 qa qb du^2/(qa + qb);
  expression = pairInvariantFD3[
    type, qa, qb, dotDummy, outputDummy
  ];
  (
    expression /. {
      dotDummy -> dotBase,
      outputDummy -> outputBase
    }
  ) + rhoFD3^2 (
    dotSecond (
      D[expression, dotDummy] /. {
        dotDummy -> dotBase,
        outputDummy -> outputBase
      }
    )
    + outputSecond (
      D[expression, outputDummy] /. {
        dotDummy -> dotBase,
        outputDummy -> outputBase
      }
    )
  )
];
pairEtaJetFD3[first_, second_] := pairEtaJetFD3[first, second] =
  pairChi1FD3["eta", first, second];
pairPhiJetFD3[first_, second_] := pairPhiJetFD3[first, second] =
  pairChi1FD3["phi", first, second];

verticalFactorFD3[q_, order_Integer] := Module[
  {q0, q2},
  q0 = Coefficient[q, rhoFD3, 0];
  q2 = Coefficient[q, rhoFD3, 2];
  If[
    EvenQ[order],
    q0^order
      + rhoFD3^2 q2 order q0^(order - 1),
    q0^order Tanh[q0]
      + rhoFD3^2 q2 (
        order q0^(order - 1) Tanh[q0]
        + q0^order Sech[q0]^2
      )
  ]
];

(* A mode record is {coefficient,kx,ky,vertical-rate,frequency}. *)
etaModesFD3[1] = Table[
  {
    amplitudesFD3[[index]],
    kxFD3[index],
    kyFD3[index],
    0,
    nuFD3[magnitudesFD3[[index]]]
  },
  {index, 3}
];
etaModesFD3[2] = Map[
  Function[pair,
    With[{first = pair[[1]], second = pair[[2]]},
      {
        pairEtaJetFD3[first, second]
          amplitudesFD3[[first]] amplitudesFD3[[second]],
        kxFD3[first] + kxFD3[second],
        kyFD3[first] + kyFD3[second],
        0,
        nuFD3[magnitudesFD3[[first]]]
          + nuFD3[magnitudesFD3[[second]]]
      }
    ]
  ],
  pairIndicesFD3
];
etaModesFD3[3] = {
  {
    etaGroupedCoefficientFD3 Times @@ amplitudesFD3,
    outputKxFD3,
    outputKyFD3,
    0,
    Total[nuFD3 /@ magnitudesFD3]
  }
};

phiModesFD3[1] = Table[
  {
    -I amplitudesFD3[[index]]
      /nuFD3[magnitudesFD3[[index]]],
    kxFD3[index],
    kyFD3[index],
    magnitudesFD3[[index]],
    nuFD3[magnitudesFD3[[index]]]
  },
  {index, 3}
];
phiModesFD3[2] = Map[
  Function[pair,
    With[{first = pair[[1]], second = pair[[2]]},
      {
        pairPhiJetFD3[first, second]
          amplitudesFD3[[first]] amplitudesFD3[[second]],
        kxFD3[first] + kxFD3[second],
        kyFD3[first] + kyFD3[second],
        pairMagnitudeFD3[first, second],
        nuFD3[magnitudesFD3[[first]]]
          + nuFD3[magnitudesFD3[[second]]]
      }
    ]
  ],
  pairIndicesFD3
];
phiModesFD3[3] = {
  {
    phiGroupedCoefficientFD3 Times @@ amplitudesFD3,
    outputKxFD3,
    outputKyFD3,
    outputMagnitudeFD3,
    Total[nuFD3 /@ magnitudesFD3]
  }
};

etaDerivativeFD3[order_, nx_, ny_, nt_] := Total[
  Function[mode,
    mode[[1]]
      (I mode[[2]])^nx
      (I mode[[3]])^ny
      (-I mode[[5]])^nt
  ] /@ etaModesFD3[order]
];
phiDerivativeFD3[order_, nx_, ny_, nz_, nt_] := Total[
  Function[mode,
    mode[[1]]
      (I mode[[2]])^nx
      (I mode[[3]])^ny
      verticalFactorFD3[mode[[4]], nz]
      (-I mode[[5]])^nt
  ] /@ phiModesFD3[order]
];

(* Evaluate the generated Euler residual in an exact labelled-polynomial
   algebra. Bit masks 1,2,4 identify the three parents. Products with a
   repeated parent label cannot contribute to A1 A2 A3 and are discarded
   immediately, preventing expression swell without copying the residual. *)
modeMasksFD3[1] = {1, 2, 4};
modeMasksFD3[2] = {3, 6, 5};
modeMasksFD3[3] = {7};

chiFromExpressionFD3[expression_] := {
  expression /. rhoFD3 -> 0,
  D[expression, rhoFD3] /. rhoFD3 -> 0,
  D[expression, {rhoFD3, 2}]/2 /. rhoFD3 -> 0
};
chiToExpressionFD3[coefficients_List] := (
  coefficients[[1]]
  + rhoFD3 coefficients[[2]]
  + rhoFD3^2 coefficients[[3]]
);
chiAddFD3[coefficients__List] := Total[{coefficients}];
chiMultiplyFD3[first_List, second_List] := Table[
  Sum[
    first[[index + 1]] second[[order - index + 1]],
    {index, 0, order}
  ],
  {order, 0, 2}
];

etaPolyDerivativeFD3[order_, nx_, ny_, nt_] := Association[
  MapThread[
    Function[{mask, mode},
      mask -> chiFromExpressionFD3[
        (mode[[1]] /. Thread[amplitudesFD3 -> 1])
        (I mode[[2]])^nx
        (I mode[[3]])^ny
        (-I mode[[5]])^nt
      ]
    ],
    {modeMasksFD3[order], etaModesFD3[order]}
  ]
];
phiPolyDerivativeFD3[order_, nx_, ny_, nz_, nt_] := Association[
  MapThread[
    Function[{mask, mode},
      mask -> chiFromExpressionFD3[
        (mode[[1]] /. Thread[amplitudesFD3 -> 1])
        (I mode[[2]])^nx
        (I mode[[3]])^ny
        verticalFactorFD3[mode[[4]], nz]
        (-I mode[[5]])^nt
      ]
    ],
    {modeMasksFD3[order], phiModesFD3[order]}
  ]
];
polyAddFD3[polynomials__Association] := Association[
  Table[
    mask -> Apply[
      chiAddFD3,
      Lookup[#, mask, {0, 0, 0}] & /@ {polynomials}
    ],
    {mask, 0, 7}
  ]
];
polyMultiplyFD3[first_Association, second_Association] := Association[
  Table[
    mask -> Total[
      Flatten[
        Table[
          If[
            BitAnd[left, right] == 0
              && BitOr[left, right] == mask,
            chiMultiplyFD3[
              Lookup[first, left, {0, 0, 0}],
              Lookup[second, right, {0, 0, 0}]
            ],
            {0, 0, 0}
          ],
          {left, 0, 7},
          {right, 0, 7}
        ],
        1
      ]
    ],
    {mask, 0, 7}
  ]
];
polyPowerFD3[base_Association, exponent_Integer?NonNegative] :=
  Nest[
    polyMultiplyFD3[#, base] &,
    <|0 -> {1, 0, 0}|>,
    exponent
  ];

polyEvaluateFD3[
  HoldPattern[
    Derivative[nx_, ny_, nt_][etaG2[order_]][_, _, _]
  ]
] := etaPolyDerivativeFD3[order, nx, ny, nt];
polyEvaluateFD3[HoldPattern[etaG2[order_][_, _, _]]] :=
  etaPolyDerivativeFD3[order, 0, 0, 0];
polyEvaluateFD3[
  HoldPattern[
    Derivative[nx_, ny_, nz_, nt_][phiG2[order_]][_, _, _, _]
  ]
] := phiPolyDerivativeFD3[order, nx, ny, nz, nt];
polyEvaluateFD3[HoldPattern[phiG2[order_][_, _, _, _]]] :=
  phiPolyDerivativeFD3[order, 0, 0, 0, 0];
polyEvaluateFD3[expression_Plus] :=
  Apply[polyAddFD3, polyEvaluateFD3 /@ List @@ expression];
polyEvaluateFD3[expression_Times] :=
  Fold[
    polyMultiplyFD3,
    <|0 -> {1, 0, 0}|>,
    polyEvaluateFD3 /@ List @@ expression
  ];
polyEvaluateFD3[Power[base_, exponent_Integer?NonNegative]] :=
  polyPowerFD3[polyEvaluateFD3[base], exponent];
polyEvaluateFD3[expression_] :=
  <|0 -> chiFromExpressionFD3[expression]|>;

modalCoefficientFD3[expression_] := With[
  {coefficients = Lookup[
    polyEvaluateFD3[expression /. gravityG2 -> 1],
    7,
    {0, 0, 0}
  ]},
  coefficients[[1]]
    + rhoFD3 coefficients[[2]]
    + rhoFD3^2 coefficients[[3]]
];

Print["Extracting finite-depth directional order-three chi1 residual."];
rKDirectionalFD3 = modalCoefficientFD3[kinematicExpectedG2[[3]]];
Print["  kinematic residual extracted."];
rDDirectionalFD3 = modalCoefficientFD3[dynamicExpectedG2[[3]]];
Print["  dynamic residual extracted."];

kinematicEtaLinearFD3 = Coefficient[
  rKDirectionalFD3, etaGroupedCoefficientFD3
];
kinematicPhiLinearFD3 = Coefficient[
  rKDirectionalFD3, phiGroupedCoefficientFD3
];
kinematicForcingFD3 = rKDirectionalFD3 /. {
  etaGroupedCoefficientFD3 -> 0,
  phiGroupedCoefficientFD3 -> 0
};
dynamicEtaLinearFD3 = Coefficient[
  rDDirectionalFD3, etaGroupedCoefficientFD3
];
dynamicPhiLinearFD3 = Coefficient[
  rDDirectionalFD3, phiGroupedCoefficientFD3
];
dynamicForcingFD3 = rDDirectionalFD3 /. {
  etaGroupedCoefficientFD3 -> 0,
  phiGroupedCoefficientFD3 -> 0
};

(* Eliminate the current-order flat bulk potential exactly:
     R_E = B_K R_D - B_D R_K. *)
elevationLinearFD3 = chiToExpressionFD3[
  chiMultiplyFD3[
    chiFromExpressionFD3[kinematicPhiLinearFD3],
    chiFromExpressionFD3[dynamicEtaLinearFD3]
  ]
  - chiMultiplyFD3[
    chiFromExpressionFD3[dynamicPhiLinearFD3],
    chiFromExpressionFD3[kinematicEtaLinearFD3]
  ]
];
elevationForcingFD3 = chiToExpressionFD3[
  chiMultiplyFD3[
    chiFromExpressionFD3[kinematicPhiLinearFD3],
    chiFromExpressionFD3[dynamicForcingFD3]
  ]
  - chiMultiplyFD3[
    chiFromExpressionFD3[dynamicPhiLinearFD3],
    chiFromExpressionFD3[kinematicForcingFD3]
  ]
];
affinePotentialCheckFD3 = PossibleZeroQ[
  kinematicPhiLinearFD3 dynamicPhiLinearFD3
  - dynamicPhiLinearFD3 kinematicPhiLinearFD3
];

diagonalRulesFD3 = {
  q1ParentFD3 -> qDiagonalFD3,
  q2ParentFD3 -> qDiagonalFD3,
  q3ParentFD3 -> qDiagonalFD3,
  rhoFD3 -> 0
};
diagonalSolveFD3 = Solve[
  Numerator[
    Together[
      (
        elevationLinearFD3 etaGroupedCoefficientFD3
        + elevationForcingFD3
      ) /. diagonalRulesFD3
    ]
  ] == 0,
  etaGroupedCoefficientFD3
];
If[Length[diagonalSolveFD3] =!= 1,
  Print["FAIL: finite-depth directional diagonal was not uniquely solved."];
  Exit[1]
];
etaDiagonalGroupedFD3 = FullSimplify[
  etaGroupedCoefficientFD3 /. First[diagonalSolveFD3],
  Assumptions -> qDiagonalFD3 > 1/2
];
etaDiagonalOrderedFD3 = etaDiagonalGroupedFD3/24;

pairCollinearCheckFD3 = Together[
  (pairEtaJetFD3[1, 2] /. rhoFD3 -> 0)
  - pairInvariantFD3[
      "eta",
      q1ParentFD3,
      q2ParentFD3,
      q1ParentFD3 q2ParentFD3,
      q1ParentFD3 + q2ParentFD3
    ]
];

certificationRulesFD3 = {
  q1ParentFD3 -> 9/10,
  q2ParentFD3 -> 21/20,
  q3ParentFD3 -> 6/5,
  rhoFD3 -> 1/25,
  u1FD3 -> -1,
  u2FD3 -> 1/3,
  u3FD3 -> 4/5,
  etaGroupedCoefficientFD3 -> 2/3,
  phiGroupedCoefficientFD3 -> I/7
};
certificationValuesFD3 = N[
  {
    elevationLinearFD3,
    elevationForcingFD3,
    kinematicEtaLinearFD3,
    kinematicPhiLinearFD3,
    dynamicEtaLinearFD3,
    dynamicPhiLinearFD3
  } /. certificationRulesFD3,
  30
];

projectRootFD3 = DirectoryName[$InputFileName, 3];
artifactPathFD3 = FileNameJoin[{
  projectRootFD3,
  "artifacts",
  "order3_finite_depth_directional_euler_result.json"
}];
inputFormStringFD3[expression_] := ToString[InputForm[expression]];
Export[
  artifactPathFD3,
  <|
    "schema_version" -> 1,
    "status" -> "exact_directional_chi1_elevation_residual",
    "scope" ->
      "finite-depth forward-directional positive pure-sum eta33 q_i>0.5",
    "candidate_input_fields" -> {"eta11"},
    "lower_order_fields_at_candidate_boundary" -> False,
    "oracle_imported" -> False,
    "direction_parameterization" ->
      "e(t)=((1-t^2)/(1+t^2),2t/(1+t^2)), t_i=rho u_i",
    "jet_orders" -> <|"rho" -> 2, "radial_q" -> "exact"|>,
    "elimination_identity" -> "R_E=B_K R_D-B_D R_K",
    "elevation_linear" -> inputFormStringFD3[elevationLinearFD3],
    "elevation_forcing" -> inputFormStringFD3[elevationForcingFD3],
    "eta3_grouped_diagonal" ->
      inputFormStringFD3[etaDiagonalGroupedFD3],
    "eta3_ordered_diagonal" ->
      inputFormStringFD3[etaDiagonalOrderedFD3],
    "affine_leaf_counts" -> <|
      "elevation_linear" -> LeafCount[elevationLinearFD3],
      "elevation_forcing" -> LeafCount[elevationForcingFD3]
    |>,
    "certification_values" ->
      inputFormStringFD3[certificationValuesFD3]
  |>,
  "RawJSON"
];

gateFD3[name_, residual_] := Module[
  {pass = TrueQ[PossibleZeroQ[residual]]},
  Print[name <> " = " <> If[pass, "PASS", "FAIL"]];
  If[!pass, Print["  residual: ", InputForm[residual]]];
  pass
];
gateResultsFD3 = {
  TrueQ[affinePotentialCheckFD3],
  gateFD3["finite-depth pair jet has the exact collinear reduction",
    pairCollinearCheckFD3],
  And @@ (NumericQ[#] & /@ certificationValuesFD3)
};
overallPassFD3 = And @@ gateResultsFD3;
Print["=== Finite-depth directional eta33 Euler residual ==="];
Print["diagonal ordered eta33 = ", InputForm[etaDiagonalOrderedFD3]];
Print["affine leaf counts = ",
  InputForm[{
    LeafCount[elevationLinearFD3],
    LeafCount[elevationForcingFD3]
  }]];
Print["artifact = ", artifactPathFD3];
Print[
  "OVERALL_ORDER3_FINITE_DEPTH_DIRECTIONAL_EULER = ",
  If[overallPassFD3, "PASS", "FAIL"]
];
If[!overallPassFD3, Exit[1]];
