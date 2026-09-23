(* ::Package:: *)
(* Exact two-horizontal-dimensional extraction of the epsilon^1--epsilon^4
   residuals of the Euler free-surface boundary conditions.

   The unknowns phiG2[n][x,y,z,t] are BULK velocity-potential coefficients.
   Every surface quantity is Taylor evaluated from z = eta to z = 0 before
   perturbation coefficients are collected.  The certified order-one through
   order-three source is imported only for regression.  No order-four kernel,
   directional candidate, tuple sample, external field, or validation oracle is
   imported.
*)

lowerOrder2DScriptG24 = FileNameJoin[{
  DirectoryName[$InputFileName],
  "order1_to_order3_symbolic_expansion_2d.wl"
}];
If[!FileExistsQ[lowerOrder2DScriptG24],
  Print["FAIL: lower-order 2D residual source not found: ",
    lowerOrder2DScriptG24];
  Exit[1]
];
Get[lowerOrder2DScriptG24];

kinematicLowerReferenceG24 = kinematicActualG2;
dynamicLowerReferenceG24 = dynamicActualG2;
combinedLowerReferenceG24 = combinedActualG2;

ClearAll[
  epsilonG2, xG2, yG2, zG2, tG2, gravityG2, etaG2, phiG2,
  etaSeriesG2, phiSeriesG2, truncateG2, surfaceTaylorG2,
  phiXSurfaceG2, phiYSurfaceG2, phiZSurfaceG2, phiTSurfaceG2,
  kinematicSeriesG2, dynamicSeriesG2, qBulkG2, combinedBulkG2,
  combinedSeriesG2, coefficientG2, atZeroG2, horizontalDotG2,
  kinematicActualG2, dynamicActualG2, combinedActualG2,
  lowerRegressionChecksG24, currentOrderLinearChecksG24,
  eta4G2, phi4G2, lowerZeroRulesG24,
  kinematicLinearOrder4G2, dynamicLinearOrder4G2,
  combinedLinearOrder4G2, kinematicForcingOrder4G2,
  dynamicForcingOrder4G2, combinedForcingOrder4G2,
  forcingSplitChecksG24, forcingBoundaryChecksG24,
  gateG24, gateResultsG24, overallPassG24
];

etaSeriesG2 = Sum[
  epsilonG2^n etaG2[n][xG2, yG2, tG2],
  {n, 1, 4}
];
phiSeriesG2 = Sum[
  epsilonG2^n phiG2[n][xG2, yG2, zG2, tG2],
  {n, 1, 4}
];

truncateG2[expression_] := Normal[
  Series[expression, {epsilonG2, 0, 4}]
];

surfaceTaylorG2[expression_] := truncateG2[
  Sum[
    etaSeriesG2^m/m!
      (D[expression, {zG2, m}] /. zG2 -> 0),
    {m, 0, 4}
  ]
];

horizontalDotG2[first_, second_] := (
  D[first, xG2] D[second, xG2]
  + D[first, yG2] D[second, yG2]
);

phiXSurfaceG2 = surfaceTaylorG2[D[phiSeriesG2, xG2]];
phiYSurfaceG2 = surfaceTaylorG2[D[phiSeriesG2, yG2]];
phiZSurfaceG2 = surfaceTaylorG2[D[phiSeriesG2, zG2]];
phiTSurfaceG2 = surfaceTaylorG2[D[phiSeriesG2, tG2]];

kinematicSeriesG2 = truncateG2[
  D[etaSeriesG2, tG2]
  + phiXSurfaceG2 D[etaSeriesG2, xG2]
  + phiYSurfaceG2 D[etaSeriesG2, yG2]
  - phiZSurfaceG2
];

dynamicSeriesG2 = truncateG2[
  phiTSurfaceG2
  + 1/2 (
      phiXSurfaceG2^2
      + phiYSurfaceG2^2
      + phiZSurfaceG2^2
    )
  + gravityG2 etaSeriesG2
];

qBulkG2 = (
  D[phiSeriesG2, xG2]^2
  + D[phiSeriesG2, yG2]^2
  + D[phiSeriesG2, zG2]^2
);
combinedBulkG2 = (
  D[phiSeriesG2, {tG2, 2}]
  + gravityG2 D[phiSeriesG2, zG2]
  + D[qBulkG2, tG2]
  + 1/2 (
      D[phiSeriesG2, xG2] D[qBulkG2, xG2]
      + D[phiSeriesG2, yG2] D[qBulkG2, yG2]
      + D[phiSeriesG2, zG2] D[qBulkG2, zG2]
    )
);
combinedSeriesG2 = surfaceTaylorG2[combinedBulkG2];

coefficientG2[series_, order_Integer] := Expand[
  Coefficient[series, epsilonG2, order]
];
atZeroG2[expression_] := expression /. zG2 -> 0;

kinematicActualG2 = Table[
  coefficientG2[kinematicSeriesG2, n],
  {n, 1, 4}
];
dynamicActualG2 = Table[
  coefficientG2[dynamicSeriesG2, n],
  {n, 1, 4}
];
combinedActualG2 = Table[
  coefficientG2[combinedSeriesG2, n],
  {n, 1, 4}
];

lowerRegressionChecksG24 = Join[
  Thread[Take[kinematicActualG2, 3] - kinematicLowerReferenceG24],
  Thread[Take[dynamicActualG2, 3] - dynamicLowerReferenceG24],
  Thread[Take[combinedActualG2, 3] - combinedLowerReferenceG24]
];

eta4G2 = etaG2[4][xG2, yG2, tG2];
phi4G2 = phiG2[4][xG2, yG2, zG2, tG2];

kinematicLinearOrder4G2 = atZeroG2[
  D[eta4G2, tG2] - D[phi4G2, zG2]
];
dynamicLinearOrder4G2 = atZeroG2[
  D[phi4G2, tG2] + gravityG2 eta4G2
];
combinedLinearOrder4G2 = atZeroG2[
  D[phi4G2, {tG2, 2}] + gravityG2 D[phi4G2, zG2]
];

kinematicForcingOrder4G2 = Expand[
  kinematicActualG2[[4]] - kinematicLinearOrder4G2
];
dynamicForcingOrder4G2 = Expand[
  dynamicActualG2[[4]] - dynamicLinearOrder4G2
];
combinedForcingOrder4G2 = Expand[
  combinedActualG2[[4]] - combinedLinearOrder4G2
];

lowerZeroRulesG24 = Join[
  Table[
    etaG2[n] -> Function[{xx, yy, tt}, 0],
    {n, 1, 3}
  ],
  Table[
    phiG2[n] -> Function[{xx, yy, zz, tt}, 0],
    {n, 1, 3}
  ]
];

currentOrderLinearChecksG24 = {
  FullSimplify[
    (kinematicActualG2[[4]] /. lowerZeroRulesG24)
      - kinematicLinearOrder4G2
  ],
  FullSimplify[
    (dynamicActualG2[[4]] /. lowerZeroRulesG24)
      - dynamicLinearOrder4G2
  ],
  FullSimplify[
    (combinedActualG2[[4]] /. lowerZeroRulesG24)
      - combinedLinearOrder4G2
  ]
};

forcingSplitChecksG24 = {
  kinematicActualG2[[4]]
    - kinematicLinearOrder4G2 - kinematicForcingOrder4G2,
  dynamicActualG2[[4]]
    - dynamicLinearOrder4G2 - dynamicForcingOrder4G2,
  combinedActualG2[[4]]
    - combinedLinearOrder4G2 - combinedForcingOrder4G2
};

forcingBoundaryChecksG24 = {
  If[
    FreeQ[kinematicForcingOrder4G2, etaG2[4] | phiG2[4]],
    0,
    1
  ],
  If[
    FreeQ[dynamicForcingOrder4G2, etaG2[4] | phiG2[4]],
    0,
    1
  ],
  If[
    FreeQ[combinedForcingOrder4G2, etaG2[4] | phiG2[4]],
    0,
    1
  ]
};

gateG24[name_, checks_List] := Module[
  {pass = And @@ (TrueQ[PossibleZeroQ[FullSimplify[#]]] & /@ checks)},
  Print[name <> " = " <> If[pass, "PASS", "FAIL"]];
  If[!pass, Print[InputForm[FullSimplify /@ checks]]];
  pass
];

Print[
  "2D kinematic residual coefficient leaf counts through order four: ",
  InputForm[LeafCount /@ kinematicActualG2]
];
Print[
  "2D dynamic residual coefficient leaf counts through order four: ",
  InputForm[LeafCount /@ dynamicActualG2]
];
Print[
  "2D combined residual coefficient leaf counts through order four: ",
  InputForm[LeafCount /@ combinedActualG2]
];
Print[
  "2D order-four forcing leaf counts {K,D,C}: ",
  InputForm[LeafCount /@ {
    kinematicForcingOrder4G2,
    dynamicForcingOrder4G2,
    combinedForcingOrder4G2
  }]
];

gateResultsG24 = {
  gateG24[
    "2D orders one through three match the certified compiler",
    lowerRegressionChecksG24
  ],
  gateG24[
    "2D fourth-order current-unknown linear blocks",
    currentOrderLinearChecksG24
  ],
  gateG24[
    "2D fourth-order forcing split reproduces original residuals",
    forcingSplitChecksG24
  ],
  gateG24[
    "2D fourth-order forcing contains lower orders only",
    forcingBoundaryChecksG24
  ]
};

overallPassG24 = And @@ gateResultsG24;
Print[
  "OVERALL_ORDER1_TO_ORDER4_RESIDUAL_EXPANSION_2D = ",
  If[overallPassG24, "PASS", "FAIL"]
];
If[!overallPassG24, Exit[1]];
