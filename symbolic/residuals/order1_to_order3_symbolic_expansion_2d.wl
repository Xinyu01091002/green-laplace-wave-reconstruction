(* ::Package:: *)
(* Exact two-horizontal-dimensional extraction of the epsilon^1--epsilon^3
   residuals of the Euler free-surface boundary conditions.

   The unknowns phiG2[n][x,y,z,t] are BULK velocity-potential coefficients.
   Every surface quantity is Taylor evaluated from z = eta to z = 0 before
   perturbation coefficients are collected.  This file contains no modal
   kernel, directional candidate, or validation oracle.
*)

ClearAll[
  epsilonG2, xG2, yG2, zG2, tG2, gravityG2, etaG2, phiG2,
  etaSeriesG2, phiSeriesG2, truncateG2, surfaceTaylorG2,
  phiXSurfaceG2, phiYSurfaceG2, phiZSurfaceG2, phiTSurfaceG2,
  kinematicSeriesG2, dynamicSeriesG2, qBulkG2, combinedBulkG2,
  combinedSeriesG2, coefficientG2, atZeroG2, horizontalDotG2,
  eta1G2, eta2G2, eta3G2, phi1G2, phi2G2, phi3G2,
  kinematicExpectedG2, dynamicExpectedG2, linearCombinedG2,
  q2G2, q3G2, combinedExpectedG2, kinematicActualG2,
  dynamicActualG2, combinedActualG2, residualChecksG2
];

etaSeriesG2 = Sum[
  epsilonG2^n etaG2[n][xG2, yG2, tG2],
  {n, 1, 3}
];
phiSeriesG2 = Sum[
  epsilonG2^n phiG2[n][xG2, yG2, zG2, tG2],
  {n, 1, 3}
];

truncateG2[expression_] := Normal[
  Series[expression, {epsilonG2, 0, 3}]
];

surfaceTaylorG2[expression_] := truncateG2[
  Sum[
    etaSeriesG2^m/m!
      (D[expression, {zG2, m}] /. zG2 -> 0),
    {m, 0, 3}
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

eta1G2 = etaG2[1][xG2, yG2, tG2];
eta2G2 = etaG2[2][xG2, yG2, tG2];
eta3G2 = etaG2[3][xG2, yG2, tG2];
phi1G2 = phiG2[1][xG2, yG2, zG2, tG2];
phi2G2 = phiG2[2][xG2, yG2, zG2, tG2];
phi3G2 = phiG2[3][xG2, yG2, zG2, tG2];

kinematicExpectedG2 = {
  atZeroG2[D[eta1G2, tG2] - D[phi1G2, zG2]],
  atZeroG2[
    D[eta2G2, tG2] - D[phi2G2, zG2]
    + horizontalDotG2[phi1G2, eta1G2]
    - eta1G2 D[phi1G2, {zG2, 2}]
  ],
  atZeroG2[
    D[eta3G2, tG2] - D[phi3G2, zG2]
    + horizontalDotG2[phi1G2, eta2G2]
    + horizontalDotG2[phi2G2, eta1G2]
    + eta1G2 horizontalDotG2[D[phi1G2, zG2], eta1G2]
    - eta1G2 D[phi2G2, {zG2, 2}]
    - eta2G2 D[phi1G2, {zG2, 2}]
    - eta1G2^2/2 D[phi1G2, {zG2, 3}]
  ]
};

dynamicExpectedG2 = {
  atZeroG2[D[phi1G2, tG2] + gravityG2 eta1G2],
  atZeroG2[
    D[phi2G2, tG2] + gravityG2 eta2G2
    + eta1G2 D[phi1G2, tG2, zG2]
    + 1/2 (
        D[phi1G2, xG2]^2
        + D[phi1G2, yG2]^2
        + D[phi1G2, zG2]^2
      )
  ],
  atZeroG2[
    D[phi3G2, tG2] + gravityG2 eta3G2
    + eta1G2 D[phi2G2, tG2, zG2]
    + eta2G2 D[phi1G2, tG2, zG2]
    + eta1G2^2/2 D[phi1G2, tG2, {zG2, 2}]
    + horizontalDotG2[phi1G2, phi2G2]
    + D[phi1G2, zG2] D[phi2G2, zG2]
    + eta1G2 (
        horizontalDotG2[phi1G2, D[phi1G2, zG2]]
        + D[phi1G2, zG2] D[phi1G2, {zG2, 2}]
      )
  ]
};

linearCombinedG2[potential_] := (
  D[potential, {tG2, 2}]
  + gravityG2 D[potential, zG2]
);
q2G2 = (
  D[phi1G2, xG2]^2
  + D[phi1G2, yG2]^2
  + D[phi1G2, zG2]^2
);
q3G2 = 2 (
  horizontalDotG2[phi1G2, phi2G2]
  + D[phi1G2, zG2] D[phi2G2, zG2]
);

combinedExpectedG2 = {
  atZeroG2[linearCombinedG2[phi1G2]],
  atZeroG2[
    linearCombinedG2[phi2G2]
    + eta1G2 D[linearCombinedG2[phi1G2], zG2]
    + D[q2G2, tG2]
  ],
  atZeroG2[
    linearCombinedG2[phi3G2]
    + eta1G2 D[linearCombinedG2[phi2G2], zG2]
    + eta2G2 D[linearCombinedG2[phi1G2], zG2]
    + eta1G2^2/2
      D[linearCombinedG2[phi1G2], {zG2, 2}]
    + D[q3G2, tG2]
    + eta1G2 D[q2G2, tG2, zG2]
    + 1/2 (
        D[phi1G2, xG2] D[q2G2, xG2]
        + D[phi1G2, yG2] D[q2G2, yG2]
        + D[phi1G2, zG2] D[q2G2, zG2]
      )
  ]
};

kinematicActualG2 = Table[
  coefficientG2[kinematicSeriesG2, n],
  {n, 1, 3}
];
dynamicActualG2 = Table[
  coefficientG2[dynamicSeriesG2, n],
  {n, 1, 3}
];
combinedActualG2 = Table[
  coefficientG2[combinedSeriesG2, n],
  {n, 1, 3}
];

residualChecksG2 = Join[
  Thread[kinematicActualG2 - kinematicExpectedG2],
  Thread[dynamicActualG2 - dynamicExpectedG2],
  Thread[combinedActualG2 - combinedExpectedG2]
];

Print[
  "2D kinematic residual coefficient leaf counts: ",
  InputForm[LeafCount /@ kinematicActualG2]
];
Print[
  "2D dynamic residual coefficient leaf counts: ",
  InputForm[LeafCount /@ dynamicActualG2]
];
Print[
  "2D combined residual coefficient leaf counts: ",
  InputForm[LeafCount /@ combinedActualG2]
];

If[
  And @@ (
    TrueQ[PossibleZeroQ[FullSimplify[#]]] & /@ residualChecksG2
  ),
  Print[
    "PASS: exact two-dimensional order-1--order-3 boundary-residual expansions."
  ],
  Print[
    "FAIL: a two-dimensional residual coefficient does not match the direct expansion."
  ];
  Print[InputForm[FullSimplify /@ residualChecksG2]];
  Exit[1]
];
