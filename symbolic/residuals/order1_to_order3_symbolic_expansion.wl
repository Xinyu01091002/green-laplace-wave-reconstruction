(* Independent symbolic extraction of the epsilon^1--epsilon^3 residuals of
   the exact Euler free-surface boundary conditions.

   The unknowns phi[n][x,z,t] are BULK velocity-potential coefficients.  All
   boundary expressions are Taylor evaluated from z=eta to z=0 before their
   perturbation coefficients are extracted.
*)

ClearAll[epsilon, x, z, t, gravity, eta, phi];

etaSeries = Sum[epsilon^n eta[n][x, t], {n, 1, 3}];
phiSeries = Sum[epsilon^n phi[n][x, z, t], {n, 1, 3}];

truncate[expression_] := Normal[Series[expression, {epsilon, 0, 3}]];

surfaceTaylor[expression_] := truncate[
  Sum[
    etaSeries^m/m! (D[expression, {z, m}] /. z -> 0),
    {m, 0, 3}
  ]
];

phiXSurface = surfaceTaylor[D[phiSeries, x]];
phiZSurface = surfaceTaylor[D[phiSeries, z]];
phiTSurface = surfaceTaylor[D[phiSeries, t]];

kinematicSeries = truncate[
  D[etaSeries, t] + phiXSurface D[etaSeries, x] - phiZSurface
];

dynamicSeries = truncate[
  phiTSurface
    + 1/2 (phiXSurface^2 + phiZSurface^2)
    + gravity etaSeries
];

qBulk = D[phiSeries, x]^2 + D[phiSeries, z]^2;
combinedBulk = (
  D[phiSeries, {t, 2}] + gravity D[phiSeries, z]
    + D[qBulk, t]
    + 1/2 (
        D[phiSeries, x] D[qBulk, x]
        + D[phiSeries, z] D[qBulk, z]
      )
);
combinedSeries = surfaceTaylor[combinedBulk];

coefficient[series_, order_Integer] := Expand[Coefficient[series, epsilon, order]];
atZero[expression_] := expression /. z -> 0;

eta1 = eta[1][x, t];
eta2 = eta[2][x, t];
eta3 = eta[3][x, t];
phi1 = phi[1][x, z, t];
phi2 = phi[2][x, z, t];
phi3 = phi[3][x, z, t];

kinematicExpected = {
  atZero[D[eta1, t] - D[phi1, z]],
  atZero[
    D[eta2, t] - D[phi2, z]
      + D[phi1, x] D[eta1, x]
      - eta1 D[phi1, {z, 2}]
  ],
  atZero[
    D[eta3, t] - D[phi3, z]
      + D[phi1, x] D[eta2, x]
      + D[phi2, x] D[eta1, x]
      + eta1 D[phi1, x, z] D[eta1, x]
      - eta1 D[phi2, {z, 2}]
      - eta2 D[phi1, {z, 2}]
      - eta1^2/2 D[phi1, {z, 3}]
  ]
};

dynamicExpected = {
  atZero[D[phi1, t] + gravity eta1],
  atZero[
    D[phi2, t] + gravity eta2
      + eta1 D[phi1, t, z]
      + 1/2 (D[phi1, x]^2 + D[phi1, z]^2)
  ],
  atZero[
    D[phi3, t] + gravity eta3
      + eta1 D[phi2, t, z]
      + eta2 D[phi1, t, z]
      + eta1^2/2 D[phi1, t, {z, 2}]
      + D[phi1, x] D[phi2, x]
      + D[phi1, z] D[phi2, z]
      + eta1 (
          D[phi1, x] D[phi1, x, z]
          + D[phi1, z] D[phi1, {z, 2}]
        )
  ]
};

linearCombined[potential_] := (
  D[potential, {t, 2}] + gravity D[potential, z]
);
q2 = D[phi1, x]^2 + D[phi1, z]^2;
q3 = 2 (
  D[phi1, x] D[phi2, x]
  + D[phi1, z] D[phi2, z]
);

combinedExpected = {
  atZero[linearCombined[phi1]],
  atZero[
    linearCombined[phi2]
      + eta1 D[linearCombined[phi1], z]
      + D[q2, t]
  ],
  atZero[
    linearCombined[phi3]
      + eta1 D[linearCombined[phi2], z]
      + eta2 D[linearCombined[phi1], z]
      + eta1^2/2 D[linearCombined[phi1], {z, 2}]
      + D[q3, t]
      + eta1 D[q2, t, z]
      + 1/2 (
          D[phi1, x] D[q2, x]
          + D[phi1, z] D[q2, z]
        )
  ]
};

kinematicActual = Table[coefficient[kinematicSeries, n], {n, 1, 3}];
dynamicActual = Table[coefficient[dynamicSeries, n], {n, 1, 3}];
combinedActual = Table[coefficient[combinedSeries, n], {n, 1, 3}];

residualChecks = Join[
  Thread[kinematicActual - kinematicExpected],
  Thread[dynamicActual - dynamicExpected],
  Thread[combinedActual - combinedExpected]
];

Print["Kinematic residual coefficient leaf counts: ",
  InputForm[LeafCount /@ kinematicActual]];
Print["Dynamic residual coefficient leaf counts: ",
  InputForm[LeafCount /@ dynamicActual]];
Print["Combined residual coefficient leaf counts: ",
  InputForm[LeafCount /@ combinedActual]];

If[
  And @@ (TrueQ[PossibleZeroQ[FullSimplify[#]]] & /@ residualChecks),
  Print["PASS: exact order-1--order-3 boundary-residual expansions."],
  Print["FAIL: a hand-written residual coefficient does not match the direct expansion."];
  Print[InputForm[FullSimplify /@ residualChecks]];
  Exit[1]
];
