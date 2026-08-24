(* ::Package:: *)
(* Exact finite-depth order-two Euler residual in two horizontal dimensions.

   This source derives the labelled positive-positive pair coefficients from
   the original two-dimensional free-surface conditions. No transfer kernel,
   VWA formula, sampled tuple, or validation oracle is imported.

   Dimensionless convention: g = h = 1,
     nu(q)^2 = q Tanh[q],
     q_i > 0, -1 < cosine12 < 1.
*)

expansionSourceFD2 = FileNameJoin[{
  DirectoryName[$InputFileName],
  "order1_to_order3_symbolic_expansion_2d.wl"
}];
If[!FileExistsQ[expansionSourceFD2],
  Print["FAIL: two-dimensional Euler expansion not found."];
  Exit[1]
];
Get[expansionSourceFD2];

ClearAll[
  amplitude1FD2, amplitude2FD2, q1FD2, q2FD2, cosineFD2,
  sineFD2, nuFD2, dnoFD2, qx1FD2, qy1FD2, qx2FD2, qy2FD2,
  dot12FD2, qOutputFD2, nu1FD2, nu2FD2, nuSumFD2,
  etaGroupedCoefficientFD2, phiGroupedCoefficientFD2,
  verticalJetFD2, eta1ModalFD2, eta2ModalFD2,
  phi1ModalFD2, phi2ModalFD2, zeroEtaFD2, zeroPhiFD2,
  modalRulesFD2, modalCoefficientFD2, rKFD2, rDFD2,
  kinematicEtaLinearFD2, kinematicPhiLinearFD2,
  kinematicForcingFD2, dynamicEtaLinearFD2,
  dynamicPhiLinearFD2, dynamicForcingFD2, determinantFD2,
  etaDirectionalCompactFD2, phiDirectionalCompactFD2,
  etaDirectionalFiniteFD2, phiDirectionalFiniteFD2,
  sourceKFD2, sourceDFD2, detuningFD2,
  rKCompactFD2, rDCompactFD2, compactChecksFD2,
  recoveredChecksFD2, swapResidualsFD2, diagonalResidualsFD2,
  deepWaterEtaFD2, deepWaterPhiFD2, projectRootFD2,
  artifactPathFD2, inputFormStringFD2, gateFD2,
  gateResultsFD2, overallPassFD2
];

$Assumptions = (
  q1FD2 > 0 && q2FD2 > 0 && -1 < cosineFD2 < 1
);

nuFD2[q_] := Sqrt[q Tanh[q]];
dnoFD2[q_] := q Tanh[q];
sineFD2 = Sqrt[1 - cosineFD2^2];
qx1FD2 = q1FD2;
qy1FD2 = 0;
qx2FD2 = q2FD2 cosineFD2;
qy2FD2 = q2FD2 sineFD2;
dot12FD2 = q1FD2 q2FD2 cosineFD2;
qOutputFD2 = Sqrt[
  q1FD2^2 + q2FD2^2 + 2 dot12FD2
];
nu1FD2 = nuFD2[q1FD2];
nu2FD2 = nuFD2[q2FD2];
nuSumFD2 = nu1FD2 + nu2FD2;

verticalJetFD2[q_, zz_] := (
  1 + dnoFD2[q] zz + q^2 zz^2/2
);

eta1ModalFD2[xx_, yy_, tt_] := (
  amplitude1FD2 Exp[I (qx1FD2 xx + qy1FD2 yy - nu1FD2 tt)]
  + amplitude2FD2 Exp[I (qx2FD2 xx + qy2FD2 yy - nu2FD2 tt)]
);
phi1ModalFD2[xx_, yy_, zz_, tt_] := (
  -I amplitude1FD2/nu1FD2 verticalJetFD2[q1FD2, zz]
    Exp[I (qx1FD2 xx + qy1FD2 yy - nu1FD2 tt)]
  -I amplitude2FD2/nu2FD2 verticalJetFD2[q2FD2, zz]
    Exp[I (qx2FD2 xx + qy2FD2 yy - nu2FD2 tt)]
);
eta2ModalFD2[xx_, yy_, tt_] := (
  etaGroupedCoefficientFD2 amplitude1FD2 amplitude2FD2
    Exp[I (
      (qx1FD2 + qx2FD2) xx
      + (qy1FD2 + qy2FD2) yy
      - nuSumFD2 tt
    )]
);
phi2ModalFD2[xx_, yy_, zz_, tt_] := (
  phiGroupedCoefficientFD2 amplitude1FD2 amplitude2FD2
    verticalJetFD2[qOutputFD2, zz]
    Exp[I (
      (qx1FD2 + qx2FD2) xx
      + (qy1FD2 + qy2FD2) yy
      - nuSumFD2 tt
    )]
);
zeroEtaFD2[xx_, yy_, tt_] := 0;
zeroPhiFD2[xx_, yy_, zz_, tt_] := 0;

modalRulesFD2 = {
  etaG2[1] -> eta1ModalFD2,
  etaG2[2] -> eta2ModalFD2,
  etaG2[3] -> zeroEtaFD2,
  phiG2[1] -> phi1ModalFD2,
  phiG2[2] -> phi2ModalFD2,
  phiG2[3] -> zeroPhiFD2
};

modalCoefficientFD2[expression_] := Coefficient[
  Coefficient[
    Expand[
      expression /. modalRulesFD2 /. {
        gravityG2 -> 1, xG2 -> 0, yG2 -> 0, tG2 -> 0
      }
    ],
    amplitude1FD2,
    1
  ],
  amplitude2FD2,
  1
];

Print["Extracting finite-depth directional order-two Euler residuals."];
rKFD2 = Together[modalCoefficientFD2[kinematicExpectedG2[[2]]]];
rDFD2 = Together[modalCoefficientFD2[dynamicExpectedG2[[2]]]];

kinematicEtaLinearFD2 = Coefficient[
  rKFD2, etaGroupedCoefficientFD2
];
kinematicPhiLinearFD2 = Coefficient[
  rKFD2, phiGroupedCoefficientFD2
];
kinematicForcingFD2 = rKFD2 /. {
  etaGroupedCoefficientFD2 -> 0,
  phiGroupedCoefficientFD2 -> 0
};
dynamicEtaLinearFD2 = Coefficient[
  rDFD2, etaGroupedCoefficientFD2
];
dynamicPhiLinearFD2 = Coefficient[
  rDFD2, phiGroupedCoefficientFD2
];
dynamicForcingFD2 = rDFD2 /. {
  etaGroupedCoefficientFD2 -> 0,
  phiGroupedCoefficientFD2 -> 0
};

determinantFD2 = Together[
  kinematicEtaLinearFD2 dynamicPhiLinearFD2
  - dynamicEtaLinearFD2 kinematicPhiLinearFD2
];
etaDirectionalFiniteFD2 = Together[
  (
    kinematicPhiLinearFD2 dynamicForcingFD2
    - dynamicPhiLinearFD2 kinematicForcingFD2
  )/determinantFD2
];
phiDirectionalFiniteFD2 = Together[
  (
    dynamicEtaLinearFD2 kinematicForcingFD2
    - kinematicEtaLinearFD2 dynamicForcingFD2
  )/determinantFD2
];

(* Compact forms are derived from the extracted affine coefficients and then
   certified against the original residuals. They are retained only to keep
   the order-three substitution graph manageable. *)
sourceKFD2 = (
  (q1FD2^2 + dot12FD2)/nu1FD2
  + (q2FD2^2 + dot12FD2)/nu2FD2
);
sourceDFD2 = (
  nu1FD2^2 + nu1FD2 nu2FD2 + nu2FD2^2
  - dot12FD2/(nu1FD2 nu2FD2)
);
detuningFD2 = dnoFD2[qOutputFD2] - nuSumFD2^2;
rKCompactFD2 = (
  -I nuSumFD2 etaGroupedCoefficientFD2
  - dnoFD2[qOutputFD2] phiGroupedCoefficientFD2
  + I sourceKFD2
);
rDCompactFD2 = (
  etaGroupedCoefficientFD2
  - I nuSumFD2 phiGroupedCoefficientFD2
  - sourceDFD2
);
compactChecksFD2 = Together /@ {
  rKFD2 - rKCompactFD2,
  rDFD2 - rDCompactFD2
};

etaDirectionalCompactFD2 = (
  (
    dnoFD2[qOutputFD2] sourceDFD2
    - nuSumFD2 sourceKFD2
  )/detuningFD2
);
phiDirectionalCompactFD2 = (
  -I (
    nuSumFD2 sourceDFD2 - sourceKFD2
  )/detuningFD2
);
etaDirectionalFiniteFD2 = Together[etaDirectionalCompactFD2];
phiDirectionalFiniteFD2 = Together[phiDirectionalCompactFD2];
recoveredChecksFD2 = Together /@ {
  rKFD2 /. {
    etaGroupedCoefficientFD2 -> etaDirectionalFiniteFD2,
    phiGroupedCoefficientFD2 -> phiDirectionalFiniteFD2
  },
  rDFD2 /. {
    etaGroupedCoefficientFD2 -> etaDirectionalFiniteFD2,
    phiGroupedCoefficientFD2 -> phiDirectionalFiniteFD2
  }
};
swapResidualsFD2 = FullSimplify[
  {
    etaDirectionalFiniteFD2
      - (
        etaDirectionalFiniteFD2 /. {
          q1FD2 -> q2FD2, q2FD2 -> q1FD2
        }
      ),
    phiDirectionalFiniteFD2
      - (
        phiDirectionalFiniteFD2 /. {
          q1FD2 -> q2FD2, q2FD2 -> q1FD2
        }
      )
  },
  Assumptions -> $Assumptions
];
diagonalResidualsFD2 = FullSimplify[
  {
    etaDirectionalFiniteFD2,
    phiDirectionalFiniteFD2
  } /. {
    q1FD2 -> q1FD2,
    q2FD2 -> q1FD2,
    cosineFD2 -> 1
  },
  Assumptions -> q1FD2 > 0
];
deepWaterEtaFD2 = PowerExpand[
  Together[etaDirectionalFiniteFD2 /. Tanh[_] -> 1]
];
deepWaterPhiFD2 = PowerExpand[
  Together[phiDirectionalFiniteFD2 /. Tanh[_] -> 1]
];

projectRootFD2 = DirectoryName[$InputFileName, 3];
artifactPathFD2 = FileNameJoin[{
  projectRootFD2,
  "artifacts",
  "order2_finite_depth_directional_euler_result.json"
}];
inputFormStringFD2[expression_] := ToString[InputForm[expression]];
Export[
  artifactPathFD2,
  <|
    "schema_version" -> 1,
    "domain" ->
      "finite-depth forward-directional order-two positive pure-sum",
    "candidate_input_fields" -> {"eta11"},
    "oracle_imported" -> False,
    "dispersion" -> "nu(q)^2=q Tanh[q]",
    "grouped_to_ordered_factor" -> 4,
    "eta_grouped" -> inputFormStringFD2[etaDirectionalFiniteFD2],
    "phi_grouped" -> inputFormStringFD2[phiDirectionalFiniteFD2],
    "detuning" -> inputFormStringFD2[detuningFD2],
    "deep_water_eta_grouped" -> inputFormStringFD2[deepWaterEtaFD2],
    "deep_water_phi_grouped" -> inputFormStringFD2[deepWaterPhiFD2]
  |>,
  "RawJSON"
];

gateFD2[name_, residual_] := Module[
  {pass = TrueQ[PossibleZeroQ[residual]]},
  Print[name <> " = " <> If[pass, "PASS", "FAIL"]];
  If[!pass, Print["  residual: ", InputForm[residual]]];
  pass
];
gateResultsFD2 = {
  gateFD2["raw Euler RK equals compact directional residual",
    compactChecksFD2[[1]]],
  gateFD2["raw Euler RD equals compact directional residual",
    compactChecksFD2[[2]]],
  gateFD2["recovered pair satisfies original RK",
    recoveredChecksFD2[[1]]],
  gateFD2["recovered pair satisfies original RD",
    recoveredChecksFD2[[2]]],
  gateFD2["pair coefficients are permutation symmetric",
    Total[swapResidualsFD2]]
};
overallPassFD2 = And @@ gateResultsFD2;
Print["=== Finite-depth directional order-two Euler gate ==="];
Print["eta grouped = ", InputForm[etaDirectionalFiniteFD2]];
Print["phi grouped = ", InputForm[phiDirectionalFiniteFD2]];
Print["artifact = ", artifactPathFD2];
Print[
  "OVERALL_ORDER2_FINITE_DEPTH_DIRECTIONAL_EULER = ",
  If[overallPassFD2, "PASS", "FAIL"]
];
If[!overallPassFD2, Exit[1]];
