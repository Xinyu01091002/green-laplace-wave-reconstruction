(* ::Package:: *)
(* Finite-depth order-two mixed-sign Euler residual gate.

   This gate imports only the generic exact Euler free-surface expansion.
   It substitutes one positive parent and the internally generated conjugate
   branch of a second positive parent, then extracts A1 Conjugate[A2]
   directly. It contains no candidate approximation.

   Primary branch:
     h > 0, g > 0, k1 > k2 > 0,
     omega_j^2 = g k_j Tanh[k_j h],
     KDelta = k1-k2 > 0.

   The strict zero equations and the one-sided nonzero limit are derived and
   reported separately.
*)

baseScriptFD20 = FileNameJoin[{
  DirectoryName[$InputFileName],
  "order1_to_order3_symbolic_expansion.wl"
}];
If[!FileExistsQ[baseScriptFD20],
  Print["FAIL: generic Euler residual source not found: ", baseScriptFD20];
  Exit[1]
];
Get[baseScriptFD20];

ClearAll[
  amplitude1FD20, amplitude2BarFD20, amplitude1BarFD20, amplitude2FD20,
  k1FD20, k2FD20, hFD20, omega1FD20, omega2FD20,
  kDeltaFD20, omegaDeltaFD20, qFD20, verticalJetFD20,
  etaGroupedFD20, phiGroupedFD20, etaMirrorFD20, phiMirrorFD20,
  eta1PrimaryFD20, eta2PrimaryFD20, eta3ZeroFD20,
  phi1PrimaryFD20, phi2PrimaryFD20, phi3ZeroFD20,
  eta1MirrorFD20, eta2MirrorFD20, phi1MirrorFD20, phi2MirrorFD20,
  primaryRulesFD20, mirrorRulesFD20, labelledCoefficientFD20,
  dispersionRelationsFD20, dispersionReduceFD20,
  rKRawFD20, rDRawFD20, rCRawFD20, rKFD20, rDFD20, rCFD20,
  forcingKFD20, forcingDFD20, forcingCFD20,
  q1FD20, q2FD20, qDeltaSurfaceFD20, detuningFD20,
  linearMatrixFD20, linearChecksFD20, combinedIdentityFD20,
  forcingIdentityFD20, etaOnlyResidualFD20, etaOnlyForcingFD20,
  phiRecoveredFD20, etaRecoveredFD20, phiOrderedFD20, etaOrderedFD20,
  rKRecoveredFD20, rDRecoveredFD20, rCRecoveredFD20,
  etaRecoveredResidualFD20,
  rKMirrorRawFD20, rDMirrorRawFD20, rCMirrorRawFD20,
  rKMirrorFD20, rDMirrorFD20, rCMirrorFD20,
  mirrorForcingCFD20, phiMirrorRecoveredFD20, etaMirrorRecoveredFD20,
  rKMirrorRecoveredFD20, rDMirrorRecoveredFD20,
  rCMirrorRecoveredFD20, mirrorEtaResidualFD20, mirrorPhiResidualFD20,
  deepWaterRulesFD20, deepEtaOrderedFD20, deepPhiOrderedFD20,
  q1dFD20, q2dFD20, nu1dFD20, nu2dFD20, qDeltadFD20,
  sigmaDeltadFD20, outputDnodFD20, detuningdFD20, ratioFD20,
  forcingKdFD20, forcingDdFD20, forcingCdFD20,
  etaGroupeddFD20, phiGroupeddFD20, etaOrdereddFD20,
  epsilonFD20, carrierQFD20, oneSidedRatioFD20,
  oneSidedDetuningScaleFD20, oneSidedEtaOrderedFD20,
  strictWaveNumberFD20, strictFrequencyFD20, strictRulesFD20,
  strictZeroMatrixFD20, strictZeroRankFD20,
  strictZeroKinematicFD20, strictZeroDynamicFD20,
  strictZeroCombinedFD20, strictZeroForcingKFD20,
  strictZeroForcingDFD20, strictZeroForcingCFD20,
  exactTextFD20, finiteExactQFD20, gateFD20, booleanGateFD20,
  gateNamesFD20, gateResultsFD20, gateAssociationFD20,
  overallPassFD20, projectRootFD20, artifactDirectoryFD20,
  artifactPathFD20
];

$Assumptions = (
  gravity > 0 && hFD20 > 0
  && k1FD20 > k2FD20 && k2FD20 > 0
  && omega1FD20 > omega2FD20 && omega2FD20 > 0
);

qFD20[waveNumber_] := waveNumber Tanh[waveNumber hFD20];
verticalJetFD20[waveNumber_, zz_] := (
  1 + qFD20[waveNumber] zz + waveNumber^2 zz^2/2
);
kDeltaFD20 = k1FD20-k2FD20;
omegaDeltaFD20 = omega1FD20-omega2FD20;
q1FD20 = qFD20[k1FD20];
q2FD20 = qFD20[k2FD20];
qDeltaSurfaceFD20 = qFD20[kDeltaFD20];
detuningFD20 = gravity qDeltaSurfaceFD20-omegaDeltaFD20^2;

eta1PrimaryFD20[xx_, tt_] := (
  amplitude1FD20 Exp[I (k1FD20 xx-omega1FD20 tt)]
  +amplitude2BarFD20 Exp[-I (k2FD20 xx-omega2FD20 tt)]
);
phi1PrimaryFD20[xx_, zz_, tt_] := (
  -I gravity amplitude1FD20/omega1FD20
    verticalJetFD20[k1FD20, zz]
    Exp[I (k1FD20 xx-omega1FD20 tt)]
  +I gravity amplitude2BarFD20/omega2FD20
    verticalJetFD20[k2FD20, zz]
    Exp[-I (k2FD20 xx-omega2FD20 tt)]
);
eta2PrimaryFD20[xx_, tt_] := (
  etaGroupedFD20 amplitude1FD20 amplitude2BarFD20
    Exp[I (kDeltaFD20 xx-omegaDeltaFD20 tt)]
);
phi2PrimaryFD20[xx_, zz_, tt_] := (
  phiGroupedFD20 amplitude1FD20 amplitude2BarFD20
    verticalJetFD20[kDeltaFD20, zz]
    Exp[I (kDeltaFD20 xx-omegaDeltaFD20 tt)]
);
eta3ZeroFD20[xx_, tt_] := 0;
phi3ZeroFD20[xx_, zz_, tt_] := 0;

primaryRulesFD20 = {
  eta[1] -> eta1PrimaryFD20,
  eta[2] -> eta2PrimaryFD20,
  eta[3] -> eta3ZeroFD20,
  phi[1] -> phi1PrimaryFD20,
  phi[2] -> phi2PrimaryFD20,
  phi[3] -> phi3ZeroFD20
};

labelledCoefficientFD20[
  expression_, rules_, firstAmplitude_, secondAmplitude_
] := Coefficient[
  Coefficient[
    Expand[expression /. rules /. {x -> 0, t -> 0}],
    firstAmplitude,
    1
  ],
  secondAmplitude,
  1
];

rKRawFD20 = Together[labelledCoefficientFD20[
  kinematicActual[[2]], primaryRulesFD20,
  amplitude1FD20, amplitude2BarFD20
]];
rDRawFD20 = Together[labelledCoefficientFD20[
  dynamicActual[[2]], primaryRulesFD20,
  amplitude1FD20, amplitude2BarFD20
]];
rCRawFD20 = Together[labelledCoefficientFD20[
  combinedActual[[2]], primaryRulesFD20,
  amplitude1FD20, amplitude2BarFD20
]];
Print["Finite-depth eta20 primary Euler coefficients extracted."];

dispersionRelationsFD20 = {
  omega1FD20^2-gravity q1FD20,
  omega2FD20^2-gravity q2FD20
};
dispersionReduceFD20[expression_] := Module[
  {rational, numerator, denominator, remainder},
  rational = Together[expression];
  numerator = Expand[Numerator[rational]];
  denominator = Denominator[rational];
  remainder = Last[PolynomialReduce[
    numerator,
    dispersionRelationsFD20,
    {omega1FD20, omega2FD20}
  ]];
  Cancel[Together[remainder/denominator]]
];

rKFD20 = dispersionReduceFD20[rKRawFD20];
rDFD20 = dispersionReduceFD20[rDRawFD20];
rCFD20 = dispersionReduceFD20[rCRawFD20];
forcingKFD20 = dispersionReduceFD20[
  rKFD20 /. {etaGroupedFD20 -> 0, phiGroupedFD20 -> 0}
];
forcingDFD20 = dispersionReduceFD20[
  rDFD20 /. {etaGroupedFD20 -> 0, phiGroupedFD20 -> 0}
];
forcingCFD20 = dispersionReduceFD20[
  rCFD20 /. {etaGroupedFD20 -> 0, phiGroupedFD20 -> 0}
];

linearMatrixFD20 = {
  {-I omegaDeltaFD20, -qDeltaSurfaceFD20},
  {gravity, -I omegaDeltaFD20}
};
linearChecksFD20 = dispersionReduceFD20 /@ {
  Coefficient[rKFD20, etaGroupedFD20]+I omegaDeltaFD20,
  Coefficient[rKFD20, phiGroupedFD20]+qDeltaSurfaceFD20,
  Coefficient[rDFD20, etaGroupedFD20]-gravity,
  Coefficient[rDFD20, phiGroupedFD20]+I omegaDeltaFD20,
  Coefficient[rCFD20, etaGroupedFD20],
  Coefficient[rCFD20, phiGroupedFD20]-detuningFD20
};
combinedIdentityFD20 = dispersionReduceFD20[
  rCFD20+I omegaDeltaFD20 rDFD20+gravity rKFD20
];
forcingIdentityFD20 = dispersionReduceFD20[
  forcingCFD20+I omegaDeltaFD20 forcingDFD20
    +gravity forcingKFD20
];
Print["Finite-depth eta20 current-order linear system derived."];

(* QDelta RD-I OmegaDelta RK eliminates Phi20 without dividing by
   OmegaDelta. This is the total physical elevation residual used by the
   candidate route. *)
etaOnlyResidualFD20 = dispersionReduceFD20[
  qDeltaSurfaceFD20 rDFD20-I omegaDeltaFD20 rKFD20
];
etaOnlyForcingFD20 = dispersionReduceFD20[
  etaOnlyResidualFD20 /. etaGroupedFD20 -> 0
];

phiRecoveredFD20 = Cancel[-forcingCFD20/detuningFD20];
etaRecoveredFD20 = Cancel[-etaOnlyForcingFD20/detuningFD20];
phiOrderedFD20 = Cancel[phiRecoveredFD20/2];
etaOrderedFD20 = Cancel[etaRecoveredFD20/2];

rKRecoveredFD20 = dispersionReduceFD20[
  rKFD20 /. {
    etaGroupedFD20 -> etaRecoveredFD20,
    phiGroupedFD20 -> phiRecoveredFD20
  }
];
rDRecoveredFD20 = dispersionReduceFD20[
  rDFD20 /. {
    etaGroupedFD20 -> etaRecoveredFD20,
    phiGroupedFD20 -> phiRecoveredFD20
  }
];
rCRecoveredFD20 = dispersionReduceFD20[
  rCFD20 /. phiGroupedFD20 -> phiRecoveredFD20
];
etaRecoveredResidualFD20 = dispersionReduceFD20[
  etaOnlyResidualFD20 /. etaGroupedFD20 -> etaRecoveredFD20
];
Print["Finite-depth eta20 numerator/determinant pair recovered."];

(* Independently extract the negative-output Hermitian companion. *)
eta1MirrorFD20[xx_, tt_] := (
  amplitude1BarFD20 Exp[-I (k1FD20 xx-omega1FD20 tt)]
  +amplitude2FD20 Exp[I (k2FD20 xx-omega2FD20 tt)]
);
phi1MirrorFD20[xx_, zz_, tt_] := (
  I gravity amplitude1BarFD20/omega1FD20
    verticalJetFD20[k1FD20, zz]
    Exp[-I (k1FD20 xx-omega1FD20 tt)]
  -I gravity amplitude2FD20/omega2FD20
    verticalJetFD20[k2FD20, zz]
    Exp[I (k2FD20 xx-omega2FD20 tt)]
);
eta2MirrorFD20[xx_, tt_] := (
  etaMirrorFD20 amplitude1BarFD20 amplitude2FD20
    Exp[-I (kDeltaFD20 xx-omegaDeltaFD20 tt)]
);
phi2MirrorFD20[xx_, zz_, tt_] := (
  phiMirrorFD20 amplitude1BarFD20 amplitude2FD20
    verticalJetFD20[kDeltaFD20, zz]
    Exp[-I (kDeltaFD20 xx-omegaDeltaFD20 tt)]
);
mirrorRulesFD20 = {
  eta[1] -> eta1MirrorFD20,
  eta[2] -> eta2MirrorFD20,
  eta[3] -> eta3ZeroFD20,
  phi[1] -> phi1MirrorFD20,
  phi[2] -> phi2MirrorFD20,
  phi[3] -> phi3ZeroFD20
};
rKMirrorRawFD20 = Together[labelledCoefficientFD20[
  kinematicActual[[2]], mirrorRulesFD20,
  amplitude1BarFD20, amplitude2FD20
]];
rDMirrorRawFD20 = Together[labelledCoefficientFD20[
  dynamicActual[[2]], mirrorRulesFD20,
  amplitude1BarFD20, amplitude2FD20
]];
rCMirrorRawFD20 = Together[labelledCoefficientFD20[
  combinedActual[[2]], mirrorRulesFD20,
  amplitude1BarFD20, amplitude2FD20
]];
rKMirrorFD20 = dispersionReduceFD20[rKMirrorRawFD20];
rDMirrorFD20 = dispersionReduceFD20[rDMirrorRawFD20];
rCMirrorFD20 = dispersionReduceFD20[rCMirrorRawFD20];
mirrorForcingCFD20 = dispersionReduceFD20[
  rCMirrorFD20 /. {etaMirrorFD20 -> 0, phiMirrorFD20 -> 0}
];
phiMirrorRecoveredFD20 = Cancel[-mirrorForcingCFD20/detuningFD20];
etaMirrorRecoveredFD20 = Cancel[
  etaMirrorFD20 /. First[Solve[
    (rDMirrorFD20 /. phiMirrorFD20 -> phiMirrorRecoveredFD20) == 0,
    etaMirrorFD20
  ]]
];
rKMirrorRecoveredFD20 = dispersionReduceFD20[
  rKMirrorFD20 /. {
    etaMirrorFD20 -> etaMirrorRecoveredFD20,
    phiMirrorFD20 -> phiMirrorRecoveredFD20
  }
];
rDMirrorRecoveredFD20 = dispersionReduceFD20[
  rDMirrorFD20 /. {
    etaMirrorFD20 -> etaMirrorRecoveredFD20,
    phiMirrorFD20 -> phiMirrorRecoveredFD20
  }
];
rCMirrorRecoveredFD20 = dispersionReduceFD20[
  rCMirrorFD20 /. phiMirrorFD20 -> phiMirrorRecoveredFD20
];
mirrorEtaResidualFD20 = dispersionReduceFD20[
  etaMirrorRecoveredFD20-etaRecoveredFD20
];
mirrorPhiResidualFD20 = dispersionReduceFD20[
  phiMirrorRecoveredFD20+phiRecoveredFD20
];
Print["Finite-depth eta20 Hermitian mirror extracted independently."];

deepWaterRulesFD20 = {
  Tanh[_] -> 1,
  omega1FD20 -> Sqrt[gravity k1FD20],
  omega2FD20 -> Sqrt[gravity k2FD20]
};
deepEtaOrderedFD20 = PowerExpand[Cancel[Together[
  etaOrderedFD20 /. deepWaterRulesFD20
]]];
deepPhiOrderedFD20 = PowerExpand[Cancel[Together[
  phiOrderedFD20 /. deepWaterRulesFD20
]]];

(* Dimensionless one-sided nonzero limit. *)
nu1dFD20 = Sqrt[q1dFD20 Tanh[q1dFD20]];
nu2dFD20 = Sqrt[q2dFD20 Tanh[q2dFD20]];
qDeltadFD20 = q1dFD20-q2dFD20;
sigmaDeltadFD20 = nu1dFD20-nu2dFD20;
outputDnodFD20 = qDeltadFD20 Tanh[qDeltadFD20];
detuningdFD20 = outputDnodFD20-sigmaDeltadFD20^2;
ratioFD20 = sigmaDeltadFD20^2/outputDnodFD20;
forcingKdFD20 = forcingKFD20 /. {
  gravity -> 1, hFD20 -> 1,
  k1FD20 -> q1dFD20, k2FD20 -> q2dFD20,
  omega1FD20 -> nu1dFD20, omega2FD20 -> nu2dFD20
};
forcingDdFD20 = forcingDFD20 /. {
  gravity -> 1, hFD20 -> 1,
  k1FD20 -> q1dFD20, k2FD20 -> q2dFD20,
  omega1FD20 -> nu1dFD20, omega2FD20 -> nu2dFD20
};
forcingCdFD20 = forcingCFD20 /. {
  gravity -> 1, hFD20 -> 1,
  k1FD20 -> q1dFD20, k2FD20 -> q2dFD20,
  omega1FD20 -> nu1dFD20, omega2FD20 -> nu2dFD20
};
phiGroupeddFD20 = Cancel[-forcingCdFD20/detuningdFD20];
etaGroupeddFD20 = Cancel[
  -(
    outputDnodFD20 forcingDdFD20
    -I sigmaDeltadFD20 forcingKdFD20
  )/detuningdFD20
];
etaOrdereddFD20 = Cancel[etaGroupeddFD20/2];

oneSidedRatioFD20 = FullSimplify[
  Limit[
    ratioFD20 /. {
      q1dFD20 -> carrierQFD20+epsilonFD20,
      q2dFD20 -> carrierQFD20
    },
    epsilonFD20 -> 0,
    Direction -> "FromAbove"
  ],
  Assumptions -> carrierQFD20 > 0
];
oneSidedDetuningScaleFD20 = FullSimplify[
  Limit[
    (
      detuningdFD20/epsilonFD20^2 /. {
        q1dFD20 -> carrierQFD20+epsilonFD20,
        q2dFD20 -> carrierQFD20
      }
    ),
    epsilonFD20 -> 0,
    Direction -> "FromAbove"
  ],
  Assumptions -> carrierQFD20 > 0
];
oneSidedEtaOrderedFD20 = TimeConstrained[
  FullSimplify[
    Limit[
      etaOrdereddFD20 /. {
        q1dFD20 -> carrierQFD20+epsilonFD20,
        q2dFD20 -> carrierQFD20
      },
      epsilonFD20 -> 0,
      Direction -> "FromAbove"
    ],
    Assumptions -> carrierQFD20 > 0
  ],
  90,
  $Failed
];

(* The strict zero is substituted into the original extracted equations,
   not obtained by taking the nonzero solution limit. *)
strictRulesFD20 = {
  k1FD20 -> strictWaveNumberFD20,
  k2FD20 -> strictWaveNumberFD20,
  omega1FD20 -> strictFrequencyFD20,
  omega2FD20 -> strictFrequencyFD20
};
strictZeroMatrixFD20 = FullSimplify[
  linearMatrixFD20 /. strictRulesFD20,
  Assumptions -> (
    gravity > 0 && hFD20 > 0 && strictWaveNumberFD20 > 0
    && strictFrequencyFD20^2
      == gravity qFD20[strictWaveNumberFD20]
  )
];
strictZeroRankFD20 = MatrixRank[strictZeroMatrixFD20];
strictZeroKinematicFD20 = FullSimplify[
  rKFD20 /. strictRulesFD20,
  Assumptions -> (
    gravity > 0 && hFD20 > 0 && strictWaveNumberFD20 > 0
    && strictFrequencyFD20^2
      == gravity qFD20[strictWaveNumberFD20]
  )
];
strictZeroDynamicFD20 = FullSimplify[
  rDFD20 /. strictRulesFD20,
  Assumptions -> (
    gravity > 0 && hFD20 > 0 && strictWaveNumberFD20 > 0
    && strictFrequencyFD20^2
      == gravity qFD20[strictWaveNumberFD20]
  )
];
strictZeroCombinedFD20 = FullSimplify[
  rCFD20 /. strictRulesFD20,
  Assumptions -> (
    gravity > 0 && hFD20 > 0 && strictWaveNumberFD20 > 0
    && strictFrequencyFD20^2
      == gravity qFD20[strictWaveNumberFD20]
  )
];
strictZeroForcingKFD20 = strictZeroKinematicFD20 /. {
  etaGroupedFD20 -> 0, phiGroupedFD20 -> 0
};
strictZeroForcingDFD20 = strictZeroDynamicFD20 /. {
  etaGroupedFD20 -> 0, phiGroupedFD20 -> 0
};
strictZeroForcingCFD20 = strictZeroCombinedFD20 /. {
  etaGroupedFD20 -> 0, phiGroupedFD20 -> 0
};
Print["Finite-depth eta20 strict zero and one-sided limit analyzed."];

exactTextFD20[expression_] := ToString[InputForm[expression]];
finiteExactQFD20[expression_] := (
  expression =!= $Failed
  && FreeQ[expression, Indeterminate]
  && FreeQ[expression, ComplexInfinity]
  && FreeQ[expression, _DirectedInfinity]
);
gateFD20[name_, residual_] := Module[
  {pass = TrueQ[PossibleZeroQ[residual]]},
  Print[name <> " = " <> If[pass, "PASS", "FAIL"]];
  If[!pass, Print["  residual: ", InputForm[residual]]];
  pass
];
booleanGateFD20[name_, condition_] := Module[
  {pass = TrueQ[condition]},
  Print[name <> " = " <> If[pass, "PASS", "FAIL"]];
  pass
];

gateNamesFD20 = {
  "kinematic eta linear block",
  "kinematic phi linear block",
  "dynamic eta linear block",
  "dynamic phi linear block",
  "combined eta cancels",
  "combined determinant is finite-depth difference detuning",
  "combined residual identity",
  "combined forcing identity",
  "total elevation residual after exact recovery",
  "original kinematic residual after exact recovery",
  "original dynamic residual after exact recovery",
  "combined residual after exact recovery",
  "mirror kinematic residual after exact recovery",
  "mirror dynamic residual after exact recovery",
  "mirror combined residual after exact recovery",
  "mirror elevation is Hermitian",
  "mirror flat potential is Hermitian",
  "deep-water ordered elevation regression",
  "deep-water ordered flat-potential regression",
  "strict zero current-order rank is one",
  "one-sided nonzero ratio is exact and finite",
  "one-sided nonzero elevation is exact and finite"
};
gateResultsFD20 = {
  gateFD20[gateNamesFD20[[1]], linearChecksFD20[[1]]],
  gateFD20[gateNamesFD20[[2]], linearChecksFD20[[2]]],
  gateFD20[gateNamesFD20[[3]], linearChecksFD20[[3]]],
  gateFD20[gateNamesFD20[[4]], linearChecksFD20[[4]]],
  gateFD20[gateNamesFD20[[5]], linearChecksFD20[[5]]],
  gateFD20[gateNamesFD20[[6]], linearChecksFD20[[6]]],
  gateFD20[gateNamesFD20[[7]], combinedIdentityFD20],
  gateFD20[gateNamesFD20[[8]], forcingIdentityFD20],
  gateFD20[gateNamesFD20[[9]], etaRecoveredResidualFD20],
  gateFD20[gateNamesFD20[[10]], rKRecoveredFD20],
  gateFD20[gateNamesFD20[[11]], rDRecoveredFD20],
  gateFD20[gateNamesFD20[[12]], rCRecoveredFD20],
  gateFD20[gateNamesFD20[[13]], rKMirrorRecoveredFD20],
  gateFD20[gateNamesFD20[[14]], rDMirrorRecoveredFD20],
  gateFD20[gateNamesFD20[[15]], rCMirrorRecoveredFD20],
  gateFD20[gateNamesFD20[[16]], mirrorEtaResidualFD20],
  gateFD20[gateNamesFD20[[17]], mirrorPhiResidualFD20],
  gateFD20[
    gateNamesFD20[[18]],
    PowerExpand[Cancel[
      Together[deepEtaOrderedFD20+(k1FD20-k2FD20)/2]
    ]]
  ],
  gateFD20[
    gateNamesFD20[[19]],
    PowerExpand[Cancel[
      Together[deepPhiOrderedFD20-I Sqrt[gravity k1FD20]]
    ]]
  ],
  booleanGateFD20[gateNamesFD20[[20]], strictZeroRankFD20 == 1],
  booleanGateFD20[
    gateNamesFD20[[21]],
    finiteExactQFD20[oneSidedRatioFD20]
    && finiteExactQFD20[oneSidedDetuningScaleFD20]
  ],
  booleanGateFD20[
    gateNamesFD20[[22]],
    finiteExactQFD20[oneSidedEtaOrderedFD20]
  ]
};
gateAssociationFD20 = AssociationThread[gateNamesFD20, gateResultsFD20];
overallPassFD20 = And @@ gateResultsFD20;

projectRootFD20 = DirectoryName[$InputFileName, 3];
artifactDirectoryFD20 = FileNameJoin[{projectRootFD20, "artifacts"}];
If[!DirectoryQ[artifactDirectoryFD20],
  CreateDirectory[
    artifactDirectoryFD20,
    CreateIntermediateDirectories -> True
  ]
];
artifactPathFD20 = FileNameJoin[{
  artifactDirectoryFD20,
  "order2_finite_depth_eta20_exact.json"
}];
Export[
  artifactPathFD20,
  <|
    "schema_version" -> 1,
    "status" -> "exact-residual-only-no-candidate-frozen",
    "source" -> "generic exact Euler free-surface expansion",
    "candidate_stage" -> False,
    "candidate_input_fields" -> {"eta11"},
    "candidate_selection_or_fitting" -> False,
    "root_result_used" -> False,
    "oracle_or_mf12_used" -> False,
    "grouped_to_ordered_factor" -> 2,
    "symbolic_domain" ->
      "h>0, g>0, k1>k2>0; no component-depth cutoff",
    "field_applicability" ->
      "peak k_p h >= 3/10; k h < 3/10 tails retained and reported",
    "residuals" -> <|
      "kinematic_grouped" -> exactTextFD20[rKFD20],
      "dynamic_grouped" -> exactTextFD20[rDFD20],
      "combined_grouped" -> exactTextFD20[rCFD20],
      "elevation_total_grouped" -> exactTextFD20[etaOnlyResidualFD20],
      "forcing_kinematic" -> exactTextFD20[forcingKFD20],
      "forcing_dynamic" -> exactTextFD20[forcingDFD20],
      "forcing_combined" -> exactTextFD20[forcingCFD20]
    |>,
    "linear_system" -> <|
      "matrix" -> exactTextFD20[linearMatrixFD20],
      "determinant" -> exactTextFD20[detuningFD20],
      "identity" -> "R_C = -(I OmegaDelta R_D + g R_K)"
    |>,
    "exact_numerator_determinant" -> <|
      "Phi20_grouped" -> exactTextFD20[phiRecoveredFD20],
      "E20_grouped" -> exactTextFD20[etaRecoveredFD20],
      "Phi20_ordered" -> exactTextFD20[phiOrderedFD20],
      "E20_ordered" -> exactTextFD20[etaOrderedFD20]
    |>,
    "one_sided_nonzero" -> <|
      "r_h" -> exactTextFD20[oneSidedRatioFD20],
      "detuning_over_qDelta_squared" ->
        exactTextFD20[oneSidedDetuningScaleFD20],
      "E20_ordered" -> exactTextFD20[oneSidedEtaOrderedFD20]
    |>,
    "strict_zero" -> <|
      "linear_matrix" -> exactTextFD20[strictZeroMatrixFD20],
      "rank" -> strictZeroRankFD20,
      "kinematic_equation" -> exactTextFD20[strictZeroKinematicFD20],
      "dynamic_equation" -> exactTextFD20[strictZeroDynamicFD20],
      "combined_equation" -> exactTextFD20[strictZeroCombinedFD20],
      "kinematic_forcing" -> exactTextFD20[strictZeroForcingKFD20],
      "dynamic_forcing" -> exactTextFD20[strictZeroForcingDFD20],
      "combined_forcing" -> exactTextFD20[strictZeroForcingCFD20],
      "candidate_decision" ->
        "excluded pending mean-volume-return-flow and time-dependent potential-gauge convention",
      "physical_mean_claim" -> False
    |>,
    "gates" -> gateAssociationFD20,
    "overall_pass" -> overallPassFD20
  |>,
  "RawJSON"
];

Print["=== Finite-depth order-two eta20 Euler gate ==="];
Print["total elevation residual = ", InputForm[etaOnlyResidualFD20]];
Print["determinant = ", InputForm[detuningFD20]];
Print["strict zero matrix = ", InputForm[strictZeroMatrixFD20]];
Print["strict zero rank = ", strictZeroRankFD20];
Print["strict zero RK = ", InputForm[strictZeroKinematicFD20]];
Print["strict zero RD = ", InputForm[strictZeroDynamicFD20]];
Print["strict zero RC = ", InputForm[strictZeroCombinedFD20]];
Print["one-sided r_h = ", InputForm[oneSidedRatioFD20]];
Print[
  "one-sided detuning/qDelta^2 = ",
  InputForm[oneSidedDetuningScaleFD20]
];
Print["one-sided E20 ordered = ", InputForm[oneSidedEtaOrderedFD20]];
Print["deep-water E20 ordered = ", InputForm[deepEtaOrderedFD20]];
Print["deep-water Phi20 ordered = ", InputForm[deepPhiOrderedFD20]];
Print["result artifact = ", artifactPathFD20];
Print[
  "OVERALL_ORDER2_FINITE_DEPTH_ETA20_EULER = ",
  If[overallPassFD20, "PASS", "FAIL"]
];
If[!overallPassFD20, Exit[1]];
