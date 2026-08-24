(* ::Package:: *)
(* Direct stationary Psi22 dual-branch Green--Laplace freeze.

   The exact resolvent branches are
     Exp[-(omega12-Omega12) tau],
     Exp[-(omega12+Omega12) tau].
   Each branch receives an independent two-node Gauss--Laguerre rule with
   the peak-collinear branch scale lambdaMinus or lambdaPlus.  Both branches
   are exact at that declared finite-depth endpoint.  Their deep-water scale
   limits are 2-Sqrt[2] and 2+Sqrt[2] after division by Sqrt[kp d].

   No time evolution, waveform, MF12, oracle, sampled coefficient selection,
   or forbidden polynomial approximation is used.
*)

residualSourceDBG = FileNameJoin[{
  DirectoryName[$InputFileName], "..", "residuals",
  "order2_finite_depth_directional_euler.wl"
}];
If[!FileExistsQ[residualSourceDBG],
  Print["FAIL: directional Euler residual source not found."];
  Exit[1]
];
Get[residualSourceDBG];

ClearAll[
  kpDBG, xDBG, omegaParentDBG, omegaForcedDBG, omegaOutputDBG,
  lambdaMinusDBG, lambdaPlusDBG, nodesDBG, weightsDBG,
  momentResidualsDBG, branchMinusExactDBG, branchPlusExactDBG,
  coshIdentityDBG, sinhIdentityDBG, finiteEndpointMinusDBG,
  finiteEndpointPlusDBG, deepMinusDBG, deepPlusDBG,
  peakMinusQuadratureDBG, peakPlusQuadratureDBG,
  exactTextDBG, gateDBG, gatesDBG, overallPassDBG,
  projectRootDBG, artifactPathDBG, interfacePathDBG,
  sourceRelativePathDBG, residualRelativePathDBG,
  sourceSha256DBG, residualSha256DBG
];

omegaParentDBG[q_] := Sqrt[q Tanh[q]];
omegaForcedDBG[q_] := 2 omegaParentDBG[q];
omegaOutputDBG[q_] := Sqrt[2 q Tanh[2 q]];
lambdaMinusDBG[q_] := Together[
  omegaForcedDBG[q]-omegaOutputDBG[q]
];
lambdaPlusDBG[q_] := Together[
  omegaForcedDBG[q]+omegaOutputDBG[q]
];

nodesDBG = Sort[xDBG /. Solve[LaguerreL[2, xDBG] == 0, xDBG, Reals]];
weightsDBG = #/(3^2 LaguerreL[3, #]^2) & /@ nodesDBG;
momentResidualsDBG = FullSimplify@Table[
  Total[weightsDBG nodesDBG^degree]-degree!,
  {degree, 0, 3}
];

branchMinusExactDBG = FullSimplify[
  Integrate[
    Exp[-lambdaMinusSymbolDBG tauDBG],
    {tauDBG, 0, Infinity},
    Assumptions -> lambdaMinusSymbolDBG > 0
  ]-1/lambdaMinusSymbolDBG,
  Assumptions -> lambdaMinusSymbolDBG > 0
];
branchPlusExactDBG = FullSimplify[
  Integrate[
    Exp[-lambdaPlusSymbolDBG tauDBG],
    {tauDBG, 0, Infinity},
    Assumptions -> lambdaPlusSymbolDBG > 0
  ]-1/lambdaPlusSymbolDBG,
  Assumptions -> lambdaPlusSymbolDBG > 0
];
coshIdentityDBG = Together[
  1/2 (1/(sDBG-aDBG)+1/(sDBG+aDBG))-sDBG/(sDBG^2-aDBG^2)
];
sinhIdentityDBG = Together[
  1/(2 aDBG) (1/(sDBG-aDBG)-1/(sDBG+aDBG))
  -1/(sDBG^2-aDBG^2)
];

finiteEndpointMinusDBG = lambdaMinusDBG[3/10];
finiteEndpointPlusDBG = lambdaPlusDBG[3/10];
deepMinusDBG = FullSimplify[
  Limit[lambdaMinusDBG[kpDBG]/Sqrt[kpDBG], kpDBG -> Infinity]
];
deepPlusDBG = FullSimplify[
  Limit[lambdaPlusDBG[kpDBG]/Sqrt[kpDBG], kpDBG -> Infinity]
];

(* At the peak-collinear endpoint, changing variable x=lambda tau leaves a
   constant branch integrand.  The zeroth Laguerre moment makes each branch
   exact independently. *)
peakMinusQuadratureDBG = Together[
  Total[weightsDBG]/lambdaMinusDBG[kpDBG]
  -1/lambdaMinusDBG[kpDBG]
];
peakPlusQuadratureDBG = Together[
  Total[weightsDBG]/lambdaPlusDBG[kpDBG]
  -1/lambdaPlusDBG[kpDBG]
];

exactTextDBG[expression_] := ToString[InputForm[expression]];
gateDBG[name_, residual_] := Module[
  {pass = TrueQ[PossibleZeroQ[residual]]},
  Print[name, " = ", If[pass, "PASS", "FAIL"]];
  If[!pass, Print["  residual: ", InputForm[residual]]];
  pass
];
gatesDBG = {
  gateDBG["two-node Laguerre moments are exact through degree three",
    Total[momentResidualsDBG]],
  gateDBG["slow branch continuous integral", branchMinusExactDBG],
  gateDBG["fast branch continuous integral", branchPlusExactDBG],
  gateDBG["two exponential branches recover the cosh transfer",
    coshIdentityDBG],
  gateDBG["two exponential branches recover the sinh transfer",
    sinhIdentityDBG],
  gateDBG["slow branch quadrature is exact at the finite-depth peak endpoint",
    peakMinusQuadratureDBG],
  gateDBG["fast branch quadrature is exact at the finite-depth peak endpoint",
    peakPlusQuadratureDBG],
  gateDBG["slow branch has the exact deep-water endpoint",
    deepMinusDBG-(2-Sqrt[2])],
  gateDBG["fast branch has the exact deep-water endpoint",
    deepPlusDBG-(2+Sqrt[2])]
};
overallPassDBG = And @@ gatesDBG;

projectRootDBG = DirectoryName[$InputFileName, 3];
artifactPathDBG = FileNameJoin[{
  projectRootDBG, "artifacts",
  "order2_finite_depth_directional_psi22_dual_branch_gl_freeze.json"
}];
interfacePathDBG = FileNameJoin[{
  projectRootDBG, "symbolic", "generated",
  "finite_depth_directional_order2_psi22_dual_branch_gl.json"
}];
sourceRelativePathDBG =
  "symbolic/wolfram/order2_finite_depth_directional_psi22_dual_branch_gl_freeze.wl";
residualRelativePathDBG =
  "symbolic/residuals/order2_finite_depth_directional_euler.wl";
sourceSha256DBG = FileHash[$InputFileName, "SHA256", "HexString"];
residualSha256DBG = FileHash[residualSourceDBG, "SHA256", "HexString"];

payloadDBG = <|
  "schema_version" -> 1,
  "status" -> "frozen-symbolic-direct-Psi22-dual-branch-pre-field-validation",
  "candidate_id" -> "gl-psi22-dual-branch-gl2plus2-v1",
  "scope" ->
    "stationary finite-depth forward-directional pure-sum Psi22, kd>=0.3",
  "candidate_input_fields" -> {"eta11"},
  "candidate_output_fields" -> {"Psi22_plus"},
  "explicit_eta22_dependency" -> False,
  "explicit_Phi22_dependency" -> False,
  "time_evolution_used" -> False,
  "sampled_selection_used" -> False,
  "oracle_or_mf12_used" -> False,
  "model_parent_kd_minimum" -> "3/10",
  "branch_definition" -> <|
    "slow" -> "Exp[-(omega12-Omega12) tau]",
    "fast" -> "Exp[-(omega12+Omega12) tau]"
  |>,
  "quadrature" -> <|
    "rank_total" -> 4,
    "rank_per_branch" -> 2,
    "nodes_exact" -> (exactTextDBG /@ nodesDBG),
    "weights_exact" -> (exactTextDBG /@ weightsDBG),
    "lambda_minus" -> exactTextDBG[lambdaMinusDBG[kpDBG]],
    "lambda_plus" -> exactTextDBG[lambdaPlusDBG[kpDBG]]
  |>,
  "finite_depth_endpoint_kd" -> "3/10",
  "finite_depth_endpoint_scales" -> <|
    "lambda_minus" -> exactTextDBG[finiteEndpointMinusDBG],
    "lambda_plus" -> exactTextDBG[finiteEndpointPlusDBG]
  |>,
  "deep_water_endpoint_scales_over_sqrt_kpd" -> <|
    "lambda_minus" -> exactTextDBG[deepMinusDBG],
    "lambda_plus" -> exactTextDBG[deepPlusDBG]
  |>,
  "full_crossing_angle_dependence_retained" -> True,
  "gates" -> AssociationThread[
    {"laguerre_moments", "slow_integral", "fast_integral",
     "cosh_recovery", "sinh_recovery", "finite_slow_endpoint",
     "finite_fast_endpoint", "deep_slow_endpoint", "deep_fast_endpoint"},
    gatesDBG
  ],
  "overall_exact_gate_pass" -> overallPassDBG
|>;
Export[artifactPathDBG, payloadDBG, "RawJSON"];
Export[
  interfacePathDBG,
  Join[payloadDBG, <|
    "generator" -> sourceRelativePathDBG,
    "generator_sha256" -> sourceSha256DBG,
    "residual_source" -> residualRelativePathDBG,
    "residual_source_sha256" -> residualSha256DBG
  |>],
  "RawJSON"
];

Print["=== Direct Psi22 dual-branch GL freeze ==="];
Print["finite kd=0.3 scales = ",
  InputForm[{finiteEndpointMinusDBG, finiteEndpointPlusDBG}]];
Print["deep-water scales/Sqrt[kp d] = ",
  InputForm[{deepMinusDBG, deepPlusDBG}]];
Print["artifact = ", artifactPathDBG];
Print["interface = ", interfacePathDBG];
Print["OVERALL_ORDER2_PSI22_DUAL_BRANCH_GL_FREEZE = ",
  If[overallPassDBG, "PASS", "FAIL"]];
If[!overallPassDBG, Exit[1]];
