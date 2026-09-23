(* ::Package:: *)
(* Freeze the independent finite-depth directional eta33 Two-Scale and
   Shared-Scale Green--Laplace representation.  This file imports no
   candidate, tuple, waveform, MF12 field, or validation oracle. *)

ClearAll[
  omegaTSG, glRuleTSG, scaledMomentTSG, momentGateTSG,
  omegaSumTSG, omegaOutputTSG, resolventTSG, splitTSG,
  qStarTSG, pairSlowTSG, pairFastTSG, tripleSlowTSG,
  tripleFastTSG, pairSharedTSG, tripleSharedTSG,
  primaryCountsTSG, sharedCountsTSG, projectRootTSG,
  sourceFilesTSG, sourceHashesTSG, exportPathTSG,
  exactStringTSG, numericStringTSG, ruleRecordTSG,
  gateTSG, gatesTSG, overallTSG
];

omegaTSG[q_] := Sqrt[q Tanh[q]];

(* Exact Golub--Welsch-free construction.  The roots remain algebraic
   numbers and the moment gates are evaluated with exact Root objects. *)
glRuleTSG[n_Integer?Positive] := Module[{x, roots, weights},
  roots = x /. Solve[LaguerreL[n, x] == 0, x, Reals];
  roots = SortBy[roots, N];
  weights = (#/((n + 1)^2 LaguerreL[n + 1, #]^2)) & /@ roots;
  Transpose[{roots, weights}]
];

scaledMomentTSG[rule_, degree_Integer?NonNegative, mu_] :=
  Total[(#[[2]] #[[1]]^degree) & /@ rule]/mu^(degree + 1);

momentGateTSG[n_Integer?Positive] := Module[{rule = glRuleTSG[n], mu},
  And @@ Table[
    TrueQ[RootReduce[
      Total[(#[[2]] #[[1]]^degree) & /@ rule] - degree!
    ] === 0],
    {degree, 0, 2 n - 1}
  ]
];

omegaSumTSG = Symbol["Omega"];
omegaOutputTSG = Symbol["omegaK"];
resolventTSG = 1/(omegaSumTSG^2 - omegaOutputTSG^2);
splitTSG = 1/(2 omegaOutputTSG) (
  1/(omegaSumTSG - omegaOutputTSG)
  - 1/(omegaSumTSG + omegaOutputTSG)
);

qStarTSG = 1;
pairSlowTSG = 2 omegaTSG[qStarTSG] - omegaTSG[2 qStarTSG];
pairFastTSG = 2 omegaTSG[qStarTSG] + omegaTSG[2 qStarTSG];
tripleSlowTSG = 3 omegaTSG[qStarTSG] - omegaTSG[3 qStarTSG];
tripleFastTSG = 3 omegaTSG[qStarTSG] + omegaTSG[3 qStarTSG];
pairSharedTSG = Sqrt[pairSlowTSG pairFastTSG];
tripleSharedTSG = Sqrt[tripleSlowTSG tripleFastTSG];

primaryCountsTSG = <|
  "order2" -> <|"slow" -> 4, "fast" -> 2|>,
  "order3" -> <|"slow" -> 4, "fast" -> 2|>
|>;
sharedCountsTSG = <|"order2" -> 3, "order3" -> 3|>;

gateTSG[name_, value_] := Module[{pass = TrueQ[value]},
  Print[name <> " = " <> If[pass, "PASS", "FAIL"]];
  pass
];

gatesTSG = {
  gateTSG["exact slow-fast resolvent factorization",
    Together[resolventTSG - splitTSG] === 0],
  gateTSG["analytic-signal cubic normalization",
    2 (1/2)^3 == 1/4 && 6/24 == 1/4],
  gateTSG["GL n=2 moments through degree 3", momentGateTSG[2]],
  gateTSG["GL n=3 moments through degree 5", momentGateTSG[3]],
  gateTSG["GL n=4 moments through degree 7", momentGateTSG[4]],
  gateTSG["pair branch anchors are positive",
    And @@ Thread[N[{pairSlowTSG, pairFastTSG}, 50] > 0]],
  gateTSG["triple branch anchors are positive",
    And @@ Thread[N[{tripleSlowTSG, tripleFastTSG}, 50] > 0]]
};
overallTSG = And @@ gatesTSG;

projectRootTSG = DirectoryName[$InputFileName, 3];
sourceFilesTSG = {
  FileNameJoin[{projectRootTSG, "symbolic", "residuals",
    "order1_to_order3_symbolic_expansion_2d.wl"}],
  FileNameJoin[{projectRootTSG, "symbolic", "residuals",
    "order2_finite_depth_directional_euler.wl"}],
  FileNameJoin[{projectRootTSG, "symbolic", "residuals",
    "order3_finite_depth_directional_euler.wl"}]
};
sourceHashesTSG = AssociationThread[
  FileNameTake /@ sourceFilesTSG,
  IntegerString[FileHash[#, "SHA256"], 16, 64] & /@ sourceFilesTSG
];

exactStringTSG[x_] := ToString[InputForm[x]];
numericStringTSG[x_] := ToString[InputForm[N[x, 34]]];
ruleRecordTSG[n_] := Map[
  <|
    "node" -> numericStringTSG[#[[1]]],
    "weight" -> numericStringTSG[#[[2]]]
  |> &,
  glRuleTSG[n]
];

exportPathTSG = FileNameJoin[{
  projectRootTSG, "symbolic", "generated",
  "finite_depth_directional_order3_two_scale_gl.json"
}];
If[!DirectoryQ[DirectoryName[exportPathTSG]],
  CreateDirectory[DirectoryName[exportPathTSG], CreateIntermediateDirectories -> True]
];
Export[exportPathTSG, <|
  "schema_version" -> 1,
  "status" -> "frozen_pre_validation",
  "candidate" -> "finite_depth_directional_eta33_two_scale_gl",
  "scope" -> "stationary strict-forward positive-pure-sum eta33 q_parent>0.5",
  "candidate_input_fields" -> {"eta11"},
  "external_lower_order_fields" -> False,
  "oracle_imported" -> False,
  "mf12_imported" -> False,
  "q_star" -> 1,
  "primary_node_counts" -> primaryCountsTSG,
  "shared_control_node_counts" -> sharedCountsTSG,
  "scales" -> <|
    "order2_slow_exact" -> exactStringTSG[pairSlowTSG],
    "order2_fast_exact" -> exactStringTSG[pairFastTSG],
    "order2_shared_exact" -> exactStringTSG[pairSharedTSG],
    "order3_slow_exact" -> exactStringTSG[tripleSlowTSG],
    "order3_fast_exact" -> exactStringTSG[tripleFastTSG],
    "order3_shared_exact" -> exactStringTSG[tripleSharedTSG],
    "order2_slow" -> numericStringTSG[pairSlowTSG],
    "order2_fast" -> numericStringTSG[pairFastTSG],
    "order2_shared" -> numericStringTSG[pairSharedTSG],
    "order3_slow" -> numericStringTSG[tripleSlowTSG],
    "order3_fast" -> numericStringTSG[tripleFastTSG],
    "order3_shared" -> numericStringTSG[tripleSharedTSG]
  |>,
  "gauss_laguerre" -> <|
    "n2" -> ruleRecordTSG[2],
    "n3" -> ruleRecordTSG[3],
    "n4" -> ruleRecordTSG[4]
  |>,
  "operator_graph" -> <|
    "outer_eta_numerator" -> "G(K) F_D + i d_tau F_K",
    "outer_phi_numerator" -> "i d_tau F_D - F_K",
    "two_scale" -> "(slow integral-fast integral)/(2 omega_K)",
    "shared_scale" -> "integral sinh(omega_K tau)/omega_K"
  |>,
  "analytic_signal_normalization" -> <|
    "input" -> "U=2 eta11_positive",
    "internal_input_factor" -> "1/2",
    "output_factor" -> "2",
    "net_cubic_factor" -> "1/4",
    "ordered_permutations_over_grouped_factor" -> "6/24=1/4"
  |>,
  "source_sha256" -> sourceHashesTSG,
  "moment_certification" -> <|"n2" -> 3, "n3" -> 5, "n4" -> 7|>,
  "claim_boundary" ->
    "finite GL moments only; no uniform shallow-water or all-angle convergence claim"
|>, "RawJSON"];

Print["artifact = ", exportPathTSG];
Print["source sha256 = ", InputForm[sourceHashesTSG]];
Print["branch scales = ", InputForm[N[{
  pairSlowTSG, pairFastTSG, pairSharedTSG,
  tripleSlowTSG, tripleFastTSG, tripleSharedTSG
}, 20]]];
Print["OVERALL_ORDER3_TWO_SCALE_GL_FREEZE = ",
  If[overallTSG, "PASS", "FAIL"]];
If[!overallTSG, Exit[1]];
