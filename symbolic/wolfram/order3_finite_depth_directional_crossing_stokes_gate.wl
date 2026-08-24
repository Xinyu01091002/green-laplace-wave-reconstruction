(* ::Package:: *)
(* Exact crossing-aware Stokes endpoint gate for eta33.

   This file freezes only the permutation-symmetric angular lift used by a
   later nested-resolvent core. It does not construct the core, choose a
   quadrature rank or scale, or use MF12, tuple samples, input-chi formulas,
   or an eta33 oracle.

   For gij=(1+ei.ej)/2, the three-parent gate is
     G3 = g12 g23 g31.
   It preserves the co-directional Stokes endpoint, vanishes if any parent
   pair is strictly opposed, and reduces to g(gamma)^2 for AAB/ABB two-lobe
   crossing triples.
*)

ClearAll[
  c12FD3G, c23FD3G, c31FD3G, cFD3G,
  x1FD3G, y1FD3G, x2FD3G, y2FD3G, x3FD3G, y3FD3G,
  variablesFD3G, abstractGateFD3G, componentGateFD3G,
  separatedTermsFD3G, termRecordFD3G, permutationRulesFD3G,
  permutationResidualsFD3G, stokesTraceFD3G, qFD3G,
  exactTextFD3G, gateFD3G, booleanGateFD3G, gatesFD3G, overallPassFD3G,
  projectRootFD3G, artifactPathFD3G, generatedDirectoryFD3G,
  interfacePathFD3G, sourceRelativePathFD3G, sourceSha256FD3G
];

abstractGateFD3G = Times[
  (1 + c12FD3G)/2,
  (1 + c23FD3G)/2,
  (1 + c31FD3G)/2
];
componentGateFD3G = Expand[abstractGateFD3G /. {
  c12FD3G -> x1FD3G x2FD3G + y1FD3G y2FD3G,
  c23FD3G -> x2FD3G x3FD3G + y2FD3G y3FD3G,
  c31FD3G -> x3FD3G x1FD3G + y3FD3G y1FD3G
}];
variablesFD3G = {
  x1FD3G, y1FD3G, x2FD3G, y2FD3G, x3FD3G, y3FD3G
};
separatedTermsFD3G = CoefficientRules[componentGateFD3G, variablesFD3G];
termRecordFD3G[rule_] := With[
  {powers = First[rule], coefficient = Last[rule]},
  <|
    "coefficient_exact" -> exactTextFD3G[coefficient],
    "coefficient_numerator" -> Numerator[coefficient],
    "coefficient_denominator" -> Denominator[coefficient],
    "leaf1_exponents" -> powers[[1 ;; 2]],
    "leaf2_exponents" -> powers[[3 ;; 4]],
    "leaf3_exponents" -> powers[[5 ;; 6]]
  |>
];

(* Pair-dot permutations induced by all six parent permutations. *)
permutationRulesFD3G = {
  {},
  {c12FD3G -> c12FD3G, c23FD3G -> c31FD3G, c31FD3G -> c23FD3G},
  {c12FD3G -> c23FD3G, c23FD3G -> c12FD3G, c31FD3G -> c31FD3G},
  {c12FD3G -> c31FD3G, c23FD3G -> c23FD3G, c31FD3G -> c12FD3G},
  {c12FD3G -> c23FD3G, c23FD3G -> c31FD3G, c31FD3G -> c12FD3G},
  {c12FD3G -> c31FD3G, c23FD3G -> c12FD3G, c31FD3G -> c23FD3G}
};
permutationResidualsFD3G = Expand[
  abstractGateFD3G - (abstractGateFD3G /. #)
] & /@ permutationRulesFD3G;

stokesTraceFD3G = Times[
  3 qFD3G^2/256,
  14 + 15 Cosh[2 qFD3G] + 6 Cosh[4 qFD3G] + Cosh[6 qFD3G],
  Csch[qFD3G]^6
];

exactTextFD3G[expression_] := ToString[InputForm[expression]];
gateFD3G[name_, residual_] := Module[
  {pass = TrueQ[PossibleZeroQ[residual]]},
  Print[name, " = ", If[pass, "PASS", "FAIL"]];
  If[!pass, Print["  residual: ", InputForm[residual]]];
  pass
];
booleanGateFD3G[name_, value_] := Module[{pass = TrueQ[value]},
  Print[name, " = ", If[pass, "PASS", "FAIL"]];
  pass
];
gatesFD3G = {
  booleanGateFD3G[
    "all parent permutations preserve G3",
    And @@ (PossibleZeroQ /@ permutationResidualsFD3G)
  ],
  gateFD3G[
    "co-directional Stokes endpoint is preserved",
    (abstractGateFD3G /. {
      c12FD3G -> 1, c23FD3G -> 1, c31FD3G -> 1
    }) - 1
  ],
  gateFD3G[
    "AAB crossing reduces to the order-two g squared gate",
    Together[
      (abstractGateFD3G /. {
        c12FD3G -> 1, c23FD3G -> cFD3G, c31FD3G -> cFD3G
      }) - ((1 + cFD3G)/2)^2
    ]
  ],
  gateFD3G[
    "any strictly opposed parent pair removes the endpoint repair",
    (abstractGateFD3G /. c12FD3G -> -1)
  ],
  gateFD3G[
    "component expansion exactly reproduces G3",
    Expand[
      componentGateFD3G - (abstractGateFD3G /. {
        c12FD3G -> x1FD3G x2FD3G + y1FD3G y2FD3G,
        c23FD3G -> x2FD3G x3FD3G + y2FD3G y3FD3G,
        c31FD3G -> x3FD3G x1FD3G + y3FD3G y1FD3G
      })
    ]
  ]
};
overallPassFD3G = And @@ gatesFD3G;

projectRootFD3G = DirectoryName[$InputFileName, 3];
artifactPathFD3G = FileNameJoin[{
  projectRootFD3G, "artifacts",
  "order3_finite_depth_directional_crossing_stokes_gate.json"
}];
sourceRelativePathFD3G =
  "symbolic/wolfram/order3_finite_depth_directional_crossing_stokes_gate.wl";
sourceSha256FD3G = FileHash[$InputFileName, "SHA256", "HexString"];
Export[
  artifactPathFD3G,
  <|
    "schema_version" -> 1,
    "status" -> "frozen-angular-gate-structure-core-amplitude-pending",
    "scope" -> "finite-depth forward pure-sum eta33 crossing endpoint repair",
    "candidate_constructed" -> False,
    "input_chi_route_used" -> False,
    "oracle_or_mf12_used" -> False,
    "sampled_selection_used" -> False,
    "gate_exact" -> exactTextFD3G[abstractGateFD3G],
    "physical_form" ->
      "G3=Product_{i<j} (1+unit_ki dot unit_kj)/2",
    "two_lobe_AAB_or_ABB_reduction" -> "((1+cos(gamma))/2)^2",
    "ordered_stokes_trace" -> exactTextFD3G[stokesTraceFD3G],
    "future_correction_amplitude" ->
      "C3(q)=D_eta3(q)-K_nested_core_diagonal(q)",
    "separated_term_count" -> Length[separatedTermsFD3G],
    "maximum_direction_degree_per_leaf" -> 2,
    "separated_terms" -> (termRecordFD3G /@ separatedTermsFD3G),
    "gates" -> AssociationThread[
      {"permutation", "codirectional", "two_lobe", "opposed", "separated"},
      gatesFD3G
    ],
    "overall_exact_gate_pass" -> overallPassFD3G
  |>,
  "RawJSON"
];

generatedDirectoryFD3G = FileNameJoin[{
  projectRootFD3G, "symbolic", "generated"
}];
If[!DirectoryQ[generatedDirectoryFD3G],
  CreateDirectory[generatedDirectoryFD3G, CreateIntermediateDirectories -> True]
];
interfacePathFD3G = FileNameJoin[{
  generatedDirectoryFD3G,
  "finite_depth_directional_order3_crossing_stokes_gate.json"
}];
Export[
  interfacePathFD3G,
  <|
    "schema_version" -> 1,
    "generator" -> sourceRelativePathFD3G,
    "generator_sha256" -> sourceSha256FD3G,
    "status" -> "frozen-angular-gate-structure-core-amplitude-pending",
    "input_chi_route_used" -> False,
    "oracle_or_mf12_used" -> False,
    "gate_exact" -> exactTextFD3G[abstractGateFD3G],
    "separated_term_count" -> Length[separatedTermsFD3G],
    "maximum_direction_degree_per_leaf" -> 2,
    "separated_terms" -> (termRecordFD3G /@ separatedTermsFD3G),
    "overall_exact_gate_pass" -> overallPassFD3G
  |>,
  "RawJSON"
];

Print["separated term count = ", Length[separatedTermsFD3G]];
Print["artifact = ", artifactPathFD3G];
Print["interface = ", interfacePathFD3G];
Print[
  "OVERALL_ORDER3_CROSSING_STOKES_GATE = ",
  If[overallPassFD3G, "PASS", "FAIL"]
];
If[!overallPassFD3G, Exit[1]];
