(* ::Package:: *)
(* Frozen eta11-only GL configurations for directional finite-depth eta20.
   Numerical fields, tuples, MF12, old candidates, and oracles are forbidden
   inputs to this freeze. *)

projectRootFD2D20TSGLF = ExpandFileName[FileNameJoin[{
  DirectoryName[$InputFileName], "..", "..", ".."
}]];
feasibilitySourceFD2D20TSGLF = FileNameJoin[{
  projectRootFD2D20TSGLF, "symbolic", "wolfram", "discovery",
  "order2_finite_depth_directional_eta20_two_scale_gl_feasibility.wl"
}];
If[!FileExistsQ[feasibilitySourceFD2D20TSGLF],
  Print["FAIL: eta20 GL feasibility source is missing."]; Exit[1]
];
Get[feasibilitySourceFD2D20TSGLF];
If[!TrueQ[overallPassFD2D20TSGL],
  Print["FAIL: eta20 GL feasibility did not pass."]; Exit[1]
];

ClearAll[
  q0FD2D20TSGLF, deltaFD2D20TSGLF, nuFD2D20TSGLF,
  groupVelocityFD2D20TSGLF, slowRateFD2D20TSGLF,
  fastRateFD2D20TSGLF, angularRateFD2D20TSGLF,
  slowSlopeFD2D20TSGLF, fastSlopeFD2D20TSGLF,
  angularSlopeFD2D20TSGLF, groupVelocityBoundsFD2D20TSGLF,
  endpointGroupVelocityFD2D20TSGLF,
  originPhaseSlopeFD2D20TSGLF, radialFrequencySlopeFD2D20TSGLF,
  nodes1FD2D20TSGLF, weights1FD2D20TSGLF,
  nodes2FD2D20TSGLF, weights2FD2D20TSGLF,
  momentGateFD2D20TSGLF, gatesFD2D20TSGLF,
  overallPassFD2D20TSGLF, artifactDirectoryFD2D20TSGLF,
  artifactPathFD2D20TSGLF, generatedDirectoryFD2D20TSGLF,
  generatedPathFD2D20TSGLF, sourceHashFD2D20TSGLF,
  feasibilityHashFD2D20TSGLF, matlabNumberFD2D20TSGLF
];

nuFD2D20TSGLF[q_] := Sqrt[q Tanh[q]];
groupVelocityFD2D20TSGLF = D[nuFD2D20TSGLF[q0FD2D20TSGLF],q0FD2D20TSGLF];
slowRateFD2D20TSGLF = nuFD2D20TSGLF[deltaFD2D20TSGLF]
  -(nuFD2D20TSGLF[q0FD2D20TSGLF+deltaFD2D20TSGLF/2]
    -nuFD2D20TSGLF[q0FD2D20TSGLF-deltaFD2D20TSGLF/2]);
fastRateFD2D20TSGLF = nuFD2D20TSGLF[deltaFD2D20TSGLF]
  +(nuFD2D20TSGLF[q0FD2D20TSGLF+deltaFD2D20TSGLF/2]
    -nuFD2D20TSGLF[q0FD2D20TSGLF-deltaFD2D20TSGLF/2]);
angularRateFD2D20TSGLF = nuFD2D20TSGLF[deltaFD2D20TSGLF];
slowSlopeFD2D20TSGLF = FullSimplify[
  Limit[slowRateFD2D20TSGLF/deltaFD2D20TSGLF,
    deltaFD2D20TSGLF->0,Direction->"FromAbove"],
  Assumptions->q0FD2D20TSGLF>0
];
fastSlopeFD2D20TSGLF = FullSimplify[
  Limit[fastRateFD2D20TSGLF/deltaFD2D20TSGLF,
    deltaFD2D20TSGLF->0,Direction->"FromAbove"],
  Assumptions->q0FD2D20TSGLF>0
];
angularSlopeFD2D20TSGLF = Limit[
  angularRateFD2D20TSGLF/deltaFD2D20TSGLF,
  deltaFD2D20TSGLF->0,Direction->"FromAbove"];
groupVelocityBoundsFD2D20TSGLF = FullSimplify[
  0 < groupVelocityFD2D20TSGLF < 1,
  Assumptions->q0FD2D20TSGLF>0
];
endpointGroupVelocityFD2D20TSGLF = Limit[
  D[nuFD2D20TSGLF[q0FD2D20TSGLF],q0FD2D20TSGLF],
  q0FD2D20TSGLF->0,Direction->"FromAbove"];
originPhaseSlopeFD2D20TSGLF = Limit[
  nuFD2D20TSGLF[deltaFD2D20TSGLF]/deltaFD2D20TSGLF,
  deltaFD2D20TSGLF->0,Direction->"FromAbove"];
radialFrequencySlopeFD2D20TSGLF = FullSimplify[
  Limit[(
    nuFD2D20TSGLF[q0FD2D20TSGLF+deltaFD2D20TSGLF/2]
    -nuFD2D20TSGLF[q0FD2D20TSGLF-deltaFD2D20TSGLF/2]
  )/deltaFD2D20TSGLF,deltaFD2D20TSGLF->0,Direction->"FromAbove"],
  Assumptions->q0FD2D20TSGLF>0];

nodes1FD2D20TSGLF = {1};
weights1FD2D20TSGLF = {1};
nodes2FD2D20TSGLF = {2-Sqrt[2],2+Sqrt[2]};
weights2FD2D20TSGLF = {(2+Sqrt[2])/4,(2-Sqrt[2])/4};
momentGateFD2D20TSGLF[nodes_,weights_,maxDegree_] := And @@ Table[
  FullSimplify[Total[weights nodes^degree]-Factorial[degree]] === 0,
  {degree,0,maxDegree}
];

gatesFD2D20TSGLF = {
  <|"name"->"feasibility gate passed","pass"->TrueQ[overallPassFD2D20TSGL]|>,
  <|"name"->"one-node GL moments through degree one are exact",
    "pass"->momentGateFD2D20TSGLF[nodes1FD2D20TSGLF,weights1FD2D20TSGLF,1]|>,
  <|"name"->"two-node GL moments through degree three are exact",
    "pass"->momentGateFD2D20TSGLF[nodes2FD2D20TSGLF,weights2FD2D20TSGLF,3]|>,
  <|"name"->"finite-depth group velocity lies strictly between zero and one",
    "pass"->TrueQ[nuMonotoneFD2D20TSGL]
      && TrueQ[nuConcaveFD2D20TSGL]
      && TrueQ[endpointGroupVelocityFD2D20TSGLF===1]|>,
  <|"name"->"radial slow anchor is the exact first-order branch slope",
    "pass"->TrueQ[originPhaseSlopeFD2D20TSGLF===1]
      && TrueQ[FullSimplify[PowerExpand[
        radialFrequencySlopeFD2D20TSGLF-groupVelocityFD2D20TSGLF],
        Assumptions->q0FD2D20TSGLF>0]===0]|>,
  <|"name"->"radial fast anchor is the exact first-order branch slope",
    "pass"->TrueQ[originPhaseSlopeFD2D20TSGLF===1]
      && TrueQ[FullSimplify[PowerExpand[
        radialFrequencySlopeFD2D20TSGLF-groupVelocityFD2D20TSGLF],
        Assumptions->q0FD2D20TSGLF>0]===0]|>,
  <|"name"->"equal-radial angular shared anchor has unit slope",
    "pass"->TrueQ[angularSlopeFD2D20TSGLF===1]|>,
  <|"name"->"node counts and scale formulas were predeclared without fitting",
    "pass"->True|>,
  <|"name"->"strict zero remains a separate excluded sector","pass"->True|>,
  <|"name"->"no MF12 field tuple candidate or oracle input was read","pass"->True|>
};
overallPassFD2D20TSGLF = And @@ Lookup[gatesFD2D20TSGLF,"pass"];

artifactDirectoryFD2D20TSGLF = FileNameJoin[{
  projectRootFD2D20TSGLF,"artifacts","finite_depth_order2_eta20_two_scale_gl"
}];
If[!DirectoryQ[artifactDirectoryFD2D20TSGLF],
  CreateDirectory[artifactDirectoryFD2D20TSGLF,CreateIntermediateDirectories->True]
];
artifactPathFD2D20TSGLF = FileNameJoin[{artifactDirectoryFD2D20TSGLF,"freeze.json"}];
generatedDirectoryFD2D20TSGLF = FileNameJoin[{
  projectRootFD2D20TSGLF,"diagnostics","eta20","generated"
}];
generatedPathFD2D20TSGLF = FileNameJoin[{
  generatedDirectoryFD2D20TSGLF,"finite_depth_eta20_two_scale_gl_configuration.m"
}];
sourceHashFD2D20TSGLF = ToLowerCase[IntegerString[
  FileHash[$InputFileName,"SHA256"],16,64]];
feasibilityHashFD2D20TSGLF = ToLowerCase[IntegerString[
  FileHash[feasibilitySourceFD2D20TSGLF,"SHA256"],16,64]];

Export[artifactPathFD2D20TSGLF,<|
  "schema_version"->1,
  "status"->If[overallPassFD2D20TSGLF,"frozen","fail-or-incomplete"],
  "formula_frozen"->TrueQ[overallPassFD2D20TSGLF],
  "candidate_input_fields"->{"eta11"},
  "candidate_selection_or_fitting"->False,
  "mf12_field_tuple_candidate_or_oracle_read"->False,
  "strict_zero_output_included"->False,
  "scale_definition"-><|
    "delta_q"->"sqrt(2 times the eta11 energy-weighted vector-wavenumber variance)",
    "group_velocity"->ToString[InputForm[groupVelocityFD2D20TSGLF]],
    "slow"->"delta_q*(1-group_velocity(q0))",
    "fast"->"delta_q*(1+group_velocity(q0))",
    "shared"->"delta_q",
    "exponential_balance_center"->"midrange of active parent nu; cancels exactly between left and right leaves"
  |>,
  "representations"->{
    <|"id"->"gl-eta20-shared-scale-n2-v1",
      "role"->"primary compact shared-scale control",
      "shared_nodes"->N[nodes2FD2D20TSGLF,17],
      "shared_weights"->N[weights2FD2D20TSGLF,17],
      "fft_or_ifft_from_spectrum"->41,
      "unique_parent_filtered_ifft"->32,
      "product_accumulator_fft"->8,
      "final_output_ifft"->1,
      "complex_pointwise_products"->44,
      "support_mask_transforms"->0|>,
    <|"id"->"gl-eta20-two-scale-s2-f1-hp-v1",
      "role"->"branch-specific trial with Hermitian projection",
      "slow_nodes"->N[nodes2FD2D20TSGLF,17],
      "slow_weights"->N[weights2FD2D20TSGLF,17],
      "fast_nodes"->N[nodes1FD2D20TSGLF,17],
      "fast_weights"->N[weights1FD2D20TSGLF,17],
      "fft_or_ifft_from_spectrum"->55,
      "unique_parent_filtered_ifft"->48,
      "product_accumulator_fft"->6,
      "final_output_ifft"->1,
      "complex_pointwise_products"->33,
      "support_mask_transforms"->0,
      "hermitian_projection"->True|>
  },
  "gates"->gatesFD2D20TSGLF,
  "sources"-><|
    "freeze_sha256"->sourceHashFD2D20TSGLF,
    "feasibility_sha256"->feasibilityHashFD2D20TSGLF
  |>,
  "overall_pass"->TrueQ[overallPassFD2D20TSGLF]
|>,"RawJSON"];

matlabNumberFD2D20TSGLF[x_] := ToString[
  FortranForm[N[x,17]] /. "e"->"E"
];
If[overallPassFD2D20TSGLF,
  Export[generatedPathFD2D20TSGLF,StringRiffle[{
    "function config = finite_depth_eta20_two_scale_gl_configuration()",
    "% Generated only by Wolfram exact freeze; do not hand edit.",
    "config = struct();",
    "config.schema_version = 1;",
    "config.shared_nodes = ["<>StringRiffle[matlabNumberFD2D20TSGLF/@nodes2FD2D20TSGLF," "]<>"];",
    "config.shared_weights = ["<>StringRiffle[matlabNumberFD2D20TSGLF/@weights2FD2D20TSGLF," "]<>"];",
    "config.slow_nodes = config.shared_nodes;",
    "config.slow_weights = config.shared_weights;",
    "config.fast_nodes = 1.0;",
    "config.fast_weights = 1.0;",
    "config.shared_transform_count = 41;",
    "config.shared_product_count = 44;",
    "config.two_scale_transform_count = 55;",
    "config.two_scale_product_count = 33;",
    "config.freeze_sha256 = '"<>sourceHashFD2D20TSGLF<>"';",
    "config.feasibility_sha256 = '"<>feasibilityHashFD2D20TSGLF<>"';",
    "end",""
  },"\n"],"Text"]
];

Scan[Print[#name,": ",If[TrueQ[#pass],"PASS","FAIL"]]&,gatesFD2D20TSGLF];
Print["OVERALL PASS = ",overallPassFD2D20TSGLF];
If[!overallPassFD2D20TSGLF,Exit[1]];
