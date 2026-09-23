(* ::Package:: *)
(* Exact named-field freeze of the order-four lower-state forcing graph.

   This file obtains the forcing from the authoritative Taylor-expanded Euler
   residual and proves that the declared 1+1+1+1, 2+1+1, 2+2, and 3+1
   MATLAB operator blocks sum to it exactly.  It exports no target eta44
   kernel and reads no tuple, MF12 field, waveform, or validation oracle.
*)

residualSourceFD4FG = FileNameJoin[{
  DirectoryName[$InputFileName], "..", "residuals",
  "order1_to_order4_symbolic_expansion_2d.wl"
}];
If[!FileExistsQ[residualSourceFD4FG],
  Print["FAIL: authoritative order-four residual source missing."];
  Exit[1]
];
Get[residualSourceFD4FG];

ClearAll[
  eFD4FG, pFD4FG, dot2FD4FG, dot3FD4FG,
  k1111FD4FG, k211FD4FG, k22FD4FG, k31FD4FG,
  d1111FD4FG, d211FD4FG, d22FD4FG, d31FD4FG,
  kDeclaredFD4FG, dDeclaredFD4FG, exactTextFD4FG,
  gateFD4FG, gatesFD4FG, overallPassFD4FG,
  projectRootFD4FG, generatedDirectoryFD4FG, interfacePathFD4FG,
  sourceRelativePathFD4FG, sourceSha256FD4FG, residualRelativePathFD4FG,
  residualSha256FD4FG
];

eFD4FG[n_, ax_:0, ay_:0, at_:0] := D[
  etaG2[n][xG2, yG2, tG2], {xG2, ax}, {yG2, ay}, {tG2, at}
];
pFD4FG[n_, ax_:0, ay_:0, az_:0, at_:0] := (
  D[phiG2[n][xG2, yG2, zG2, tG2],
    {xG2, ax}, {yG2, ay}, {zG2, az}, {tG2, at}] /. zG2 -> 0
);
dot2FD4FG[ax_, ay_, bx_, by_] := ax bx + ay by;
dot3FD4FG[ax_, ay_, az_, bx_, by_, bz_] := ax bx + ay by + az bz;

(* 3+1 partition. *)
k31FD4FG = Expand[
  dot2FD4FG[pFD4FG[1,1], pFD4FG[1,0,1],
    eFD4FG[3,1], eFD4FG[3,0,1]]
  + dot2FD4FG[pFD4FG[3,1], pFD4FG[3,0,1],
    eFD4FG[1,1], eFD4FG[1,0,1]]
  - eFD4FG[1] pFD4FG[3,0,0,2]
  - eFD4FG[3] pFD4FG[1,0,0,2]
];
d31FD4FG = Expand[
  eFD4FG[1] pFD4FG[3,0,0,1,1]
  + eFD4FG[3] pFD4FG[1,0,0,1,1]
  + dot3FD4FG[
      pFD4FG[1,1], pFD4FG[1,0,1], pFD4FG[1,0,0,1],
      pFD4FG[3,1], pFD4FG[3,0,1], pFD4FG[3,0,0,1]
    ]
];

(* 2+2 partition. *)
k22FD4FG = Expand[
  dot2FD4FG[pFD4FG[2,1], pFD4FG[2,0,1],
    eFD4FG[2,1], eFD4FG[2,0,1]]
  - eFD4FG[2] pFD4FG[2,0,0,2]
];
d22FD4FG = Expand[
  eFD4FG[2] pFD4FG[2,0,0,1,1]
  + 1/2 dot3FD4FG[
      pFD4FG[2,1], pFD4FG[2,0,1], pFD4FG[2,0,0,1],
      pFD4FG[2,1], pFD4FG[2,0,1], pFD4FG[2,0,0,1]
    ]
];

(* 2+1+1 partition. *)
k211FD4FG = Expand[
  eFD4FG[1] dot2FD4FG[
    pFD4FG[1,1,0,1], pFD4FG[1,0,1,1],
    eFD4FG[2,1], eFD4FG[2,0,1]]
  + dot2FD4FG[
      eFD4FG[1] pFD4FG[2,1,0,1]
        + eFD4FG[2] pFD4FG[1,1,0,1],
      eFD4FG[1] pFD4FG[2,0,1,1]
        + eFD4FG[2] pFD4FG[1,0,1,1],
      eFD4FG[1,1], eFD4FG[1,0,1]
    ]
  - eFD4FG[1]^2 pFD4FG[2,0,0,3]/2
  - eFD4FG[1] eFD4FG[2] pFD4FG[1,0,0,3]
];
d211FD4FG = Expand[
  eFD4FG[1]^2 pFD4FG[2,0,0,2,1]/2
  + eFD4FG[1] eFD4FG[2] pFD4FG[1,0,0,2,1]
  + dot3FD4FG[
      pFD4FG[1,1], pFD4FG[1,0,1], pFD4FG[1,0,0,1],
      eFD4FG[1] pFD4FG[2,1,0,1]
        + eFD4FG[2] pFD4FG[1,1,0,1],
      eFD4FG[1] pFD4FG[2,0,1,1]
        + eFD4FG[2] pFD4FG[1,0,1,1],
      eFD4FG[1] pFD4FG[2,0,0,2]
        + eFD4FG[2] pFD4FG[1,0,0,2]
    ]
  + dot3FD4FG[
      pFD4FG[2,1], pFD4FG[2,0,1], pFD4FG[2,0,0,1],
      eFD4FG[1] pFD4FG[1,1,0,1],
      eFD4FG[1] pFD4FG[1,0,1,1],
      eFD4FG[1] pFD4FG[1,0,0,2]
    ]
];

(* 1+1+1+1 partition. *)
k1111FD4FG = Expand[
  eFD4FG[1]^2/2 dot2FD4FG[
    pFD4FG[1,1,0,2], pFD4FG[1,0,1,2],
    eFD4FG[1,1], eFD4FG[1,0,1]]
  - eFD4FG[1]^3 pFD4FG[1,0,0,4]/6
];
d1111FD4FG = Expand[
  eFD4FG[1]^3 pFD4FG[1,0,0,3,1]/6
  + eFD4FG[1]^2/2 dot3FD4FG[
      pFD4FG[1,1], pFD4FG[1,0,1], pFD4FG[1,0,0,1],
      pFD4FG[1,1,0,2], pFD4FG[1,0,1,2], pFD4FG[1,0,0,3]
    ]
  + eFD4FG[1]^2/2 dot3FD4FG[
      pFD4FG[1,1,0,1], pFD4FG[1,0,1,1], pFD4FG[1,0,0,2],
      pFD4FG[1,1,0,1], pFD4FG[1,0,1,1], pFD4FG[1,0,0,2]
    ]
];

kDeclaredFD4FG = Expand[k1111FD4FG + k211FD4FG + k22FD4FG + k31FD4FG];
dDeclaredFD4FG = Expand[d1111FD4FG + d211FD4FG + d22FD4FG + d31FD4FG];
exactTextFD4FG[expression_] := ToString[InputForm[expression]];
gateFD4FG[name_, residual_] := Module[
  {pass = TrueQ[PossibleZeroQ[Expand[residual]]]},
  Print[name, " = ", If[pass, "PASS", "FAIL"]];
  If[!pass, Print["  residual: ", InputForm[Expand[residual]]]];
  pass
];
gatesFD4FG = {
  gateFD4FG["declared kinematic graph equals Euler forcing",
    kDeclaredFD4FG - kinematicForcingOrder4G2],
  gateFD4FG["declared dynamic graph equals Euler forcing",
    dDeclaredFD4FG - dynamicForcingOrder4G2]
};
overallPassFD4FG = And @@ gatesFD4FG;

projectRootFD4FG = DirectoryName[$InputFileName, 3];
generatedDirectoryFD4FG = FileNameJoin[{
  projectRootFD4FG, "symbolic", "generated"
}];
If[!DirectoryQ[generatedDirectoryFD4FG],
  CreateDirectory[generatedDirectoryFD4FG, CreateIntermediateDirectories -> True]
];
interfacePathFD4FG = FileNameJoin[{
  generatedDirectoryFD4FG,
  "finite_depth_directional_order4_forcing_graph.json"
}];
sourceRelativePathFD4FG =
  "symbolic/wolfram/order4_finite_depth_directional_forcing_graph_freeze.wl";
residualRelativePathFD4FG =
  "symbolic/residuals/order1_to_order4_symbolic_expansion_2d.wl";
sourceSha256FD4FG = FileHash[$InputFileName, "SHA256", "HexString"];
residualSha256FD4FG = FileHash[residualSourceFD4FG, "SHA256", "HexString"];

Export[
  interfacePathFD4FG,
  <|
    "schema_version" -> 1,
    "generator" -> sourceRelativePathFD4FG,
    "generator_sha256" -> sourceSha256FD4FG,
    "residual_source" -> residualRelativePathFD4FG,
    "residual_source_sha256" -> residualSha256FD4FG,
    "status" -> "exact-order4-lower-state-forcing-graph-frozen",
    "oracle_or_mf12_used" -> False,
    "target_eta44_kernel_used" -> False,
    "partition_blocks" -> {"1+1+1+1", "2+1+1", "2+2", "3+1"},
    "analytic_positive_field_weights" -> <|
      "1+1+1+1" -> "1/8", "2+1+1" -> "1/4",
      "2+2" -> "1/2", "3+1" -> "1/2"
    |>,
    "kinematic_blocks_exact" -> <|
      "1+1+1+1" -> exactTextFD4FG[k1111FD4FG],
      "2+1+1" -> exactTextFD4FG[k211FD4FG],
      "2+2" -> exactTextFD4FG[k22FD4FG],
      "3+1" -> exactTextFD4FG[k31FD4FG]
    |>,
    "dynamic_blocks_exact" -> <|
      "1+1+1+1" -> exactTextFD4FG[d1111FD4FG],
      "2+1+1" -> exactTextFD4FG[d211FD4FG],
      "2+2" -> exactTextFD4FG[d22FD4FG],
      "3+1" -> exactTextFD4FG[d31FD4FG]
    |>,
    "exact_gates" -> gatesFD4FG,
    "overall_exact_gate_pass" -> overallPassFD4FG
  |>,
  "RawJSON"
];

Print["interface = ", interfacePathFD4FG];
Print["OVERALL_ORDER4_FORCING_GRAPH_FREEZE = ",
  If[overallPassFD4FG, "PASS", "FAIL"]];
If[!overallPassFD4FG, Exit[1]];
