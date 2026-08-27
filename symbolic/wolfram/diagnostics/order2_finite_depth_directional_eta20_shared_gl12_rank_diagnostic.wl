(* ::Package:: *)
(* Post-baseline pure shared-scale GL12 rank diagnostic for eta20.
   The frozen forcing graph and scale are unchanged. GL12 is the standard
   Gauss--Laguerre rule and is not selected or fitted from MF12 fields. *)

ClearAll["Global`*"];
projectRootFD2D20GL12 = DirectoryName[$InputFileName,4];
baseSourceFD2D20GL12 = FileNameJoin[{projectRootFD2D20GL12,"symbolic","wolfram",
  "discovery","order2_finite_depth_directional_eta20_two_scale_gl_freeze.wl"}];
If[!FileExistsQ[baseSourceFD2D20GL12],
  Print["FAIL: base eta20 shared GL freeze is missing."]; Exit[1]];
Get[baseSourceFD2D20GL12];
If[!TrueQ[overallPassFD2D20TSGLF],
  Print["FAIL: base eta20 shared GL freeze did not pass."]; Exit[1]];

nFD2D20GL12 = 12;
xFD2D20GL12 = Unique["x"];
polyFD2D20GL12 = LaguerreL[nFD2D20GL12,xFD2D20GL12];
orthogonalityFD2D20GL12 = And[
  Exponent[polyFD2D20GL12,xFD2D20GL12]===nFD2D20GL12,
  And@@Table[Integrate[Exp[-xFD2D20GL12] xFD2D20GL12^degree
    polyFD2D20GL12,{xFD2D20GL12,0,Infinity},GenerateConditions->False]===0,
    {degree,0,nFD2D20GL12-1}]
];
nodesFD2D20GL12 = Sort[N[xFD2D20GL12 /. NSolve[
  polyFD2D20GL12==0,xFD2D20GL12,Reals,WorkingPrecision->90],75]];
weightsFD2D20GL12 = N[nodesFD2D20GL12/((nFD2D20GL12+1)^2
  (LaguerreL[nFD2D20GL12+1,#]&/@nodesFD2D20GL12)^2),75];
momentResidualFD2D20GL12 = Max[Abs@Table[
  Total[weightsFD2D20GL12 nodesFD2D20GL12^degree]-Factorial[degree],
  {degree,0,2nFD2D20GL12-1}]];
positiveFD2D20GL12 = (
  And@@Positive[nodesFD2D20GL12]
  && And@@Positive[weightsFD2D20GL12]
);
overallFD2D20GL12 = (
  TrueQ[orthogonalityFD2D20GL12]
  && TrueQ[positiveFD2D20GL12]
  && momentResidualFD2D20GL12<10^-10
);

sourceHashFD2D20GL12 = ToLowerCase[IntegerString[
  FileHash[$InputFileName,"SHA256"],16,64]];
baseHashFD2D20GL12 = ToLowerCase[IntegerString[
  FileHash[baseSourceFD2D20GL12,"SHA256"],16,64]];
artifactDirectoryFD2D20GL12 = FileNameJoin[{projectRootFD2D20GL12,"artifacts",
  "finite_depth_order2_eta20_two_scale_gl"}];
If[!DirectoryQ[artifactDirectoryFD2D20GL12],CreateDirectory[
  artifactDirectoryFD2D20GL12,CreateIntermediateDirectories->True]];
artifactPathFD2D20GL12 = FileNameJoin[{artifactDirectoryFD2D20GL12,
  "shared_gl12_rank_diagnostic.json"}];
generatedPathFD2D20GL12 = FileNameJoin[{projectRootFD2D20GL12,"diagnostics",
  "eta20","generated","finite_depth_eta20_shared_gl12_diagnostic_configuration.m"}];

matlabNumberFD2D20GL12[value_] := ToString[FortranForm[N[value,17]]];
matlabVectorFD2D20GL12[values_] := "["<>StringRiffle[
  matlabNumberFD2D20GL12/@values," "]<>"]";
linesFD2D20GL12 = {
  "function config = finite_depth_eta20_shared_gl12_diagnostic_configuration()",
  "% Wolfram-generated standard GL12 rank diagnostic; do not hand edit.",
  "config = struct();",
  "config.schema_version = 1;",
  "config.gl12_nodes = "<>matlabVectorFD2D20GL12[nodesFD2D20GL12]<>";",
  "config.gl12_weights = "<>matlabVectorFD2D20GL12[weightsFD2D20GL12]<>";",
  "config.gl12_transform_count = 435;",
  "config.gl12_product_count = 264;",
  "config.freeze_sha256 = '"<>sourceHashFD2D20GL12<>"';",
  "config.base_freeze_sha256 = '"<>baseHashFD2D20GL12<>"';",
  "config.field_selection_or_fitting = false;",
  "end",""
};
If[overallFD2D20GL12,Export[generatedPathFD2D20GL12,
  StringRiffle[linesFD2D20GL12,"\n"],"Text"]];
Export[artifactPathFD2D20GL12,<|
  "schema_version"->1,
  "status"->"post-baseline-standard-gl12-rank-diagnostic",
  "formula_frozen"->TrueQ[overallFD2D20GL12],
  "candidate_input_fields"->{"eta11"},
  "shared_scale"->"delta_q_rms_pair unchanged",
  "quadrature_rank"->12,
  "field_selection_or_fitting"->False,
  "mf12_field_tuple_candidate_or_oracle_read"->False,
  "nodes"->N[nodesFD2D20GL12,18],
  "weights"->N[weightsFD2D20GL12,18],
  "exact_moment_degree_by_gauss_theorem"->23,
  "high_precision_max_moment_residual"->ToString[InputForm[momentResidualFD2D20GL12]],
  "source_sha256"->sourceHashFD2D20GL12,
  "base_freeze_sha256"->baseHashFD2D20GL12,
  "overall_pass"->TrueQ[overallFD2D20GL12]
|>,"RawJSON"];
Print["GL12 exact Laguerre orthogonality: ",orthogonalityFD2D20GL12];
Print["GL12 positive nodes and weights: ",positiveFD2D20GL12];
Print["GL12 max moment residual: ",InputForm[momentResidualFD2D20GL12]];
Print["OVERALL_ETA20_SHARED_GL12_RANK_DIAGNOSTIC = ",
  If[overallFD2D20GL12,"PASS","FAIL"]];
If[!overallFD2D20GL12,Exit[1]];
