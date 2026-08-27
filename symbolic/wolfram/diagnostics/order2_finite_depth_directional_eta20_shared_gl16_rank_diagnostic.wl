(* ::Package:: *)
(* Post-baseline standard GL16 rank diagnostic for the frozen eta20 graph. *)

ClearAll["Global`*"];
projectRootFD2D20GL16 = DirectoryName[$InputFileName,4];
baseSourceFD2D20GL16 = FileNameJoin[{projectRootFD2D20GL16,"symbolic","wolfram",
  "discovery","order2_finite_depth_directional_eta20_two_scale_gl_freeze.wl"}];
If[!FileExistsQ[baseSourceFD2D20GL16],Print["FAIL: base freeze missing."];Exit[1]];
Get[baseSourceFD2D20GL16];
If[!TrueQ[overallPassFD2D20TSGLF],Print["FAIL: base freeze failed."];Exit[1]];

nFD2D20GL16=16;xFD2D20GL16=Unique["x"];
polyFD2D20GL16=LaguerreL[nFD2D20GL16,xFD2D20GL16];
orthogonalityFD2D20GL16=And[
  Exponent[polyFD2D20GL16,xFD2D20GL16]===nFD2D20GL16,
  And@@Table[Integrate[Exp[-xFD2D20GL16] xFD2D20GL16^degree
    polyFD2D20GL16,{xFD2D20GL16,0,Infinity},GenerateConditions->False]===0,
    {degree,0,nFD2D20GL16-1}]
];
nodesFD2D20GL16=Sort[N[xFD2D20GL16/.NSolve[polyFD2D20GL16==0,
  xFD2D20GL16,Reals,WorkingPrecision->100],80]];
weightsFD2D20GL16=N[nodesFD2D20GL16/((nFD2D20GL16+1)^2
  (LaguerreL[nFD2D20GL16+1,#]&/@nodesFD2D20GL16)^2),80];
momentResidualFD2D20GL16=Max[Abs@Table[
  Total[weightsFD2D20GL16 nodesFD2D20GL16^degree]-Factorial[degree],
  {degree,0,2nFD2D20GL16-1}]];
positiveFD2D20GL16=(And@@Positive[nodesFD2D20GL16]
  &&And@@Positive[weightsFD2D20GL16]);
overallFD2D20GL16=(TrueQ[orthogonalityFD2D20GL16]
  &&TrueQ[positiveFD2D20GL16]&&momentResidualFD2D20GL16<10^-10);
sourceHashFD2D20GL16=ToLowerCase[IntegerString[FileHash[$InputFileName,"SHA256"],16,64]];
baseHashFD2D20GL16=ToLowerCase[IntegerString[FileHash[baseSourceFD2D20GL16,"SHA256"],16,64]];
artifactDirectoryFD2D20GL16=FileNameJoin[{projectRootFD2D20GL16,"artifacts",
  "finite_depth_order2_eta20_two_scale_gl"}];
If[!DirectoryQ[artifactDirectoryFD2D20GL16],CreateDirectory[
  artifactDirectoryFD2D20GL16,CreateIntermediateDirectories->True]];
artifactPathFD2D20GL16=FileNameJoin[{artifactDirectoryFD2D20GL16,
  "shared_gl16_rank_diagnostic.json"}];
generatedPathFD2D20GL16=FileNameJoin[{projectRootFD2D20GL16,"diagnostics",
  "eta20","generated","finite_depth_eta20_shared_gl16_diagnostic_configuration.m"}];
matlabNumberFD2D20GL16[value_]:=ToString[FortranForm[N[value,17]]];
matlabVectorFD2D20GL16[values_]:="["<>StringRiffle[
  matlabNumberFD2D20GL16/@values," "]<>"]";
linesFD2D20GL16={
  "function config = finite_depth_eta20_shared_gl16_diagnostic_configuration()",
  "% Wolfram-generated standard GL16 rank diagnostic; do not hand edit.",
  "config = struct();","config.schema_version = 1;",
  "config.gl16_nodes = "<>matlabVectorFD2D20GL16[nodesFD2D20GL16]<>";",
  "config.gl16_weights = "<>matlabVectorFD2D20GL16[weightsFD2D20GL16]<>";",
  "config.gl16_transform_count = 579;","config.gl16_product_count = 352;",
  "config.freeze_sha256 = '"<>sourceHashFD2D20GL16<>"';",
  "config.base_freeze_sha256 = '"<>baseHashFD2D20GL16<>"';",
  "config.field_selection_or_fitting = false;","end",""};
If[overallFD2D20GL16,Export[generatedPathFD2D20GL16,
  StringRiffle[linesFD2D20GL16,"\n"],"Text"]];
Export[artifactPathFD2D20GL16,<|
  "schema_version"->1,"status"->"post-baseline-standard-gl16-rank-diagnostic",
  "formula_frozen"->TrueQ[overallFD2D20GL16],"candidate_input_fields"->{"eta11"},
  "shared_scale"->"delta_q_rms_pair unchanged","quadrature_rank"->16,
  "field_selection_or_fitting"->False,
  "mf12_field_tuple_candidate_or_oracle_read"->False,
  "nodes"->N[nodesFD2D20GL16,18],"weights"->N[weightsFD2D20GL16,18],
  "exact_moment_degree_by_gauss_theorem"->31,
  "high_precision_max_moment_residual"->ToString[InputForm[momentResidualFD2D20GL16]],
  "source_sha256"->sourceHashFD2D20GL16,"base_freeze_sha256"->baseHashFD2D20GL16,
  "overall_pass"->TrueQ[overallFD2D20GL16]|>,"RawJSON"];
Print["GL16 exact Laguerre orthogonality: ",orthogonalityFD2D20GL16];
Print["GL16 positive nodes and weights: ",positiveFD2D20GL16];
Print["GL16 max moment residual: ",InputForm[momentResidualFD2D20GL16]];
Print["OVERALL_ETA20_SHARED_GL16_RANK_DIAGNOSTIC = ",
  If[overallFD2D20GL16,"PASS","FAIL"]];
If[!overallFD2D20GL16,Exit[1]];
