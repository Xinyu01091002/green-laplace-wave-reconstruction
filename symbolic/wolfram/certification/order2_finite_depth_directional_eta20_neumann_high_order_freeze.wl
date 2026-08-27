(* ::Package:: *)
(* Exact extension of the frozen eta20 Neumann-resolvent grammar to R4/R6.

   This source imports the certified Euler residual, R1/R2 grammar, and fixed
   product compiler.  It changes no coefficient from numerical validation.
   For a truncation retaining powers through r^p, let L=p+1 and

     G_p = -Sum[B^j/A^(j+1), {j,0,p}],
     K_p = G_p N + v^(2L) c0 z^L,

   so the cross-multiplied inverse residual is -(B/A)^(p+1). *)

ClearAll["Global`*"];
compilerSourceFD2D20NRH = FileNameJoin[{DirectoryName[$InputFileName,2],
  "discovery","order2_finite_depth_directional_eta20_neumann_resolvent_compiler.wl"}];
If[!FileExistsQ[compilerSourceFD2D20NRH],
  Print["FAIL: certified Neumann compiler is missing."]; Exit[1]];
Get[compilerSourceFD2D20NRH];

ClearAll[
  neumannInverseFD2D20NRH, candidateFD2D20NRH,
  compiledTermsFD2D20NRH, compiledCandidateFD2D20NRH,
  parentTermsFD2D20NRH, countsFD2D20NRH,
  inverseIdentityFD2D20NRH, residualIdentityFD2D20NRH,
  angularEndpointFD2D20NRH, radialEndpointFD2D20NRH,
  differencePowerIdentityFD2D20NRH, radialPowerIdentityFD2D20NRH,
  abstractCompilerIdentityFD2D20NRH, compilerIdentityFD2D20NRH,
  gateFD2D20NRH, powersFD2D20NRH,
  recordsFD2D20NRH, gatesFD2D20NRH, overallFD2D20NRH,
  projectRootFD2D20NRH, artifactPathFD2D20NRH,
  interfacePathFD2D20NRH, exactTextFD2D20NRH
];

powersFD2D20NRH = {2,4,6};
neumannInverseFD2D20NRH[rank_Integer] := -Sum[
  bExactFD2D20NR^j/aExactFD2D20NR^(j+1), {j,0,rank-1}];
candidateFD2D20NRH[rank_Integer] :=
  neumannInverseFD2D20NRH[rank] nExactFD2D20NR
  +vExactFD2D20NR^(2 rank) c0ExactFD2D20NR zExactFD2D20NR^rank;

compiledTermsFD2D20NRH[rank_Integer] := Module[
  {fdTerms,fkTerms,correctionTerms},
  fdTerms = Table[
    compileDifferencePowerFD2D20NR[forcingDTermsCompilerFD2D20NR,2 index],
    {index,0,rank-1}];
  fkTerms = Table[
    compileDifferencePowerFD2D20NR[forcingKTermsCompilerFD2D20NR,2 index+1],
    {index,0,rank-1}];
  correctionTerms = compileRadialGapPowerFD2D20NR[2 rank];
  <|"fd"->fdTerms,"fk"->fkTerms,"correction"->correctionTerms|>
];

compiledCandidateFD2D20NRH[rank_Integer] := Module[
  {terms,result,index},
  terms = compiledTermsFD2D20NRH[rank];
  result = -kernelFromTermsFD2D20NR[terms["fd"][[1]]]/2;
  For[index=1,index<=rank-1,index++,
    result = result+Cosh[qDeltaFD2D20]^index/aExactFD2D20NR^index (
      I kernelFromTermsFD2D20NR[terms["fk"][[index]]]/2
      -kernelFromTermsFD2D20NR[terms["fd"][[index+1]]]/2);
  ];
  result = result+I Cosh[qDeltaFD2D20]^rank
    /(2 aExactFD2D20NR^rank)
    kernelFromTermsFD2D20NR[terms["fk"][[rank]]]
    +vExactFD2D20NR^(2 rank) c0ExactFD2D20NR/qDeltaFD2D20^(2 rank)
      kernelFromTermsFD2D20NR[terms["correction"]];
  result /. compilerRulesFD2D20NR
];

parentTermsFD2D20NRH[rank_Integer] := Module[{terms},
  terms=compiledTermsFD2D20NRH[rank];
  Join[Flatten[terms["fd"],1],Flatten[terms["fk"],1],terms["correction"]]
];
countsFD2D20NRH[rank_Integer] := Module[
  {parentFilters,products,productFfts,transforms},
  parentFilters=16 rank-1;
  products=11 rank^2+8 rank+1;
  productFfts=rank+2;
  transforms=17 rank+2;
  <|"parent_filters"->parentFilters,"products"->products,
    "product_ffts"->productFfts,"transforms"->transforms,
    "output_filters"->rank+2|>
];

inverseIdentityFD2D20NRH[rank_Integer] := Module[{aa,bb},
  Together[(bb-aa)(-Sum[bb^j/aa^(j+1),{j,0,rank-1}])-1
    +(bb/aa)^rank]
];
residualIdentityFD2D20NRH[rank_Integer] := Module[{aa,bb,nn,cc},
  Together[(bb-aa)(-Sum[bb^j/aa^(j+1),{j,0,rank-1}]nn+cc)-nn
    -(-(bb/aa)^rank nn+(bb-aa)cc)]
];
angularEndpointFD2D20NRH[rank_Integer] := Module[{aa,nn},
  (* Equal radial magnitude gives B=0 and z=0 for every nonzero angle. *)
  Together[(-aa)(-nn/aa)-nn]
];
radialEndpointFD2D20NRH[rank_Integer] := Module[{vv,cc},
  (* The certified one-sided limit is B/A=v^2 and z=1. *)
  Together[-(vv^2)^rank cc+vv^(2 rank)cc]
];
differencePowerIdentityFD2D20NRH[terms_,power_Integer] :=
  Module[{xx,yy},
    Expand[Sum[(-1)^right Binomial[power,right]
      xx^(power-right)yy^right,{right,0,power}]-(xx-yy)^power] === 0
    && Length[compileDifferencePowerFD2D20NR[terms,power]]
      ==Length[terms](power+1)
  ];
radialPowerIdentityFD2D20NRH[power_Integer] := Module[{xx,yy},
  Expand[Sum[(-1)^right Binomial[power,right]
    xx^(power-right)yy^right,{right,0,power}]-(xx-yy)^power] === 0
  && Length[compileRadialGapPowerFD2D20NR[power]]==power+1
];
abstractCompilerIdentityFD2D20NRH[rank_Integer] := Module[
  {aa,bb,cc,ss,fd,fk,nn,target,compiled,index},
  bb=ss^2 cc;
  nn=(aa fd-I ss cc fk)/2;
  target=-Sum[bb^j/aa^(j+1),{j,0,rank-1}] nn;
  compiled=-fd/2;
  For[index=1,index<=rank-1,index++,
    compiled=compiled+cc^index/aa^index (
      I ss^(2 index-1)fk/2-ss^(2 index)fd/2)];
  compiled=compiled+I cc^rank ss^(2 rank-1)fk/(2 aa^rank);
  Together[compiled-target] === 0
];
compilerIdentityFD2D20NRH[rank_Integer] := And[
  abstractCompilerIdentityFD2D20NRH[rank],
  And@@Table[differencePowerIdentityFD2D20NRH[
    forcingDTermsCompilerFD2D20NR,2 index],{index,0,rank-1}],
  And@@Table[differencePowerIdentityFD2D20NRH[
    forcingKTermsCompilerFD2D20NR,2 index+1],{index,0,rank-1}],
  radialPowerIdentityFD2D20NRH[2 rank]
];

gateFD2D20NRH[name_,value_] := Module[{pass=TrueQ[value]},
  Print[name," = ",If[pass,"PASS","FAIL"]];pass];
Print["Building exact R2/R4/R6 even-power records and fixed-graph counts."];
Do[With[{layers=power+1},Print["R",power," compiler subgates = ",InputForm[{
  abstractCompilerIdentityFD2D20NRH[layers],
  And@@Table[differencePowerIdentityFD2D20NRH[
    forcingDTermsCompilerFD2D20NR,2 index],{index,0,layers-1}],
  And@@Table[differencePowerIdentityFD2D20NRH[
    forcingKTermsCompilerFD2D20NR,2 index+1],{index,0,layers-1}],
  radialPowerIdentityFD2D20NRH[2 layers]}]]],{power,powersFD2D20NRH}];
recordsFD2D20NRH = Association@Table[With[{layers=power+1},
  ToString[power] -> <|
    "candidate_id"->("neumann-eta20-r"<>ToString[power]<>"-series-v1"),
    "neumann_max_power"->power,
    "neumann_layers"->layers,
    "kernel"->ToString[InputForm[candidateFD2D20NRH[layers]]],
    "counts"->countsFD2D20NRH[layers],
    "inverse_identity_pass"->TrueQ[PossibleZeroQ[inverseIdentityFD2D20NRH[layers]]],
    "residual_identity_pass"->TrueQ[PossibleZeroQ[residualIdentityFD2D20NRH[layers]]],
    "angular_endpoint_pass"->TrueQ[PossibleZeroQ[angularEndpointFD2D20NRH[layers]]],
    "radial_endpoint_pass"->TrueQ[PossibleZeroQ[radialEndpointFD2D20NRH[layers]]],
    "compiler_identity_pass"->TrueQ[compilerIdentityFD2D20NRH[layers]]
  |>], {power,powersFD2D20NRH}];
gatesFD2D20NRH = Table[
  gateFD2D20NRH["R"<>ToString[power]<>" exact and compiler gates",
    And@@Lookup[recordsFD2D20NRH[ToString[power]],{
      "inverse_identity_pass","residual_identity_pass",
      "angular_endpoint_pass","radial_endpoint_pass","compiler_identity_pass"}]],
  {power,powersFD2D20NRH}];
overallFD2D20NRH = And@@gatesFD2D20NRH;

projectRootFD2D20NRH=DirectoryName[$InputFileName,4];
artifactPathFD2D20NRH=FileNameJoin[{projectRootFD2D20NRH,"artifacts",
  "order2_finite_depth_directional_eta20_neumann_high_order_freeze.json"}];
interfacePathFD2D20NRH=FileNameJoin[{projectRootFD2D20NRH,"symbolic","generated",
  "eta20_neumann_r_series.json"}];
payloadFD2D20NRH=<|
  "schema_version"->1,
  "status"->"frozen-symbolic-neumann-even-power-R2-R4-R6-pre-validation",
  "scope"->"strict-nonzero finite-depth directional eta20",
  "candidate_input_fields"->{"eta11"},
  "oracle_or_mf12_used"->False,
  "sampled_selection_used"->False,
  "strict_zero_mode"->"separate and set to zero",
  "source_compiler"->FileNameDrop[compilerSourceFD2D20NRH,4],
  "source_compiler_sha256"->FileHash[compilerSourceFD2D20NRH,"SHA256","HexString"],
  "generator"->FileNameDrop[$InputFileName,4],
  "generator_sha256"->FileHash[$InputFileName,"SHA256","HexString"],
  "records"->recordsFD2D20NRH,
  "overall_exact_gate_pass"->overallFD2D20NRH|>;
Export[artifactPathFD2D20NRH,payloadFD2D20NRH,"RawJSON"];
Export[interfacePathFD2D20NRH,payloadFD2D20NRH,"RawJSON"];
Do[Print["R",power," counts = ",
  InputForm[countsFD2D20NRH[power+1]]],{power,powersFD2D20NRH}];
Print["interface = ",interfacePathFD2D20NRH];
Print["OVERALL_ORDER2_ETA20_NEUMANN_HIGH_ORDER_FREEZE = ",
  If[overallFD2D20NRH,"PASS","FAIL"]];
If[!overallFD2D20NRH,Exit[1]];
