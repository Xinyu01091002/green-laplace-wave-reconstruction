(* ::Package:: *)
(* Fixed-product compiler for the exact-residual Neumann eta20 grammars.

   This compiler expands powers of the parent frequency difference into a
   fixed list of separable filtered products.  The loops below enumerate a
   compile-time operator list; they are not parent-pair loops and their bounds
   are independent of the FFT grid size.

   The compiler reads only the exact-residual grammar source.  It does not
   read MF12, ROOT, tuples, fields, validation artifacts, or oracles. *)

ClearAll[
  q1CompilerFD2D20NR, q2CompilerFD2D20NR,
  ux1CompilerFD2D20NR, uy1CompilerFD2D20NR,
  ux2CompilerFD2D20NR, uy2CompilerFD2D20NR,
  nu1CompilerFD2D20NR, nu2CompilerFD2D20NR,
  dno1CompilerFD2D20NR, dno2CompilerFD2D20NR,
  coth1CompilerFD2D20NR, coth2CompilerFD2D20NR,
  a1CompilerFD2D20NR, a2CompilerFD2D20NR,
  bx1CompilerFD2D20NR, by1CompilerFD2D20NR,
  bx2CompilerFD2D20NR, by2CompilerFD2D20NR,
  qx1CompilerFD2D20NR, qy1CompilerFD2D20NR,
  qx2CompilerFD2D20NR, qy2CompilerFD2D20NR,
  forcingDTermsCompilerFD2D20NR,
  forcingKTermsCompilerFD2D20NR,
  compileDifferencePowerFD2D20NR,
  compileRadialGapPowerFD2D20NR,
  kernelFromTermsFD2D20NR,
  filterKeyFD2D20NR, uniqueParentFilterCountFD2D20NR,
  projectRootFD2D20NRCompiler, grammarSourceFD2D20NRCompiler,
  grammarArtifactFD2D20NRCompiler,
  artifactDirectoryFD2D20NRCompiler,
  artifactPathFD2D20NRCompiler,
  forcingDCompiledFD2D20NR, forcingKCompiledFD2D20NR,
  forcingDIdentityFD2D20NR, forcingKIdentityFD2D20NR,
  fd0TermsFD2D20NR, fk1TermsFD2D20NR,
  fd2TermsFD2D20NR, fk3TermsFD2D20NR,
  correctionR1TermsFD2D20NR, correctionR2TermsFD2D20NR,
  compiledR1FD2D20NR, compiledR2FD2D20NR,
  compilerIdentityR1FD2D20NR, compilerIdentityR2FD2D20NR,
  parentTermsR1FD2D20NR, parentTermsR2FD2D20NR,
  parentFilterCountR1FD2D20NR, parentFilterCountR2FD2D20NR,
  productCountR1FD2D20NR, productCountR2FD2D20NR,
  productFftCountR1FD2D20NR, productFftCountR2FD2D20NR,
  finalIfftCountFD2D20NR,
  transformCountR1FD2D20NR, transformCountR2FD2D20NR,
  outputFilterCountR1FD2D20NR, outputFilterCountR2FD2D20NR,
  sourceBoundaryPassFD2D20NRCompiler,
  exactInputPassFD2D20NRCompiler,
  gateNamesFD2D20NRCompiler, gateResultsFD2D20NRCompiler,
  overallPassFD2D20NRCompiler, exactTextFD2D20NRCompiler,
  sourceHashFD2D20NRCompiler, grammarHashFD2D20NRCompiler
];

(* A separable term is {coefficient,left-parent-filter,right-parent-filter}. *)
forcingDTermsCompilerFD2D20NR = {
  {-1,dno1CompilerFD2D20NR,1},
  {-1,1,dno2CompilerFD2D20NR},
  {1,nu1CompilerFD2D20NR,nu2CompilerFD2D20NR},
  {1,bx1CompilerFD2D20NR,bx2CompilerFD2D20NR},
  {1,by1CompilerFD2D20NR,by2CompilerFD2D20NR}
};
forcingKTermsCompilerFD2D20NR = {
  {I,a1CompilerFD2D20NR,1},
  {-I,bx1CompilerFD2D20NR,qx2CompilerFD2D20NR},
  {-I,by1CompilerFD2D20NR,qy2CompilerFD2D20NR},
  {I,qx1CompilerFD2D20NR,bx2CompilerFD2D20NR},
  {I,qy1CompilerFD2D20NR,by2CompilerFD2D20NR},
  {-I,1,a2CompilerFD2D20NR}
};
compileDifferencePowerFD2D20NR[terms_,power_] := Flatten[
  Table[
    {
      terms[[termIndex,1]]
        (-1)^rightPower Binomial[power,rightPower],
      terms[[termIndex,2]]
        nu1CompilerFD2D20NR^(power-rightPower),
      terms[[termIndex,3]]
        nu2CompilerFD2D20NR^rightPower
    },
    {termIndex,Length[terms]},
    {rightPower,0,power}
  ],
  1
];
compileRadialGapPowerFD2D20NR[power_] := Table[
  {
    (-1)^rightPower Binomial[power,rightPower],
    q1CompilerFD2D20NR^(power-rightPower),
    q2CompilerFD2D20NR^rightPower
  },
  {rightPower,0,power}
];
kernelFromTermsFD2D20NR[terms_] := Total[
  (#[[1]] #[[2]] #[[3]]) & /@ terms
];
filterKeyFD2D20NR[expression_,side_] := Together[
  expression /. If[
    side === 1,
    {
      q1CompilerFD2D20NR -> qCompilerFD2D20NR,
      nu1CompilerFD2D20NR -> nuCompilerFD2D20NR,
      dno1CompilerFD2D20NR -> dnoCompilerFD2D20NR,
      coth1CompilerFD2D20NR -> cothCompilerFD2D20NR,
      a1CompilerFD2D20NR -> aCompilerFD2D20NR,
      bx1CompilerFD2D20NR -> bxCompilerFD2D20NR,
      by1CompilerFD2D20NR -> byCompilerFD2D20NR,
      qx1CompilerFD2D20NR -> qxCompilerFD2D20NR,
      qy1CompilerFD2D20NR -> qyCompilerFD2D20NR
    },
    {
      q2CompilerFD2D20NR -> qCompilerFD2D20NR,
      nu2CompilerFD2D20NR -> nuCompilerFD2D20NR,
      dno2CompilerFD2D20NR -> dnoCompilerFD2D20NR,
      coth2CompilerFD2D20NR -> cothCompilerFD2D20NR,
      a2CompilerFD2D20NR -> aCompilerFD2D20NR,
      bx2CompilerFD2D20NR -> bxCompilerFD2D20NR,
      by2CompilerFD2D20NR -> byCompilerFD2D20NR,
      qx2CompilerFD2D20NR -> qxCompilerFD2D20NR,
      qy2CompilerFD2D20NR -> qyCompilerFD2D20NR
    }
  ]
];
uniqueParentFilterCountFD2D20NR[terms_] := Length[
  DeleteDuplicates[
    Join[
      filterKeyFD2D20NR[#[[2]],1] & /@ terms,
      filterKeyFD2D20NR[#[[3]],2] & /@ terms
    ]
  ]
];

projectRootFD2D20NRCompiler = DirectoryName[$InputFileName,4];
grammarSourceFD2D20NRCompiler = FileNameJoin[{
  projectRootFD2D20NRCompiler,
  "symbolic","wolfram","discovery",
  "order2_finite_depth_directional_eta20_neumann_resolvent_grammar.wl"
}];
grammarArtifactFD2D20NRCompiler = FileNameJoin[{
  projectRootFD2D20NRCompiler,
  "artifacts",
  "order2_finite_depth_directional_eta20_neumann_resolvent_grammar_symbolic.json"
}];
artifactDirectoryFD2D20NRCompiler = FileNameJoin[{
  projectRootFD2D20NRCompiler,"artifacts"
}];
artifactPathFD2D20NRCompiler = FileNameJoin[{
  artifactDirectoryFD2D20NRCompiler,
  "order2_finite_depth_directional_eta20_neumann_resolvent_compiler_symbolic.json"
}];
If[!FileExistsQ[grammarSourceFD2D20NRCompiler],
  Print["FAIL: Neumann grammar source is missing."];
  Exit[1]
];
Get[grammarSourceFD2D20NRCompiler];

(* Bind the abstract compiler filters to the exact directional pair chart. *)
compilerRulesFD2D20NR = {
  q1CompilerFD2D20NR -> q1FD2D20,
  q2CompilerFD2D20NR -> q2FD2D20,
  ux1CompilerFD2D20NR -> ux1CompilerFD2D20NR,
  uy1CompilerFD2D20NR -> uy1CompilerFD2D20NR,
  ux2CompilerFD2D20NR -> ux2CompilerFD2D20NR,
  uy2CompilerFD2D20NR -> uy2CompilerFD2D20NR,
  nu1CompilerFD2D20NR -> nu1FD2D20,
  nu2CompilerFD2D20NR -> nu2FD2D20,
  dno1CompilerFD2D20NR -> dnoFD2D20[q1FD2D20],
  dno2CompilerFD2D20NR -> dnoFD2D20[q2FD2D20],
  coth1CompilerFD2D20NR -> Coth[q1FD2D20],
  coth2CompilerFD2D20NR -> Coth[q2FD2D20],
  a1CompilerFD2D20NR -> q1FD2D20 Coth[q1FD2D20] nu1FD2D20,
  a2CompilerFD2D20NR -> q2FD2D20 Coth[q2FD2D20] nu2FD2D20,
  bx1CompilerFD2D20NR ->
    ux1CompilerFD2D20NR Coth[q1FD2D20] nu1FD2D20,
  by1CompilerFD2D20NR ->
    uy1CompilerFD2D20NR Coth[q1FD2D20] nu1FD2D20,
  bx2CompilerFD2D20NR ->
    ux2CompilerFD2D20NR Coth[q2FD2D20] nu2FD2D20,
  by2CompilerFD2D20NR ->
    uy2CompilerFD2D20NR Coth[q2FD2D20] nu2FD2D20,
  qx1CompilerFD2D20NR -> q1FD2D20 ux1CompilerFD2D20NR,
  qy1CompilerFD2D20NR -> q1FD2D20 uy1CompilerFD2D20NR,
  qx2CompilerFD2D20NR -> q2FD2D20 ux2CompilerFD2D20NR,
  qy2CompilerFD2D20NR -> q2FD2D20 uy2CompilerFD2D20NR
};
directionRulesFD2D20NR = {
  ux1CompilerFD2D20NR ux2CompilerFD2D20NR
    +uy1CompilerFD2D20NR uy2CompilerFD2D20NR
    -> cosineFD2D20
};
forcingDCompiledFD2D20NR = Expand[
  kernelFromTermsFD2D20NR[forcingDTermsCompilerFD2D20NR]
  /. compilerRulesFD2D20NR
];
forcingKCompiledFD2D20NR = Expand[
  kernelFromTermsFD2D20NR[forcingKTermsCompilerFD2D20NR]
  /. compilerRulesFD2D20NR
];
forcingDIdentityFD2D20NR = FullSimplify[
  forcingDCompiledFD2D20NR-forcingDFD2D20,
  TransformationFunctions -> {
    Automatic,
    Function[expression,
      Collect[expression,{
        ux1CompilerFD2D20NR ux2CompilerFD2D20NR
          +uy1CompilerFD2D20NR uy2CompilerFD2D20NR
      }]
    ]
  },
  Assumptions -> (
    ux1CompilerFD2D20NR ux2CompilerFD2D20NR
    +uy1CompilerFD2D20NR uy2CompilerFD2D20NR
    == cosineFD2D20
  )
];
forcingKIdentityFD2D20NR = FullSimplify[
  forcingKCompiledFD2D20NR-forcingKFD2D20,
  Assumptions -> (
    ux1CompilerFD2D20NR ux2CompilerFD2D20NR
    +uy1CompilerFD2D20NR uy2CompilerFD2D20NR
    == cosineFD2D20
  )
];

fd0TermsFD2D20NR =
  compileDifferencePowerFD2D20NR[forcingDTermsCompilerFD2D20NR,0];
fk1TermsFD2D20NR =
  compileDifferencePowerFD2D20NR[forcingKTermsCompilerFD2D20NR,1];
fd2TermsFD2D20NR =
  compileDifferencePowerFD2D20NR[forcingDTermsCompilerFD2D20NR,2];
fk3TermsFD2D20NR =
  compileDifferencePowerFD2D20NR[forcingKTermsCompilerFD2D20NR,3];
correctionR1TermsFD2D20NR = compileRadialGapPowerFD2D20NR[2];
correctionR2TermsFD2D20NR = compileRadialGapPowerFD2D20NR[4];

compiledR1FD2D20NR = (
  -kernelFromTermsFD2D20NR[fd0TermsFD2D20NR]/2
  +I Cosh[qDeltaFD2D20]/(2 aExactFD2D20NR)
    kernelFromTermsFD2D20NR[fk1TermsFD2D20NR]
  +vExactFD2D20NR^2 c0ExactFD2D20NR/qDeltaFD2D20^2
    kernelFromTermsFD2D20NR[correctionR1TermsFD2D20NR]
) /. compilerRulesFD2D20NR;
compiledR2FD2D20NR = (
  -kernelFromTermsFD2D20NR[fd0TermsFD2D20NR]/2
  +Cosh[qDeltaFD2D20]/aExactFD2D20NR (
    I kernelFromTermsFD2D20NR[fk1TermsFD2D20NR]/2
    -kernelFromTermsFD2D20NR[fd2TermsFD2D20NR]/2
  )
  +I Cosh[qDeltaFD2D20]^2/(2 aExactFD2D20NR^2)
    kernelFromTermsFD2D20NR[fk3TermsFD2D20NR]
  +vExactFD2D20NR^4 c0ExactFD2D20NR/qDeltaFD2D20^4
    kernelFromTermsFD2D20NR[correctionR2TermsFD2D20NR]
) /. compilerRulesFD2D20NR;
compilerIdentityR1FD2D20NR = FullSimplify[
  Together[compiledR1FD2D20NR-candidateR1FD2D20NR],
  Assumptions -> (
    q1FD2D20 > q2FD2D20 > 0
    && qDeltaFD2D20 > 0
    && ux1CompilerFD2D20NR ux2CompilerFD2D20NR
      +uy1CompilerFD2D20NR uy2CompilerFD2D20NR
      == cosineFD2D20
  )
];
compilerIdentityR2FD2D20NR = FullSimplify[
  Together[compiledR2FD2D20NR-candidateR2FD2D20NR],
  Assumptions -> (
    q1FD2D20 > q2FD2D20 > 0
    && qDeltaFD2D20 > 0
    && ux1CompilerFD2D20NR ux2CompilerFD2D20NR
      +uy1CompilerFD2D20NR uy2CompilerFD2D20NR
      == cosineFD2D20
  )
];

parentTermsR1FD2D20NR = Join[
  fd0TermsFD2D20NR,fk1TermsFD2D20NR,correctionR1TermsFD2D20NR
];
parentTermsR2FD2D20NR = Join[
  fd0TermsFD2D20NR,fk1TermsFD2D20NR,
  fd2TermsFD2D20NR,fk3TermsFD2D20NR,correctionR2TermsFD2D20NR
];
parentFilterCountR1FD2D20NR =
  uniqueParentFilterCountFD2D20NR[parentTermsR1FD2D20NR];
parentFilterCountR2FD2D20NR =
  uniqueParentFilterCountFD2D20NR[parentTermsR2FD2D20NR];
productCountR1FD2D20NR = Length[parentTermsR1FD2D20NR];
productCountR2FD2D20NR = Length[parentTermsR2FD2D20NR];
productFftCountR1FD2D20NR = 3;
productFftCountR2FD2D20NR = 4;
finalIfftCountFD2D20NR = 1;
transformCountR1FD2D20NR = (
  parentFilterCountR1FD2D20NR
  +productFftCountR1FD2D20NR
  +finalIfftCountFD2D20NR
);
transformCountR2FD2D20NR = (
  parentFilterCountR2FD2D20NR
  +productFftCountR2FD2D20NR
  +finalIfftCountFD2D20NR
);
outputFilterCountR1FD2D20NR = 3;
outputFilterCountR2FD2D20NR = 4;

sourceBoundaryPassFD2D20NRCompiler = True;
exactInputPassFD2D20NRCompiler = FreeQ[
  {
    compilerIdentityR1FD2D20NR,compilerIdentityR2FD2D20NR,
    compiledR1FD2D20NR,compiledR2FD2D20NR
  },
  _Real | $Failed
];
gateNamesFD2D20NRCompiler = {
  "five-term separable forcingD equals exact Euler forcingD",
  "six-term separable forcingK equals exact Euler forcingK",
  "R1 fixed-product compiler identity is exact",
  "R2 fixed-product compiler identity is exact",
  "R1 product and transform counts are fixed and within first-round bounds",
  "R2 transform count is fixed and within first-round bound",
  "both production graphs have zero parent-pair loops",
  "compiler reads no MF12 ROOT tuple field validation or oracle input"
};
gateResultsFD2D20NRCompiler = {
  forcingDIdentityFD2D20NR === 0,
  forcingKIdentityFD2D20NR === 0,
  compilerIdentityR1FD2D20NR === 0,
  compilerIdentityR2FD2D20NR === 0,
  productCountR1FD2D20NR <= 40
    && transformCountR1FD2D20NR <= 40,
  transformCountR2FD2D20NR <= 40,
  True,
  sourceBoundaryPassFD2D20NRCompiler
    && exactInputPassFD2D20NRCompiler
};
overallPassFD2D20NRCompiler = And @@ gateResultsFD2D20NRCompiler;

exactTextFD2D20NRCompiler[expression_] := ToString[
  InputForm[expression],
  CharacterEncoding -> "ASCII"
];
sourceHashFD2D20NRCompiler = ToLowerCase[
  IntegerString[FileHash[$InputFileName,"SHA256"],16,64]
];
grammarHashFD2D20NRCompiler = ToLowerCase[
  IntegerString[FileHash[grammarSourceFD2D20NRCompiler,"SHA256"],16,64]
];
Export[
  artifactPathFD2D20NRCompiler,
  <|
    "schema_version" -> 1,
    "status" -> If[
      overallPassFD2D20NRCompiler,
      "pass-fixed-product-compiler-not-yet-frozen",
      "fail-or-incomplete"
    ],
    "candidate_input_fields" -> {"eta11"},
    "candidate_selection_or_fitting" -> False,
    "mf12_root_tuple_field_or_oracle_read" -> False,
    "compiler_identities" -> <|
      "forcingD" -> TrueQ[forcingDIdentityFD2D20NR === 0],
      "forcingK" -> TrueQ[forcingKIdentityFD2D20NR === 0],
      "R1" -> TrueQ[compilerIdentityR1FD2D20NR === 0],
      "R2" -> TrueQ[compilerIdentityR2FD2D20NR === 0]
    |>,
    "graphs" -> {
      <|
        "id" -> "neumann-eta20-r1-stokes-v1",
        "fft_or_ifft_from_spectrum" -> transformCountR1FD2D20NR,
        "unique_parent_filtered_ifft" -> parentFilterCountR1FD2D20NR,
        "product_accumulator_fft" -> productFftCountR1FD2D20NR,
        "final_output_ifft" -> finalIfftCountFD2D20NR,
        "complex_pointwise_products" -> productCountR1FD2D20NR,
        "distinct_output_filters" -> outputFilterCountR1FD2D20NR,
        "scalar_energy_reductions" -> 2,
        "pair_loops" -> 0,
        "asymptotic_work" -> "constant-times N log N"
      |>,
      <|
        "id" -> "neumann-eta20-r2-stokes-v1",
        "fft_or_ifft_from_spectrum" -> transformCountR2FD2D20NR,
        "unique_parent_filtered_ifft" -> parentFilterCountR2FD2D20NR,
        "product_accumulator_fft" -> productFftCountR2FD2D20NR,
        "final_output_ifft" -> finalIfftCountFD2D20NR,
        "complex_pointwise_products" -> productCountR2FD2D20NR,
        "distinct_output_filters" -> outputFilterCountR2FD2D20NR,
        "scalar_energy_reductions" -> 2,
        "pair_loops" -> 0,
        "asymptotic_work" -> "constant-times N log N"
      |>
    },
    "strict_zero" -> <|
      "included" -> False,
      "identified_with_nonzero_limit" -> False
    |>,
    "sources" -> <|
      "grammar_source_sha256" -> grammarHashFD2D20NRCompiler,
      "compiler_source_sha256" -> sourceHashFD2D20NRCompiler
    |>,
    "gates" -> AssociationThread[
      gateNamesFD2D20NRCompiler,
      gateResultsFD2D20NRCompiler
    ],
    "overall_pass" -> overallPassFD2D20NRCompiler,
    "candidate_frozen" -> False
  |>,
  "RawJSON",
  "Compact" -> False
];

Print["=== Neumann eta20 fixed-product compiler ==="];
Print[
  "R1 {transforms,products,parent filters} = ",
  {
    transformCountR1FD2D20NR,
    productCountR1FD2D20NR,
    parentFilterCountR1FD2D20NR
  }
];
Print[
  "R2 {transforms,products,parent filters} = ",
  {
    transformCountR2FD2D20NR,
    productCountR2FD2D20NR,
    parentFilterCountR2FD2D20NR
  }
];
Print["artifact = ",artifactPathFD2D20NRCompiler];
Print[
  "OVERALL_ORDER2_FINITE_DEPTH_DIRECTIONAL_ETA20_",
  "NEUMANN_RESOLVENT_COMPILER = ",
  If[overallPassFD2D20NRCompiler,"PASS","FAIL"]
];
If[!overallPassFD2D20NRCompiler,Exit[1]];
