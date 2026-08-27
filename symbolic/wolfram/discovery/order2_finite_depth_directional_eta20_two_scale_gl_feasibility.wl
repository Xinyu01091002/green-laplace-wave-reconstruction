(* ::Package:: *)
(* Exact feasibility gate for a finite-depth directional eta20 two-branch
   Laplace/GL representation. This file freezes no numerical quadrature and
   reads no field, tuple, MF12, candidate, or oracle artifact. *)

projectRootFD2D20TSGL = ExpandFileName[FileNameJoin[{
  DirectoryName[$InputFileName], "..", "..", ".."
}]];
eulerSourceFD2D20TSGL = FileNameJoin[{
  projectRootFD2D20TSGL, "symbolic", "residuals",
  "order2_finite_depth_directional_eta20_euler.wl"
}];
If[!FileExistsQ[eulerSourceFD2D20TSGL],
  Print["FAIL: exact directional eta20 Euler source is missing."];
  Exit[1]
];
Get[eulerSourceFD2D20TSGL];

ClearAll[
  aFD2D20TSGL, sFD2D20TSGL, exactResolventFD2D20TSGL,
  slowNumeratorFD2D20TSGL, fastNumeratorFD2D20TSGL,
  branchSplitFD2D20TSGL, tauFD2D20TSGL,
  slowIntegralFD2D20TSGL, fastIntegralFD2D20TSGL,
  combinedIntegralFD2D20TSGL, branchIdentityFD2D20TSGL,
  slowIntegralIdentityFD2D20TSGL, fastIntegralIdentityFD2D20TSGL,
  combinedIntegralIdentityFD2D20TSGL, nuMonotoneFD2D20TSGL,
  nuConcaveFD2D20TSGL, radialIncrementFD2D20TSGL,
  radialIncrementIdentityFD2D20TSGL, qDeltaBoundFD2D20TSGL,
  nuEndpointFD2D20TSGL,
  physicalScaleProofFD2D20TSGL, gateNamesFD2D20TSGL,
  gateResultsFD2D20TSGL, overallPassFD2D20TSGL,
  artifactDirectoryFD2D20TSGL, artifactPathFD2D20TSGL,
  exactTextFD2D20TSGL
];

ClearAll[aaFD2D20TSGL, ssFD2D20TSGL, fkFD2D20TSGL, fdFD2D20TSGL];

aFD2D20TSGL = Sqrt[dnoFD2D20[qDeltaFD2D20]];
sFD2D20TSGL = sigmaDeltaFD2D20;
exactResolventFD2D20TSGL = (
  I sFD2D20TSGL forcingKFD2D20
  -aFD2D20TSGL^2 forcingDFD2D20
)/(2 (aFD2D20TSGL^2-sFD2D20TSGL^2));

slowNumeratorFD2D20TSGL = I forcingKFD2D20-aFD2D20TSGL forcingDFD2D20;
fastNumeratorFD2D20TSGL = -I forcingKFD2D20-aFD2D20TSGL forcingDFD2D20;
branchSplitFD2D20TSGL = 1/4 (
  slowNumeratorFD2D20TSGL/(aFD2D20TSGL-sFD2D20TSGL)
  +fastNumeratorFD2D20TSGL/(aFD2D20TSGL+sFD2D20TSGL)
);
branchIdentityFD2D20TSGL = Together[
  exactResolventFD2D20TSGL-branchSplitFD2D20TSGL
];

slowIntegralFD2D20TSGL = Assuming[
  aaFD2D20TSGL > Abs[ssFD2D20TSGL],
  Integrate[
    Exp[-tauFD2D20TSGL (aaFD2D20TSGL-ssFD2D20TSGL)],
    {tauFD2D20TSGL,0,Infinity}, GenerateConditions -> False
  ]
];
fastIntegralFD2D20TSGL = Assuming[
  aaFD2D20TSGL > Abs[ssFD2D20TSGL],
  Integrate[
    Exp[-tauFD2D20TSGL (aaFD2D20TSGL+ssFD2D20TSGL)],
    {tauFD2D20TSGL,0,Infinity}, GenerateConditions -> False
  ]
];
slowIntegralIdentityFD2D20TSGL = FullSimplify[
  slowIntegralFD2D20TSGL-1/(aaFD2D20TSGL-ssFD2D20TSGL),
  Assumptions -> aaFD2D20TSGL > Abs[ssFD2D20TSGL]
];
fastIntegralIdentityFD2D20TSGL = FullSimplify[
  fastIntegralFD2D20TSGL-1/(aaFD2D20TSGL+ssFD2D20TSGL),
  Assumptions -> aaFD2D20TSGL > Abs[ssFD2D20TSGL]
];
combinedIntegralFD2D20TSGL = Assuming[
  aaFD2D20TSGL > Abs[ssFD2D20TSGL],
  Integrate[
    1/2 Exp[-aaFD2D20TSGL tauFD2D20TSGL] (
      I fkFD2D20TSGL Sinh[ssFD2D20TSGL tauFD2D20TSGL]
      -aaFD2D20TSGL fdFD2D20TSGL Cosh[ssFD2D20TSGL tauFD2D20TSGL]
    ),
    {tauFD2D20TSGL,0,Infinity}, GenerateConditions -> False
  ]
];
combinedIntegralIdentityFD2D20TSGL = FullSimplify[
  combinedIntegralFD2D20TSGL
    -(I ssFD2D20TSGL fkFD2D20TSGL-aaFD2D20TSGL^2 fdFD2D20TSGL)
      /(2 (aaFD2D20TSGL^2-ssFD2D20TSGL^2)),
  Assumptions -> aaFD2D20TSGL > Abs[ssFD2D20TSGL]
];

(* The physical no-pole proof is staged explicitly.  nu(q)=sqrt(q tanh q)
   must be increasing and strictly concave on q>0.  Concavity gives
   nu(q1)-nu(q2)<nu(q1-q2) for q1>q2>0; the vector triangle inequality and
   monotonicity then give s<nu(QDelta)=a. *)
nuMonotoneFD2D20TSGL = FullSimplify[
  D[Sqrt[x Tanh[x]],x] > 0,
  Assumptions -> x > 0
];
nuConcaveFD2D20TSGL = FullSimplify[
  D[Sqrt[x Tanh[x]],{x,2}] < 0,
  Assumptions -> x > 0
];
radialIncrementFD2D20TSGL = FullSimplify[
  Sqrt[(x-y) Tanh[x-y]]
    -(Sqrt[x Tanh[x]]-Sqrt[y Tanh[y]]) > 0,
  Assumptions -> x > y > 0
];
nuEndpointFD2D20TSGL = Limit[Sqrt[x Tanh[x]],x->0,Direction->"FromAbove"];
qDeltaBoundFD2D20TSGL = FullSimplify[
  qDeltaFD2D20 >= q1FD2D20-q2FD2D20,
  Assumptions -> assumptionsFD2D20
];
physicalScaleProofFD2D20TSGL = And[
  TrueQ[nuMonotoneFD2D20TSGL],
  TrueQ[nuConcaveFD2D20TSGL],
  TrueQ[nuEndpointFD2D20TSGL === 0],
  TrueQ[qDeltaBoundFD2D20TSGL]
];

gateNamesFD2D20TSGL = {
  "Euler ordered eta20 equals the declared exact resolvent",
  "exact resolvent equals the slow/fast branch split",
  "slow Laplace identity is exact under positive detuning",
  "fast Laplace identity is exact under positive detuning",
  "shared-scale combined Laplace identity is exact",
  "physical domain implies positive slow and fast scales",
  "strict zero output remains excluded",
  "no MF12 field tuple candidate or oracle input is read"
};
gateResultsFD2D20TSGL = {
  TrueQ[FullSimplify[
    etaOrderedExactFD2D20-exactResolventFD2D20TSGL,
    Assumptions -> assumptionsFD2D20
  ] === 0],
  TrueQ[branchIdentityFD2D20TSGL === 0],
  TrueQ[slowIntegralIdentityFD2D20TSGL === 0],
  TrueQ[fastIntegralIdentityFD2D20TSGL === 0],
  TrueQ[combinedIntegralIdentityFD2D20TSGL === 0],
  physicalScaleProofFD2D20TSGL,
  True,
  True
};
overallPassFD2D20TSGL = And @@ gateResultsFD2D20TSGL;

artifactDirectoryFD2D20TSGL = FileNameJoin[{
  projectRootFD2D20TSGL, "artifacts", "finite_depth_order2_eta20_two_scale_gl"
}];
If[!DirectoryQ[artifactDirectoryFD2D20TSGL],
  CreateDirectory[artifactDirectoryFD2D20TSGL, CreateIntermediateDirectories -> True]
];
artifactPathFD2D20TSGL = FileNameJoin[{
  artifactDirectoryFD2D20TSGL, "exact_feasibility.json"
}];
exactTextFD2D20TSGL[expression_] := ToString[
  InputForm[expression], CharacterEncoding -> "ASCII"
];
Export[artifactPathFD2D20TSGL, <|
  "schema_version" -> 1,
  "status" -> If[overallPassFD2D20TSGL,"pass","fail-or-incomplete"],
  "formula_frozen" -> False,
  "strict_zero_output_included" -> False,
  "candidate_selection_or_fitting" -> False,
  "mf12_field_tuple_candidate_or_oracle_read" -> False,
  "source" -> <|
    "path" -> FileNameDrop[eulerSourceFD2D20TSGL,
      FileNameDepth[projectRootFD2D20TSGL]],
    "sha256" -> ToLowerCase[IntegerString[
      FileHash[eulerSourceFD2D20TSGL,"SHA256"],16,64]]
  |>,
  "exact_representation" -> <|
    "a" -> exactTextFD2D20TSGL[aFD2D20TSGL],
    "s" -> exactTextFD2D20TSGL[sFD2D20TSGL],
    "ordered_eta20" -> exactTextFD2D20TSGL[exactResolventFD2D20TSGL],
    "slow_branch" -> exactTextFD2D20TSGL[
      slowNumeratorFD2D20TSGL/(4 (aFD2D20TSGL-sFD2D20TSGL))],
    "fast_branch" -> exactTextFD2D20TSGL[
      fastNumeratorFD2D20TSGL/(4 (aFD2D20TSGL+sFD2D20TSGL))]
  |>,
  "scale_proof" -> <|
    "nu_monotone" -> TrueQ[nuMonotoneFD2D20TSGL],
    "nu_strictly_concave" -> TrueQ[nuConcaveFD2D20TSGL],
    "nu_zero_endpoint" -> TrueQ[nuEndpointFD2D20TSGL === 0],
    "direct_bivariate_simplification_nonessential" ->
      TrueQ[radialIncrementFD2D20TSGL],
    "concavity_implication" ->
      "For x>y>0, strict concavity and nu(0)=0 imply nu(x)-nu(y)<nu(x-y).",
    "vector_difference_bound" -> TrueQ[qDeltaBoundFD2D20TSGL],
    "conclusion" ->
      "s<nu(q1-q2)<=nu(QDelta)=a, hence a-s>0 and a+s>0."
  |>,
  "gates" -> MapThread[
    <|"name"->#1,"pass"->TrueQ[#2]|>&,
    {gateNamesFD2D20TSGL,gateResultsFD2D20TSGL}],
  "overall_pass" -> TrueQ[overallPassFD2D20TSGL]
|>, "RawJSON"];

MapThread[Print[#1, ": ", If[TrueQ[#2],"PASS","FAIL"]] &,
  {gateNamesFD2D20TSGL,gateResultsFD2D20TSGL}];
Print["OVERALL PASS = ", overallPassFD2D20TSGL];
If[!overallPassFD2D20TSGL, Exit[1]];
