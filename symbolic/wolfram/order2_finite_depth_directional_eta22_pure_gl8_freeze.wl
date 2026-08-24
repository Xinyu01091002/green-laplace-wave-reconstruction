(* ::Package:: *)
(* Frozen pure Green--Laplace rank-8 eta22 representation.
   The candidate contains no Stokes-diagonal or angular correction. *)

residualSourceEta22GL8 = FileNameJoin[{
  DirectoryName[$InputFileName], "..", "residuals",
  "order2_finite_depth_directional_euler.wl"
}];
Get[residualSourceEta22GL8];

ClearAll[
  aaEta22GL8, ssEta22GL8, aEta22GL8, sEta22GL8,
  exactSinhEta22GL8, exactCoshEta22GL8,
  numeratorEta22GL8, exactKernelEta22GL8, continuousKernelEta22GL8,
  jacobiEta22GL8, characteristicEta22GL8, nodesEta22GL8,
  weightsEta22GL8, numericSystemEta22GL8, numericPairsEta22GL8,
  numericNodesEta22GL8, numericWeightsEta22GL8, momentsEta22GL8,
  momentDegreeEta22GL8, momentPassEta22GL8,
  exactTextEta22GL8, gateEta22GL8, gatesEta22GL8, rootEta22GL8,
  artifactEta22GL8
];

aEta22GL8 = Sqrt[dnoFD2[qOutputFD2]];
sEta22GL8 = nuSumFD2;
numeratorEta22GL8 =
  (dnoFD2[qOutputFD2] sourceDFD2-sEta22GL8 sourceKFD2)/4;
exactKernelEta22GL8 = Together[
  numeratorEta22GL8/(aEta22GL8^2-sEta22GL8^2)];

exactSinhEta22GL8 = FullSimplify[
  Integrate[Exp[-ssEta22GL8 tEta22GL8]
    Sinh[aaEta22GL8 tEta22GL8]/aaEta22GL8,
    {tEta22GL8,0,Infinity},Assumptions->ssEta22GL8>aaEta22GL8>0],
  Assumptions->ssEta22GL8>aaEta22GL8>0];
exactCoshEta22GL8 = FullSimplify[
  Integrate[Exp[-ssEta22GL8 tEta22GL8]
    Cosh[aaEta22GL8 tEta22GL8],
    {tEta22GL8,0,Infinity},Assumptions->ssEta22GL8>aaEta22GL8>0],
  Assumptions->ssEta22GL8>aaEta22GL8>0];
continuousKernelEta22GL8 = Together[
  (-aEta22GL8^2 sourceDFD2
      (exactSinhEta22GL8/.{aaEta22GL8->aEta22GL8,
        ssEta22GL8->sEta22GL8})
    +sourceKFD2 (exactCoshEta22GL8/.{aaEta22GL8->aEta22GL8,
        ssEta22GL8->sEta22GL8}))/4];

jacobiEta22GL8 = Normal[SparseArray[{
  Band[{1,1}]->Range[1,15,2],
  Band[{1,2}]->Range[7],
  Band[{2,1}]->Range[7]
},{8,8}]];
characteristicEta22GL8 = CharacteristicPolynomial[
  jacobiEta22GL8,xEta22GL8];
nodesEta22GL8 = Table[
  Root[Function[{xEta22GL8},Evaluate[characteristicEta22GL8]],index],
  {index,1,8}];
weightsEta22GL8 = (#/(9^2 LaguerreL[9,#]^2)& /@ nodesEta22GL8);
momentsEta22GL8 = (MatrixPower[jacobiEta22GL8,#][[1,1]]&) /@
  Range[0,15];
momentPassEta22GL8 = SameQ[
  momentsEta22GL8,Factorial /@ Range[0,15]];
Print["GL8 exact moment residuals = ",InputForm[
  momentsEta22GL8-(Factorial /@ Range[0,15])]];
numericSystemEta22GL8 = Eigensystem[N[jacobiEta22GL8,40]];
numericPairsEta22GL8 = SortBy[
  Transpose[{numericSystemEta22GL8[[1]],
    numericSystemEta22GL8[[2,All,1]]^2}],First];
numericNodesEta22GL8 = numericPairsEta22GL8[[All,1]];
numericWeightsEta22GL8 = numericPairsEta22GL8[[All,2]];

exactTextEta22GL8[value_] := ToString[InputForm[value]];
gateEta22GL8[name_,residual_] := Module[{pass=TrueQ[PossibleZeroQ[residual]]},
  Print[name," = ",If[pass,"PASS","FAIL"]];pass];
gatesEta22GL8 = {
  gateEta22GL8["sinh Green integral",
    Together[exactSinhEta22GL8
      -1/(ssEta22GL8^2-aaEta22GL8^2)]],
  gateEta22GL8["cosh Green integral",
    Together[exactCoshEta22GL8
      -ssEta22GL8/(ssEta22GL8^2-aaEta22GL8^2)]],
  gateEta22GL8["continuous Green kernel equals exact Euler kernel",
    Together[continuousKernelEta22GL8-exactKernelEta22GL8]],
  momentPassEta22GL8
};
Print["Gauss-Laguerre degree-fifteen exactness = ",
  If[momentPassEta22GL8,"PASS","FAIL"]];

rootEta22GL8 = DirectoryName[$InputFileName,3];
artifactEta22GL8 = FileNameJoin[{rootEta22GL8,"symbolic","generated",
  "finite_depth_directional_order2_eta22_pure_gl8.json"}];
Export[artifactEta22GL8,<|
  "schema_version"->1,
  "status"->"frozen-pure-gl8-retained-bounded",
  "candidate_id"->"gl-eta22-pure-gl8-detroot-v1",
  "residual_source"->"symbolic/residuals/order2_finite_depth_directional_euler.wl",
  "residual_source_sha256"->FileHash[
    residualSourceEta22GL8,"SHA256","HexString"],
  "scope"->"finite-depth strict-forward positive-pure-sum eta22; certified Q/2>0.3 validation domain",
  "candidate_input_fields"->{"eta11"},
  "external_lower_order_fields"->False,
  "oracle_or_mf12_used_in_formula"->False,
  "sampled_coefficients_used"->False,
  "stokes_correction_used"->False,
  "angular_correction_used"->False,
  "quadrature_family"->"standard Gauss-Laguerre determinant-root scaled",
  "quadrature_rank"->8,
  "quadrature_exact_moment_degree"->15,
  "nodes_exact"->(exactTextEta22GL8 /@ nodesEta22GL8),
  "weights_exact"->(exactTextEta22GL8 /@ weightsEta22GL8),
  "nodes"->N[numericNodesEta22GL8,17],
  "weights"->N[numericWeightsEta22GL8,17],
  "scale"->"lambda(qp)=Sqrt[(2 nu(qp))^2-2 qp Tanh[2 qp]]",
  "continuous_identity"->"1/(s^2-A)=Integral Exp[-s t] Sinh[Sqrt[A] t]/Sqrt[A] dt",
  "dynamic_kinematic_split"->"(-A SD I_sinh + SK I_cosh)/4",
  "stable_backend"->"common exponential balance with exact two-parent compensation",
  "cost"-><|"fft_ifft"->85,"pointwise_products"->57,"pair_loops"->0|>,
  "promotion_provenance"->"rank retained after the already completed post-freeze non-selecting GL2/4/6/8 validation ladder",
  "claim_boundary"->"bounded rank-8 pure-GL result; no uniform shallow, near-opposition, exact-zero, or all-spectrum claim",
  "gates"->gatesEta22GL8,
  "overall_exact_gate_pass"->And@@gatesEta22GL8
|>,"RawJSON"];
Print["artifact = ",artifactEta22GL8];
Print["OVERALL_ORDER2_FINITE_DEPTH_DIRECTIONAL_ETA22_PURE_GL8_FREEZE = ",
  If[And@@gatesEta22GL8,"PASS","FAIL"]];
If[!And@@gatesEta22GL8,Exit[1]];
