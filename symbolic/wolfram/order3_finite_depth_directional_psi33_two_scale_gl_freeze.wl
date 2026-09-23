(* ::Package:: *)
(* Freeze the surface-trace relation for the independent finite-depth
   directional Psi33 Two-Scale/Shared-Scale graph.  No historical Psi33
   candidate, MF12 field/kernel, tuple, waveform, or oracle is imported. *)

ClearAll[
  epsPSG, zPSG, eta1PSG, eta2PSG, eta3PSG,
  phi10PSG, phi1zPSG, phi1zzPSG, phi20PSG, phi2zPSG, phi30PSG,
  etaSeriesPSG, bulkTaylorPSG, psi3ExtractedPSG, psi3ExpectedPSG,
  omegaSumPSG, omegaOutputPSG, resolventPSG, splitPSG,
  phiNumeratorPSG, gauge3PSG, gatesPSG, gatePSG, overallPSG,
  tauPSG, balancePSG, projectionSumPSG, balancedInputExponentPSG,
  projectRootPSG, etaFreezePathPSG, exportPathPSG, sourceFilesPSG,
  sourceHashesPSG, etaFreezeHashPSG
];

etaSeriesPSG = epsPSG eta1PSG + epsPSG^2 eta2PSG + epsPSG^3 eta3PSG;
bulkTaylorPSG =
  epsPSG (phi10PSG + zPSG phi1zPSG + zPSG^2 phi1zzPSG/2) +
  epsPSG^2 (phi20PSG + zPSG phi2zPSG) + epsPSG^3 phi30PSG;
psi3ExtractedPSG = Expand[
  SeriesCoefficient[bulkTaylorPSG /. zPSG -> etaSeriesPSG,
    {epsPSG, 0, 3}]
];
psi3ExpectedPSG = phi30PSG + eta1PSG phi2zPSG +
  eta2PSG phi1zPSG + eta1PSG^2 phi1zzPSG/2;

omegaSumPSG = Symbol["Omega"];
omegaOutputPSG = Symbol["omegaK"];
resolventPSG = 1/(omegaSumPSG^2 - omegaOutputPSG^2);
splitPSG = 1/(2 omegaOutputPSG) (
  1/(omegaSumPSG - omegaOutputPSG) -
  1/(omegaSumPSG + omegaOutputPSG));
phiNumeratorPSG = I Symbol["dTauFD"] - Symbol["FK"];
gauge3PSG = 0;
tauPSG = Symbol["tau"];
balancePSG = Symbol["cBalance"];
projectionSumPSG = Symbol["KParallel"];
balancedInputExponentPSG = -tauPSG (omegaSumPSG-balancePSG projectionSumPSG);

gatePSG[name_, value_] := Module[{pass = TrueQ[value]},
  Print[name <> " = " <> If[pass, "PASS", "FAIL"]];
  pass
];

gatesPSG = {
  gatePSG["surface-trace Taylor coefficient through order three",
    Expand[psi3ExtractedPSG - psi3ExpectedPSG] === 0],
  gatePSG["exact slow-fast resolvent factorization",
    Together[resolventPSG - splitPSG] === 0],
  gatePSG["flat Phi33 outer numerator convention",
    phiNumeratorPSG === I Symbol["dTauFD"] - Symbol["FK"]],
  gatePSG["no additive Psi33 gauge", gauge3PSG === 0],
  gatePSG["slow additive-projection exponential balance",
    Expand[balancedInputExponentPSG+
      tauPSG (omegaOutputPSG-balancePSG projectionSumPSG)+
      tauPSG (omegaSumPSG-omegaOutputPSG)] === 0],
  gatePSG["fast additive-projection exponential balance",
    Expand[balancedInputExponentPSG-
      tauPSG (omegaOutputPSG+balancePSG projectionSumPSG)+
      tauPSG (omegaSumPSG+omegaOutputPSG)] === 0],
  gatePSG["analytic-signal cubic normalization",
    2 (1/2)^3 == 1/4 && 6/24 == 1/4]
};
overallPSG = And @@ gatesPSG;

projectRootPSG = DirectoryName[$InputFileName, 3];
etaFreezePathPSG = FileNameJoin[{
  projectRootPSG, "symbolic", "generated",
  "finite_depth_directional_order3_two_scale_gl.json"
}];
If[!FileExistsQ[etaFreezePathPSG],
  Print["missing frozen eta33 resolvent record: ", etaFreezePathPSG];
  Exit[1]
];
etaFreezeHashPSG = IntegerString[FileHash[etaFreezePathPSG, "SHA256"], 16, 64];
sourceFilesPSG = {
  FileNameJoin[{projectRootPSG, "symbolic", "residuals",
    "order1_to_order3_symbolic_expansion_2d.wl"}],
  FileNameJoin[{projectRootPSG, "symbolic", "residuals",
    "order2_finite_depth_directional_euler.wl"}],
  FileNameJoin[{projectRootPSG, "symbolic", "residuals",
    "order3_finite_depth_directional_euler.wl"}]
};
sourceHashesPSG = AssociationThread[
  FileNameTake /@ sourceFilesPSG,
  IntegerString[FileHash[#, "SHA256"], 16, 64] & /@ sourceFilesPSG
];

exportPathPSG = FileNameJoin[{
  projectRootPSG, "symbolic", "generated",
  "finite_depth_directional_order3_psi33_two_scale_gl.json"
}];
Export[exportPathPSG, <|
  "schema_version" -> 1,
  "status" -> "frozen_pre_validation",
  "candidate" -> "finite_depth_directional_psi33_two_scale_gl",
  "scope" -> "stationary strict-forward positive-pure-sum Psi33 q_parent>0.5",
  "candidate_input_fields" -> {"eta11"},
  "target" -> "Psi33=phi(x,y,z=eta,t)|order3",
  "flat_phi33_is_internal_only" -> True,
  "external_lower_order_fields" -> False,
  "oracle_imported" -> False,
  "mf12_imported" -> False,
  "surface_trace_relation" ->
    "Psi33=Phi33+eta11 phi22_z+eta22 phi11_z+(1/2) eta11^2 phi11_zz",
  "outer_phi_numerator" -> "i d_tau F_D-F_K",
  "exponential_balance" ->
    "c=(1/2) min_active(omega/k_parallel); input exp[-tau(omega-c k_parallel)]; exact output compensation uses sum k_parallel=K_parallel",
  "gauge" -> <|"C3(t)" -> "0", "strict_K0" -> "0"|>,
  "analytic_signal_normalization" -> <|
    "input" -> "U=2 eta11_positive",
    "internal_input_factor" -> "1/2",
    "output_factor" -> "2",
    "net_cubic_factor" -> "1/4",
    "ordered_permutations_over_grouped_factor" -> "6/24=1/4"
  |>,
  "reused_eta33_resolvent_json" -> etaFreezePathPSG,
  "reused_eta33_resolvent_sha256" -> etaFreezeHashPSG,
  "source_sha256" -> sourceHashesPSG,
  "claim_boundary" ->
    "surface-trace identity and frozen finite-node GL graph only; no uniform shallow-water or all-angle claim"
|>, "RawJSON"];

Print["Psi33 Taylor coefficient = ", InputForm[psi3ExtractedPSG]];
Print["eta33 resolvent sha256 = ", etaFreezeHashPSG];
Print["artifact = ", exportPathPSG];
Print["OVERALL_ORDER3_PSI33_TWO_SCALE_GL_FREEZE = ",
  If[overallPSG, "PASS", "FAIL"]];
If[!overallPSG, Exit[1]];
