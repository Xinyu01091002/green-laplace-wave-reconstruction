(* Exact identities used by the direct joint-modal eta22 prototype. *)
ClearAll[s,a,tau,nu1,nu2,sd,sk];
assumptions=s>a>0&&nu1>0&&nu2>0;
gates={
  FullSimplify[
    Integrate[Exp[-s tau] Sinh[a tau]/a,{tau,0,Infinity},
      Assumptions->s>a>0]-1/(s^2-a^2),
    Assumptions->s>a>0]===0,
  FullSimplify[
    Integrate[Exp[-s tau] Cosh[a tau],{tau,0,Infinity},
      Assumptions->s>a>0]-s/(s^2-a^2),
    Assumptions->s>a>0]===0,
  FullSimplify[
    Exp[-nu1 tau] Exp[-nu2 tau]-Exp[-(nu1+nu2) tau],
    Assumptions->nu1>0&&nu2>0]===0,
  FullSimplify[
    (-a^2 sd/(s^2-a^2)+s sk/(s^2-a^2))/4-
      (-a^2 sd+s sk)/(4 (s^2-a^2)),
    Assumptions->s>a>0]===0
};
Print["DIRECT_MODAL_GL22_GATES = ",gates];
Print["DIRECT_MODAL_GL22_OVERALL = ",If[And@@gates,"PASS","FAIL"]];
If[!And@@gates,Exit[1]];
