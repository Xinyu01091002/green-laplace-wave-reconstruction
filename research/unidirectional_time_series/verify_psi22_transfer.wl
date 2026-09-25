ClearAll[s,a,t,sd,sk,u,v];
checks=<|
 "dual_branch_identity"->TrueQ[FullSimplify[TrigToExp[I/4 Exp[-s t](Cosh[a t]sd-Sinh[a t]sk/a)]-
 I/8 ((sd-sk/a)Exp[-(s-a)t]+(sd+sk/a)Exp[-(s+a)t])]==0],
 "surface_Taylor_symmetrization"->TrueQ[Expand[(-I v/2-I u/2)/2+I(u+v)/4]==0]
|>;
root=DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];
Export[FileNameJoin[{root,"artifacts","unidirectional_time_series","psi22_algebra.json"}],checks,"RawJSON"];
Print[checks];If[And@@Values[checks],Exit[0],Exit[1]];
