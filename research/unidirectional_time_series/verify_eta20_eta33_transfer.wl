ClearAll[u,v,p,r,sd,sk,g,pp,ee,s,a,x,t,lambda];
sd=u^2+u v+v^2-p r/(u v);
sk=(p^2+p r)/u+(r^2+p r)/v;
ds=(2 v^3+4 u v^2-2 p r/v+2 u^3+4 v u^2-2 p r/u)/2;
ks=(2 (u r^2/v+r^2)+2(p r+p v r/u)+2(v p^2/u+p^2)+2(r p+r u p/v))/2;
(* Fold pair-times-single and single-times-pair into the same ordered sum. *)
fkOriginal=(p/u) I r ee+(I r pp) I p- r^2 pp-ee*(-I p^2/u);
fkFolded=I (p r+p^2)/u ee-(p r+r^2)pp;
fdOriginal=g*(-I s)*pp-ee*u^2+(p/u)*I r pp+(-I u)*g*pp;
fdFolded=-I s g pp-u^2 ee+I p r/u pp-I u g pp;
checks=<|
 "pair_dynamic_time_derivative"->TrueQ[FullSimplify[ds-(u+v)sd]==0],
 "pair_kinematic_time_derivative"->TrueQ[FullSimplify[ks-(u+v)sk]==0],
 "cubic_pair_kinematic_fold"->TrueQ[Expand[fkOriginal-fkFolded]==0],
 "cubic_pair_dynamic_fold"->TrueQ[Expand[fdOriginal-fdFolded]==0],
 "difference_branch_balance"->TrueQ[FullSimplify[Exp[t*(u-v-a)]+Exp[-t*(u-v+a)]-2 Exp[-a t] Cosh[(u-v)t]]==0]
|>;
root=DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];
Export[FileNameJoin[{root,"artifacts","unidirectional_time_series","eta20_eta33_algebra.json"}],checks,"RawJSON"];
Print[checks];If[And@@Values[checks],Exit[0],Exit[1]];
