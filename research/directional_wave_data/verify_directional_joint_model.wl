ClearAll[b,b1,b2,b3,x1,x2,y1,y2,theta,aa,bb,ww,tt,gg];
tot=b1+b2+b3;
checks=<|
 "observed_coefficient_recovered"->TrueQ[FullSimplify[Total[b {b1,b2,b3}/tot]-b,tot!=0]==0],
 "prior_identity"->TrueQ[FullSimplify[tot {b1,b2,b3}/tot-{b1,b2,b3},tot!=0]=={0,0,0}],
 "vector_dot_rotation_invariance"->TrueQ[FullSimplify[(x1 Cos[theta]-y1 Sin[theta])(x2 Cos[theta]-y2 Sin[theta])+(x1 Sin[theta]+y1 Cos[theta])(x2 Sin[theta]+y2 Cos[theta])-x1 x2-y1 y2]==0],
 "temporal_conjugation_eta"->TrueQ[FullSimplify[ComplexExpand[Re[(aa+I bb)Exp[I ww tt]]-Re[(aa-I bb)Exp[-I ww tt]]]]==0],
 "temporal_conjugation_psi"->TrueQ[FullSimplify[ComplexExpand[Re[I gg/ww (aa+I bb)Exp[I ww tt]]-Re[-I gg/ww (aa-I bb)Exp[-I ww tt]]]]==0]
|>;
root=DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];
out=FileNameJoin[{root,"artifacts","directional_joint_input"}];If[!DirectoryQ[out],CreateDirectory[out,CreateIntermediateDirectories->True]];
Export[FileNameJoin[{out,"algebra_checks.json"}],checks,"RawJSON"];Print[checks];If[And@@Values[checks],Exit[0],Exit[1]];
