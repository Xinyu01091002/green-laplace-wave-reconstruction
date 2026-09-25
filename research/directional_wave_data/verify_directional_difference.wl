ClearAll[x1,y1,x2,y2,u,v,nu,fd,fk,t,theta];
dot=x1 x2+y1 y2;q1=x1^2+y1^2;q2=x2^2+y2^2;
checks=<|
 "difference_forcing_vector_form"->TrueQ[Expand[I(q1/u-dot/u+dot/v-q2/v)-I((q1-dot)/u+(dot-q2)/v)]==0],
 "opposite_pair_hermitian_kernel"->TrueQ[FullSimplify[((I fk-nu fd)Exp[-(nu-u+v)t]+(-I fk-nu fd)Exp[-(nu+u-v)t])-((I(-fk)-nu fd)Exp[-(nu-v+u)t]+(-I(-fk)-nu fd)Exp[-(nu+v-u)t])]==0],
 "equal_frequency_distinct_direction_nonzero_K"->TrueQ[FullSimplify[(Cos[theta]-Cos[-theta])^2+(Sin[theta]-Sin[-theta])^2-4 Sin[theta]^2]==0]
|>;
root=DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];out=FileNameJoin[{root,"artifacts","directional_sweep"}];
If[!DirectoryQ[out],CreateDirectory[out,CreateIntermediateDirectories->True]];
Export[FileNameJoin[{out,"eta20_algebra.json"}],checks,"RawJSON"];Print[checks];If[And@@Values[checks],Exit[0],Exit[1]];
