(* Algebraic transfer of the existing GL kernel, not a new physical model. *)
ClearAll[aa, ss, tau, xx, sd, sk, u, v, p, r];
old = Exp[xx-ss tau] (-aa sd Sinh[aa tau]+sk Cosh[aa tau]);
stable = ((-aa sd+sk) Exp[xx-(ss-aa)tau]+(aa sd+sk) Exp[xx-(ss+aa)tau])/2;
checks = <|
 "stable_exponential_identity" -> TrueQ[FullSimplify[TrigToExp[old]-stable]==0],
 "dynamic_source_symmetrization" -> TrueQ[Expand[(2 v^2+u v-p r/(u v)+2 u^2+u v-p r/(u v))/2-(u^2+u v+v^2-p r/(u v))]==0],
 "kinematic_source_symmetrization" -> TrueQ[Expand[(2 r^2/v+2 p r/u+2 p^2/u+2 p r/v)/2-((p^2+p r)/u+(r^2+p r)/v)]==0],
 "ordered_pair_symmetry" -> TrueQ[FullSimplify[(u^2+u v+v^2-p r/(u v))-(v^2+v u+u^2-r p/(v u))]==0]
|>;
Print[checks];
root=DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];
out=FileNameJoin[{root,"artifacts","unidirectional_time_series"}];
If[!DirectoryQ[out],CreateDirectory[out,CreateIntermediateDirectories->True]];
Export[FileNameJoin[{out,"pair_algebra.json"}],checks,"RawJSON"];
If[And@@Values[checks],Exit[0],Exit[1]];
