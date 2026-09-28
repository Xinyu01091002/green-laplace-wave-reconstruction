ClearAll[t,tau,f,g,s,a,sd,sk,c];
shift[expr_]:=expr/.t->t-I tau;
old=Exp[-s tau](-a Sinh[a tau] sd+Cosh[a tau]sk);
transferred=((-a sd+sk)Exp[-(s-a)tau]+(a sd+sk)Exp[-(s+a)tau])/2;
sigma=Sqrt[7]/3;
checks=<|
 "analytic_time_shift_preserves_products"->TrueQ[Expand[shift[f[t]g[t]]-shift[f[t]]shift[g[t]]]==0],
 "output_node_is_original_GL_node"->TrueQ[FullSimplify[TrigToExp[old]-transferred]==0],
 "inner_potential_node_transfer"->TrueQ[FullSimplify[TrigToExp[I Exp[-s tau](Cosh[a tau]sd-Sinh[a tau]sk/a)]-I((sd-sk/a)Exp[-(s-a)tau]+(sd+sk/a)Exp[-(s+a)tau])/2]==0],
 "outer_cubic_node_transfer"->TrueQ[FullSimplify[TrigToExp[Exp[-s tau](a Sinh[a tau]sd-I Cosh[a tau]sk)]-((a sd-I sk)Exp[-(s-a)tau]+(-a sd-I sk)Exp[-(s+a)tau])/2]==0],
 "balance_shift_cancels"->TrueQ[FullSimplify[Exp[-tau(s-2 c)]Exp[-2 c tau]-Exp[-s tau]]==0],
 "commensurate_native_space_time_fixture"->TrueQ[FullSimplify[4/(1+sigma^2)-9/4]==0]
|>;
Print[checks];root=DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];
Export[FileNameJoin[{root,"artifacts","unidirectional_time_series","gl_time_node_transfer.json"}],checks,"RawJSON"];
If[And@@Values[checks],Exit[0],Exit[1]];
