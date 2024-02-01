(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"];


(* ::Input::Initialization:: *)
$HistoryLength=1;


(* ::Input::Initialization:: *)
$MinPrecision=15;


(* ::Input::Initialization:: *)
$MaxExtraPrecision=50;


(* ::Section::Initialization::Closed:: *)
(*(*(*Einstein equations in Eddington-Finkelstein gauge*)*)*)


(* ::Input::Initialization:: *)
(* Coordinates and metric *)
xmu={r,t,y,x};
d=Length[xmu];
metric={{0,1,0,0},{1,-A[r,t],0,0},{0,0,E^- B[r,t] S[r,t]^2,0},{0,0,0,E^B[r,t] S[r,t]^2}};
metricinv=Inverse[metric];
g=Det[metric];
(*metric//MatrixForm*)


(* ::Input::Initialization:: *)
(* Christoffel symbol and contraction *)
Gudd=Table[Sum[1/2metricinv[[i,l]](D[metric[[k,l]],xmu[[j]]]+D[metric[[j,l]],xmu[[k]]]-D[metric[[j,k]],xmu[[l]]]),{l,1,d}],{i,1,d},{j,1,d},{k,1,d}];
GGuddd=Table[Sum[Gudd[[i,j,m]]*Gudd[[m,k,l]],{m,1,d}],{i,1,d},{j,1,d},{k,1,d},{l,1,d}];


(* ::Input::Initialization:: *)
(* Riemann tensor, Ricci tensor and Ricci scalar *)
Ruddd=Table[D[Gudd[[i,l,j]],xmu[[k]]]-D[Gudd[[i,k,j]],xmu[[l]]]+GGuddd[[i,k,l,j]]-GGuddd[[i,l,k,j]],{i,1,d},{j,1,d},{k,1,d},{l,1,d}];
Rdd=Table[Sum[Ruddd[[m,i,m,j]],{m,1,d}],{i,1,d},{j,1,d}]//Simplify;
Rud=Table[Sum[metricinv[[i,m]]*Rdd[[m,j]],{m,1,d}],{i,1,d},{j,1,d}]//Simplify;
Rscal=Sum[Rud[[m,m]],{m,1,d}]//Simplify;


(* ::Input::Initialization:: *)
(* Vacuum Einstein equations + cosmological constant *)
replL={L->1}; (* we set the AdS radius to 1 *)
Lambda=-(d-1)(d-2)/(2L^2)/.replL;
EinsteinEqns=Rdd-1/2metric*Rscal+Lambda*metric//Simplify;
EE1=Drop[Union[Flatten[EinsteinEqns]],1];
(* Depending on the Mathematica version the equations are ordered differently! *)
(*EE={EE1[[1]],EE1[[2]],EE1[[5]],EE1[[3]],EE1[[4]]};*)(* Mathematica<11 *)
EE={EE1[[1]],EE1[[3]],EE1[[4]],EE1[[2]],EE1[[5]]};(* Mathematica11 *)
(*TableForm[EEeq=#\[Equal]0&/@EE]*)


(* ::Input::Initialization:: *)
Put[EE,"ee.dat"];


(* ::Section::Initialization::Closed:: *)
(*(*(*Coordinate transformation r-> z=1/r *)*)*)


(* ::Input::Initialization:: *)
Dr[func_, times_Integer] := D[func[1/r, t], {r, times}] /. r -> 1/z;
replrtoz = {\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{
RowBox[{"Derivative", "[", 
RowBox[{"ddr_", ",", "ddt_"}], "]"}], "[", "f_", "]"}], "[", 
RowBox[{"r", ",", "t"}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) -> Dr[\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{"Derivative", "[", 
RowBox[{"0", ",", "ddt"}], "]"}], "[", "f", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\), ddr], r -> z};
EOMz=EE/.replrtoz ;


(* ::Input::Initialization:: *)
Put[EOMz,"eez.dat"];


(* ::Section::Initialization::Closed:: *)
(*(*(*Near boundary analysis*)*)*)


(* ::Input::Initialization:: *)
(* ansatz for the near boundary expansion, for order>7 the series ansatz must be extended with higher power log terms. *)
Asb[z_,t_,ord_]:=1/(z \[Epsilon])^2Sum[(Subscript[a, n][t]+Subscript[\[Alpha]1, n][t]Log[z]+Subscript[\[Alpha]2, n][t]Log[z]^2)(z \[Epsilon])^n,{n,0,ord}];
Ssb[z_,t_,ord_]:=1/(z \[Epsilon]) Sum[(Subscript[s, n][t]+Subscript[\[Sigma]1, n][t]Log[z]+Subscript[\[Sigma]2, n][t]Log[z]^2)(z \[Epsilon])^n,{n,0,ord}];
Bsb[z_,t_,ord_]:=Sum[(Subscript[b, n][t]+Subscript[\[Beta]1, n][t]Log[z]+Subscript[\[Beta]2, n][t]Log[z]^2)(z \[Epsilon])^n,{n,0,ord}];


(* ::Input::Initialization:: *)
(* boundary conditions,Subscript[a, 1]\[Rule](0&) *)
BCs={Subscript[a, 0][t]->1,Subscript[a, 1]->(0&),Subscript[s, 0][t]->1,Subscript[b, 0][t]->b0[t]};


(* ::Input::Initialization:: *)
(* replacements for putting the expansion parameter into the EOMs *)
ztoZ={f_[z,t]->f[Z,t]};
putEps={\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"n_", ",", "m_"}], ")"}],
Derivative],
MultilineFunction->None]\)[Z,t]->\[Epsilon]^-n \!\(\*SuperscriptBox[\(f\), 
TagBox[
RowBox[{"(", 
RowBox[{"n", ",", "m"}], ")"}],
Derivative],
MultilineFunction->None]\)[Z,t],z->\[Epsilon] z};
Ztoz={f_[Z,t]->f[z,t]};


(* ::Input::Initialization:: *)
(* overall factors make the EOMs have the same leading order in z *)
EOMzEps={1/z^2 EOMz[[1]],z^2EOMz[[2]],z^2EOMz[[3]],z^2EOMz[[4]],z^2EOMz[[5]]}/.ztoZ/.putEps/.Ztoz;


(* ::Input::Initialization:: *)
(* solving the Einstein equations order by order in z and save the result to a file *)
order=5;resFile="NBsolAdS4Ord5.res";datFile="NBsolAdS4Ord5.dat";
sol=BCs;putDel={Log[z]->\[Delta]};
Do[
Ab[z_,t_]=Asb[z,t,i]/.sol;
Sb[z_,t_]=Ssb[z,t,i]/.sol;
Bb[z_,t_]=Bsb[z,t,i]/.sol;
replBdry={A->Ab,S->Sb,B->Bb};
EOMnow=SeriesCoefficient[EOMzEps[[;;]]/.replBdry/.putDel,{\[Epsilon],0,i},{\[Delta],0,0}];
EOMnowLog=SeriesCoefficient[EOMzEps[[;;]]/.replBdry/.putDel,{\[Epsilon],0,i},{\[Delta],0,1}];
EOMnowLog2=SeriesCoefficient[EOMzEps[[;;]]/.replBdry/.putDel,{\[Epsilon],0,i},{\[Delta],0,2}];
solnow=Solve[Thread[{EOMnow[[;;]]==0,EOMnowLog[[;;]]==0,EOMnowLog2[[;;]]==0}],{Subscript[a, 3]'[t],Subscript[b, 3]'[t],Subscript[a, i][t],Subscript[b, i][t],Subscript[s, i][t],Subscript[\[Alpha]1, i][t],Subscript[\[Sigma]1, i][t],Subscript[\[Beta]1, i][t],Subscript[\[Alpha]2, i][t],Subscript[\[Sigma]2, i][t],Subscript[\[Beta]2, i][t]}]//Flatten//Quiet;
sol=Flatten[{sol,solnow}];
checknow=Thread[SeriesCoefficient[EOMzEps[[;;]]/.replBdry//.sol,{\[Epsilon],0,i}]==0]//Simplify;
,{i,0,order}];
nbSol=sol;
Clear[order];
Put[nbSol//Simplify,resFile];
Put[nbSol//Simplify,datFile];


(* ::Section::Initialization::Closed:: *)
(*(*(*Transformation from EF  to FG coordinates*)*)*)


(* ::Input::Initialization:: *)
(* metric in EF coordinates *)
metricEFz=({
 {0, -1/z^2, 0, 0},
 {-1/z^2, -A[z,t], 0, 0},
 {0, 0, S[z,t]^2 E^-B[z,t], 0},
 {0, 0, 0, S[z,t]^2 E^B[z,t]}
});


(* ::Input::Initialization:: *)
(* loading the near-boundary solution *)
(*SetDirectory[NotebookDirectory[]];*)
nbSol=Get["NBsolAdS4Ord5.res"];order=5;
Clear[Ab,Sb,Bb];


(* ::Input::Initialization:: *)
(* near-boundary series in EF coordinates *)
Ab[z_,t_]=Asb[z,t,order]/.nbSol/.\[Epsilon]->1;
Sb[z_,t_]=Ssb[z,t,order]/.nbSol/.\[Epsilon]->1;
Bb[z_,t_]=Bsb[z,t,order]/.nbSol/.\[Epsilon]->1;


(* ::Input::Initialization:: *)
(* line element in EF-coordinates; we need to transform only the dudv and dv^2 parts, the rest is already diagonal *)
ds2EF=-(2/z^2)dz dt-Ab[z,t]dt^2;


(* ::Input::Initialization:: *)
(* series ansatz for the EF-coords in terms of FG-coords *)
ZEFs[v_,u_,ord_]:=Sum[(Subscript[zEF, n][v]+Subscript[\[Zeta]EF, n][v]Log[u])u^n,{n,1,ord}];
TEFs[v_,u_,ord_]:=v+Sum[(Subscript[tEF, n][v]+Subscript[\[Tau]EF, n][v]Log[u])u^n,{n,1,ord}];


(* ::Input::Initialization:: *)
BCs={Subscript[tEF, 1][v]->-1,Subscript[\[Tau]EF, 1][v]->0};(* this follows from Subscript[g, tt]=-1 *)
sol=BCs;
Do[
ZEF=ZEFs[v,u,i]/.sol;
TEF=TEFs[v,u,i]/.sol;
replFG={z->ZEF,t->TEF,dz->D[ZEF,u] du+D[ZEF,v] dv,dt->D[TEF,u] du+D[TEF,v] dv};
(* derivative expansion of the line element *)
ds2FG=ds2EF//.replFG;
ds2FGnow=Coefficient[Series[ds2FG,{u,0,i-3}]/.Log[u]->0,{du^2,du dv}];
ds2FGnowLog=Coefficient[D[Series[ds2FG,{u,0,i-3}],Log[u]],{du^2,du dv}];
solNow=Solve[{ds2FGnow[[1]]==1/u^2,ds2FGnow[[2]]==0,
ds2FGnowLog[[1]]==0,ds2FGnowLog[[2]]==0},{Subscript[tEF, i][v],Subscript[zEF, i][v],Subscript[\[Zeta]EF, i][v],Subscript[\[Tau]EF, i][v]}]//Quiet;
(*Print[solNow];*)
sol=Flatten[{sol,solNow}];
,{i,1,6}];


(* ::Input::Initialization:: *)
(* FG coordinates in terms of EF series coefficients *)
ZEF=ZEFs[v,u,6]/.sol;
TEF=TEFs[v,u,6]/.sol;
(* expansion of the v-coordinate *)
(*Series[TEF,{u,0,5}]*)
(* expansion of the u-coordinate *)
(*Series[ZEF,{u,0,5}]*)


(* ::Input::Initialization:: *)
(* metric in FG coordinates *)
replFG={z->ZEF,t->TEF,dz->D[ZEF,u] du+D[ZEF,v] dv,dt->D[TEF,u] du+D[TEF,v] dv};
ds2FG=ds2EF+Sb[z,t]^2(E^-Bb[z,t] dy^2+E^ Bb[z,t] dx^2)//.replFG;
(* leading order expansion of the metric in FG coordinates *)
(*Collect[Series[ds2FG,{u,0,-2}],{du,dv,dy,dx}];*)


(* ::Input::Initialization:: *)
(* line element in u-coordinate *)
ds2u=Collect[Series[ds2FG,{u,0,2}],{du,dv,dx,u,Log[u]}];


(* ::Input::Initialization:: *)
(* metric in \[Rho]-coordinate (ds^2=d\[Rho]^2/(4\[Rho])+1/\[Rho]*Subscript[g, ij](x)dx^idx^j) such as used in (arXiv:hep-th/0002230) *)
replrho={du^2->1/(4\[Rho])d\[Rho]^2,u->Sqrt[\[Rho]]};
g\[Rho]\[Rho]=Coefficient[Series[ds2u/.replrho,{\[Rho],0,-1}],d\[Rho]^2];
gvv0=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,-1}],dv^2];
gvv2=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,0}],dv^2];
gvv3=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,1/2}],dv^2];
gvv4=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,1}],dv^2]/.Log[\[Rho]]->0;
hvv4=D[Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,1}],dv^2],Log[\[Rho]]];
gyy0=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,-1}],dy^2];
gyy2=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,0}],dy^2];
gyy3=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,1/2}],dy^2];
gyy4=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,1}],dy^2]/.Log[\[Rho]]->0;
hyy4=D[Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,1}],dy^2],Log[\[Rho]]];
gxx0=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,-1}],dx^2];
gxx2=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,0}],dx^2];
gxx3=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,1/2}],dx^2];
gxx4=Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,1}],dx^2]/.Log[\[Rho]]->0;
hxx4=D[Coefficient[SeriesCoefficient[ds2u/.replrho,{\[Rho],0,1}],dx^2],Log[\[Rho]]];
g0=DiagonalMatrix[{gvv0,gyy0,gxx0}];
g0inv=Inverse[g0];
g2=DiagonalMatrix[{gvv2,gyy2,gxx2}];
g3=DiagonalMatrix[{gvv3,gyy3,gxx3}];
g4=DiagonalMatrix[{gvv4,gyy4,gxx4}];
h14=DiagonalMatrix[{hvv4,hyy4,hxx4,hxx4}];
G\[Rho]\[Rho]=g\[Rho]\[Rho];
Gvv=1/\[Rho](gvv0+\[Rho] gvv2+ \[Rho]^(3/2) gvv3+\[Rho]^2 gvv4+\[Rho]^2 Log[\[Rho]]hvv4);
Gyy=1/\[Rho](gyy0+\[Rho] gyy2+ \[Rho]^(3/2) gyy3+\[Rho]^2 gyy4+\[Rho]^2 Log[\[Rho]]hyy4);
Gxx=1/\[Rho](gxx0+\[Rho] gxx2+ \[Rho]^(3/2) gxx3+\[Rho]^2 gxx4+\[Rho]^2 Log[\[Rho]]hxx4);
Gdd=DiagonalMatrix[{G\[Rho]\[Rho],Gvv,Gyy,Gxx}];


(* ::Input::Initialization:: *)
(* A,B and S field in FG coordinates *)
AFG[\[Rho]_,v_]:=Collect[Series[Ab[z,t]/.replFG/.replrho,{\[Rho],0,3}],{\[Rho],Log[\[Rho]]}];
BFG[\[Rho]_,v_]:=Collect[Series[Bb[z,t]/.replFG/.replrho,{\[Rho],0,3}],{\[Rho],Log[\[Rho]]}];
SFG[\[Rho]_,v_]:=Collect[Series[Sb[z,t]/.replFG/.replrho,{\[Rho],0,3}],{\[Rho],Log[\[Rho]]}];


(* ::Input::Initialization:: *)
Put[{Simplify[AFG[\[Rho],v]],Simplify[BFG[\[Rho],v]],Simplify[SFG[\[Rho],v]]},"fg.dat"];


(* ::Section::Initialization::Closed:: *)
(*(*(*Holographic EMT*)*)*)


(* ::Input::Initialization:: *)
(*Print["Holographic EMT:"]*)
(* Holographic energy-momentum tensor from Eq.(3.18) of (0002230) *)
EMTdd=g3//Simplify;
(*Print["Subscript[T, \[Mu]\[Nu]][v]=",MatrixForm[EMTdd]]*)


(* ::Input::Initialization:: *)
(*Print["Trace of the EMT:"]*)
(* the EMT is traceless *)
EMTud=g0inv . EMTdd;
(*Print["Subsuperscript[T, \[Mu], \[Mu]][v]=",Tr[EMTud]//Expand]*)


(* ::Input::Initialization:: *)
(* Coordinates and boundary metric *)
xmuB={v,y,x};
dB=Length[xmuB];
metricB=g0;
metricBinv=Inverse[g0];


(* ::Input::Initialization:: *)
(* Christoffel symbol of the boundary metric *)
GuddB=Table[Sum[1/2metricBinv[[i,l]](D[metricB[[k,l]],xmuB[[j]]]+D[metricB[[j,l]],xmuB[[k]]]-D[metricB[[j,k]],xmuB[[l]]]),{l,1,dB}],{i,1,dB},{j,1,dB},{k,1,dB}];


(* ::Input::Initialization:: *)
(* the EMT is covariantly conserved *)
EMTuu=g0inv . EMTud;
(* relation obtained in the near boundary analysis *)
a3dv=Derivative[1][Subscript[a, 3]][t]/.nbSol/.t->v;
repla3dv={Derivative[1][Subscript[a, 3]][v]->a3dv};
(*Print["Covariant derivative of the EMT:"]*)
(* covariant derivative: Subscript[D, \[Mu]]T^\[Mu]\[Nu]=\!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]
\*SuperscriptBox[\(T\), \(\[Mu]\[Nu]\)]\)+\!\(
\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Sigma]\[Mu]\), \(\[Mu]\)]
\*SuperscriptBox[\(T\), \(\[Sigma]\[Nu]\)]\)+\!\(
\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Sigma]\[Mu]\), \(\[Nu]\)]
\*SuperscriptBox[\(T\), \(\[Mu]\[Sigma]\)]\) *)
(*Print["Subscript[D, \[Mu]]T^\[Mu]\[Nu]==Subscript[\[PartialD], \[Mu]]T^\[Mu]\[Nu]+Subsuperscript[\[CapitalGamma], \[Sigma]\[Mu], \[Mu]]T^\[Sigma]\[Nu]+Subsuperscript[\[CapitalGamma], \[Sigma]\[Mu], \[Nu]]T^\[Mu]\[Sigma]=",MatrixForm[Table[Sum[D[EMTuu[[i,j]],xmuB[[i]]],{i,1,dB}],{j,1,dB}]+
Table[Sum[GuddB[[k,i,j]]EMTuu[[i,j]],{i,1,dB},{j,1,dB}],{k,1,dB}]+
Table[Sum[GuddB[[i,j,i]]EMTuu[[j,k]],{i,1,dB},{j,1,dB}],{k,1,dB}]/.repla3dv//Simplify]]*)
ward=MatrixForm[Table[Sum[D[EMTuu[[i,j]],xmuB[[i]]],{i,1,dB}],{j,1,dB}]+
Table[Sum[GuddB[[k,i,j]]EMTuu[[i,j]],{i,1,dB},{j,1,dB}],{k,1,dB}]+
Table[Sum[GuddB[[i,j,i]]EMTuu[[j,k]],{i,1,dB},{j,1,dB}],{k,1,dB}]/.repla3dv//Simplify];
Export["ward.dat","\!\(\*SubscriptBox[\(D\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(T\), \(\[Mu]\[Nu]\)]\)==\!\(\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(T\), \(\[Mu]\[Nu]\)]\)+\!\(\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Sigma]\[Mu]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(T\), \(\[Sigma]\[Nu]\)]\)+\!\(\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Sigma]\[Mu]\), \(\[Nu]\)]\)\!\(\*SuperscriptBox[\(T\), \(\[Mu]\[Sigma]\)]\)="ward];


(* ::Section::Initialization::Closed:: *)
(*(*(*Equations of motion in characteristic form*)*)*)


(* ::Input::Initialization:: *)
(* replacements for the modified derivatives: \!\(
\*SubscriptBox[\(\[PartialD]\), \(t\)]f\)\[Rule]Overscript[f, .]-1/2A*\!\(
\*SubscriptBox[\(\[PartialD]\), \(r\)]f\) *)
rule1t=Flatten[Table[{
\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->D[\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]-1/2 A[r,t] \!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],{r,ddr}],
\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->D[\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]-1/2 A[r,t] \!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],{r,ddr}],
\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->D[\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]-1/2 A[r,t] \!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],{r,ddr}],
\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->D[\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]-1/2 A[r,t] \!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],{r,ddr}],
\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->D[\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]-1/2 A[r,t] \!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],{r,ddr}],
\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->D[\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]-1/2 A[r,t] \!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],{r,ddr}]
},{ddr,0,2}]];
rule2t={
\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->D[\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]-1/2 A[r,t] \!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],t],
\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->D[\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]-1/2 A[r,t] \!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],t],
\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->D[\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]-1/2 A[r,t] \!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],t]
};


(* ::Input::Initialization:: *)
(* transforming the Einstein equations to characteristic form *)
EOM1=EE[[1]]/.rule2t/.rule1t;
replEOM1=Solve[EOM1==0,\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]]//Flatten;


(* ::Input::Initialization:: *)
EOM2=EE[[3]]/.rule2t/.rule1t//FullSimplify;
replEOM2=Solve[EOM2==0,\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]]/.replEOM1//FullSimplify//Flatten;


(* ::Input::Initialization:: *)
EOMaux=EE[[4]]/.rule2t/.rule1t//Simplify;
replEOMaux=Solve[EOMaux==0,\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]]/.replEOM2/.replEOM1//FullSimplify//Flatten;


(* ::Input::Initialization:: *)
EOM3=EE[[2]]/.rule2t/.rule1t//Simplify;
replEOM3=Solve[EOM3==0/.replEOM2/.replEOM1/.replEOMaux,\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]]//FullSimplify//Flatten;


(* ::Input::Initialization:: *)
EOM4=EOMaux;
replEOM4=replEOMaux/.replEOM2/.replEOM1/.replEOM3//FullSimplify//Flatten;


(* ::Input::Initialization:: *)
EOM5=EE[[5]]/.rule2t/.rule1t//Simplify;
replEOM5=Solve[EOM5==0/.replEOM1/.replEOM2,\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]]//FullSimplify//Flatten;


(* ::Input::Initialization:: *)
(* it is important to note that in EOMList time-derivative means dot-derivative *)
solEE=Flatten[{replEOM1,replEOM2,replEOM3,replEOM4,replEOM5}/.Rule->Equal];
replTex1={A[r,t]->A,S[r,t]->S,B[r,t]->B,\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[_,_]->f',\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[_,_]->f'',\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[_,_]->\!\(\*OverscriptBox[\(f\), \(.\)]\)',\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[_,_]->\!\(\*OverscriptBox[\(f\), \(.\)]\),\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[_,_]->\!\(\*OverscriptBox[\(f\), \(\[DoubleDot]\)]\)};
solEEprint=solEE//.replTex1//Expand//TableForm;
Export["chareqns.dat",solEEprint];


(* ::Section::Initialization::Closed:: *)
(*(*(*Removing singularities*)*)*)


(* ::Input::Initialization:: *)
(* loading the near boundary solution *)
(*SetDirectory[NotebookDirectory[]];*)
order=5;resFile="NBsolAdS4Ord5.res";
nbSol=Get[resFile];
Clear[Ab,Sb,Bb];
Ab[z_,t_]=Asb[z,t,order]/.\[Epsilon]->1/.nbSol;
Sb[z_,t_]=Ssb[z,t,order]/.\[Epsilon]->1/.nbSol;
Bb[z_,t_]=Bsb[z,t,order]/.\[Epsilon]->1/.nbSol;
replBdry={A->Ab,S->Sb,B->Bb};


(* ::Input::Initialization:: *)
(* leading parts that are known explicitly from the near boundary expansion *)
subtractA[z_,t_]=Normal[Series[A[z,t]/.replBdry,{z,0,0}]];
subtractS[z_,t_]=Normal[Series[S[z,t]/.replBdry,{z,0,2}]];
subtractB[z_,t_]=Normal[Series[B[z,t]/.replBdry,{z,0,2}]];
(* here we write the dot derivative explicitly in the z-coordinate: Sdot=dS/dt-z^2/2A*dS/dz,... *)
subtractSdot[z_,t_]=Normal[Series[D[S[z,t],t]-1/2z^2A[z,t]D[S[z,t],z]/.replBdry,{z,0,0}]];
subtractBdot[z_,t_]=Normal[Series[D[B[z,t],t]-1/2z^2A[z,t]D[B[z,t],z]/.replBdry,{z,0,1}]];
(* leading Log contributions *)
subtractLogA[z_,t_]=D[Normal[Series[A[z,t]/.replBdry,{z,0,3}]],Log[z]]Log[z];
subtractLogS[z_,t_]=D[Normal[Series[S[z,t]/.replBdry,{z,0,5}]],Log[z]]Log[z];
subtractLogB[z_,t_]=D[Normal[Series[B[z,t]/.replBdry,{z,0,5}]],Log[z]]Log[z];
subtractLogSdot[z_,t_]=D[Normal[Series[D[S[z,t],t]-1/2z^2A[z,t]D[S[z,t],z]/.replBdry,{z,0,3}]],Log[z]]Log[z];
subtractLogBdot[z_,t_]=D[Normal[Series[D[B[z,t],t]-1/2z^2A[z,t]D[B[z,t],z]/.replBdry,{z,0,4}]],Log[z]]Log[z];


(* ::Input::Initialization:: *)
(* replacements for switching to finite variables Af,Sf,Bf,... *)
replfunc={
A->(Evaluate[subtractA[#1,#2]+subtractLogA[#1,#2]+Af[#1,#2]#1^1]&),
S->(Evaluate[subtractS[#1,#2]+subtractLogS[#1,#2]+Sf[#1,#2]#1^3]&),
B->(Evaluate[subtractB[#1,#2]+subtractLogB[#1,#2]+Bf[#1,#2]#1^3]&),
Sdot->(Evaluate[subtractSdot[#1,#2]+subtractLogSdot[#1,#2]+Sdotf[#1,#2]#1^1]&),
Bdot->(Evaluate[subtractBdot[#1,#2]+subtractLogBdot[#1,#2]+Bdotf[#1,#2]#1^2]&)
};


(* ::Input::Initialization:: *)
(* define list of functions which we want to evolve numerically and assign them numbers B=1,S=2,... *)
funcs={B,S,Sdot,Bdot,A,Sdotdt};
funcsf={Bf,Sf,Sdotf,Bdotf,Af,Sdotdtf};
{nB,nS,nSdot,nBdot,nA,nSdotdt}={1,2,3,4,5,6};


(* ::Input::Initialization:: *)
(* introduce dotted variables and multiplying the Einstein equations with powers of S to make them simpler *)
repldot={\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr_", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r_,t]->\!\(\*SuperscriptBox[\(Bdot\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr_", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r_,t]->\!\(\*SuperscriptBox[\(Sdot\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddr", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t],\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]->\!\(\*SuperscriptBox[\(Sdot\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]+1/2 A[r,t] \!\(\*SuperscriptBox[\(Sdot\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r,t]};
repldotback={Sdot->(\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[#1,#2]+1/2 A[#1,#2]\!\(\*SuperscriptBox[\(S\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[#1,#2]&),Bdot->(\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[#1,#2]+1/2 A[#1,#2]\!\(\*SuperscriptBox[\(B\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[#1,#2]&)};
multp={1,S[r,t],S[r,t],S[r,t]^2,1};
newEE=multp(solEE/.Equal->Subtract)/.repldot//Expand//FullSimplify;
newEEz=newEE/.replrtoz//FullSimplify;


(* ::Input::Initialization:: *)
(* Einstein equations in terms of finite variables *)
newEEfin=newEEz/.replfunc//Simplify;


(* ::Input::Initialization:: *)
(* determining the leading order in z of the Einstein equations *)
Series[newEEfin[[1]],{z,0,5}];
Series[newEEfin[[2]],{z,0,2}];
Series[newEEfin[[3]],{z,0,2}];
Series[newEEfin[[4]],{z,0,1}];
Series[newEEfin[[5]],{z,0,0}];


(* ::Input::Initialization:: *)
(* multiply the Einstein equations with the appropriate order of z such that the equations become finite at z=0  *)
powerz={z^-5,z^-2,z^-2,z^-1,z^0};
newEEf=powerz*newEEfin//Simplify;


(* ::Input::Initialization:: *)
Export["neweef.dat",newEEf];


(* ::Section::Initialization::Closed:: *)
(*(*(*Making equations numerical*)*)*)


(* ::Input::Initialization:: *)
(* deactivating some warnings that would pop up in the following replacement rule *)
Off[General::"partd"];


(* ::Input::Initialization:: *)
(* defining array "nums" that contains the numerical values of the metric functions at the grid points *)
replnum ={
Bf[z,t]->nums[[1,1]],\!\(\*SuperscriptBox[\(Bf\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[1,2]],\!\(\*SuperscriptBox[\(Bf\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[1,3]],
Sf[z,t]->nums[[2,1]],\!\(\*SuperscriptBox[\(Sf\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[2,2]],\!\(\*SuperscriptBox[\(Sf\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[2,3]],
Sdotf[z,t]->nums[[3,1]],\!\(\*SuperscriptBox[\(Sdotf\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[3,2]],\!\(\*SuperscriptBox[\(Sdotf\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[3,3]],
Bdotf[z,t]->nums[[4,1]],\!\(\*SuperscriptBox[\(Bdotf\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[4,2]],\!\(\*SuperscriptBox[\(Bdotf\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[4,3]],
Af[z,t]->nums[[5,1]],\!\(\*SuperscriptBox[\(Af\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[5,2]],\!\(\*SuperscriptBox[\(Af\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[5,3]],
Sdotdtf[z,t]->nums[[6,1]],\!\(\*SuperscriptBox[\(Sdotdtf\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[6,2]],\!\(\*SuperscriptBox[\(Sdotdtf\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]->nums[[6,3]],
Subscript[a, 3][t]->a3,Subscript[b, 3][t]->b3,Derivative[1][Subscript[b, 3]][t]->b3dt
};
(* replacements for the boundary values of the fields *)
replnumb = replnum /. {z -> 0, nums[[a_, b_]] -> nums[[a, b, -1]]};//Quiet;
(* replacements for the boundary source *)
replb0 = Prepend[Table[\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{
RowBox[{"Derivative", "[", "i", "]"}], "[", "b0", "]"}], "[", "t", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) ->  b0num[[i + 1]], {i, 7}], b0[t] -> b0num[[1]]];


(* ::Input::Initialization:: *)
(* updating the derivatives of the numerical functions using the spectral matrices  *)
update[n_]:=(
nums[[n,2]]=DSpec . nums[[n,1]];
nums[[n,3]]=D2Spec . nums[[n,1]];
b3=nums[[nB,1,-1]];
);
SetAttributes[update,Listable];


(* ::Input::Initialization:: *)
(* collecting the variables in the order in which they appear in the nested scheme *)
solvefor={Sf,Sdotf,Bdotf,Af};


(* ::Input::Initialization:: *)
(* building numerical arrays for the coefficient functions in the Einstein equations *)
Do[
Subscript[As, i]=D[newEEf[[i]],\!\(\*SuperscriptBox[\(solvefor\[LeftDoubleBracket]i\[RightDoubleBracket]\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t] ]//Expand//Simplify;
Subscript[Asb, i]=SeriesCoefficient[Subscript[As, i],{z,0,0}]//Expand//Simplify;
Subscript[Bs, i]=D[newEEf[[i]],\!\(\*SuperscriptBox[\(solvefor\[LeftDoubleBracket]i\[RightDoubleBracket]\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]]//Expand//Simplify;
Subscript[Bsb, i]=SeriesCoefficient[Subscript[Bs, i],{z,0,0}]//Expand//Simplify;
Subscript[Cs, i]=D[newEEf[[i]],solvefor[[i]][z,t]]//Expand//Simplify;
Subscript[Csb, i]=SeriesCoefficient[Subscript[Cs, i],{z,0,0}]//Expand//Simplify;
Subscript[Ss, i]=newEEf[[i]]/.solvefor[[i]]->(0&)//Expand//Simplify;
Subscript[Ssb, i]=SeriesCoefficient[Subscript[Ss, i],{z,0,0}]//Expand//Simplify;
,{i,Length[solvefor]}];


(* ::Input::Initialization:: *)
(* Replacing the analytic expressions with numerical arrays. *)
(* For example Subscript[As, 1] contains the analytic expression in front of the second derivative term in the first Einstein equation, *)
(* Subscript[As, 1+30] contains the same expression but with the analytic expression replaced by the numerical arrays. *)
makeformulas:=(
Do[Subscript[#, ii+30]=Subscript[#, ii]/.replnum/.replnumb/.replb0;,{ii,Length[solvefor]}]&/@{As,Asb,Bs,Bsb,Cs,Csb,Ss,Ssb};
Subscript[Bsb, 2+30]=0;
Subscript[Csb, 2+30]=1;
Subscript[Ssb, 2+30]=-SeriesCoefficient[D[S[z,t],t]-1/2z^2A[z,t]D[S[z,t],z]/.replBdry,{z,0,1}]/.Log[z]->0/.replnum/.replb0;
);
makeformulas;
(* The last Einstein equation is a constraint which we use to monitor the accuracy of the time evolution. *)
constraint=newEEf[[-1]]/.{\!\(\*SuperscriptBox[\(Sdotf\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddz_", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]-> \!\(\*SuperscriptBox[\(Sdotdtf\), 
TagBox[
RowBox[{"(", 
RowBox[{"ddz", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]}//.replnum//.replb0//Simplify;


(* ::Input::Initialization:: *)
(* saving the numerical formulas to a file such that they can be loaded fast once they were generated before *)
restable=(Table[Subscript[#, ii+30],{ii,Length[solvefor]}]&/@{As,Asb,Bs,Bsb,Cs,Csb,Ss,Ssb})~Join~{constraint}//Simplify;
Put[restable,"formulaAdS4.dat"];


(* ::Input::Initialization:: *)
(* loading the formulas from file *)
restable=Get["formulaAdS4.dat"];
Do[Subscript[{As,Asb,Bs,Bsb,Cs,Csb,Ss,Ssb}[[kk]], ii+30]=restable[[kk,ii]],{ii,Length[solvefor]},{kk,8}];
constraint=restable[[-1]];
(* determining the order of the differential equations *)
solveorder=If[Subscript[As, #+30]===0,If[Subscript[Bs, #+30]===0,0,1],2]&/@Range[Length[solvefor]];
Export["solveorder.dat",solveorder];


(* ::Input::Initialization:: *)
(* expressing the t-derivatives of B and b3 in terms of Bdot *)
dBfdt=(\!\(\*SuperscriptBox[\(Bf\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]/.Flatten[Solve[D[B[z,t],t]==Bdot[z,t]+z^2/2 A[z,t]D[B[z,t],z]/.replfunc,\!\(\*SuperscriptBox[\(Bf\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[z,t]]]);
dBfdtn=dBfdt/.replnum/.replb0;


(* ::Input::Initialization:: *)
dBfdtb=SeriesCoefficient[dBfdt,{z,0,0}];
dBfdtbn=dBfdtb/.replnum/.replnumb/.replb0;
b3dt=dBfdtbn;
Export["b3dt.dat",b3dt];


(* ::Input::Initialization:: *)
(* The expression for the t-derivative of a4 follows from conservation of the stress tensor. *)
(* We have obtained this expression from the near boundary expansion. *)
a3dtn=Derivative[1][Subscript[a, 3]][t]/.nbSol/.replnum/.replb0;
Export["a3dtn.dat",a3dtn];


(* ::Input::Initialization:: *)
(* VEV of the scalar operator which is the source for the YM sector *) 
VEVnum=VEV[t]/.replnum/.replb0;
VEVdt[t_]=D[VEV[t],t];
VEVdtnum=VEVdt[t]/.replnum/.replb0;
Export["vevdtnum.dat",VEVdtnum];


(* ::Section::Initialization::Closed:: *)
(*(*(*Initializing the Chebyshev grid and time stepping*)*)*)


(* ::Input::Initialization:: *)
(* the function "initGrid" builds the Chebyshev grid and the spectral differentialization matrices DSpec and D2Spec *)
(*initGrid[points_]:=(
nz=points;
z=N[Table[Cos[(r \[Pi])/(nz-1)],{r,0,nz-1}]];
DSpec=Table[If[i\[Equal]nz&&j\[Equal]nz,-((2(nz-1)^2+1)/6),If[i\[Equal]1&&j\[Equal]1,(2(nz-1)^2+1)/6,If[i\[Equal]j,-z[[i]]/(2(1-z[[i]]^2)),(-1)^(i+j)/(z[[i]]-z[[j]])If[i\[Equal]1||i\[Equal]nz,2,1]/If[j\[Equal]1||j\[Equal]nz,2,1]]]],{i,nz},{j,nz}];
z=z/2+1/2;
DSpec*=2;
D2Spec=DSpec.DSpec;z);*)
initGrid[points_,zMin_,zMax_]:=(
nz=points;
z=SetPrecision[Table[(1/2)((zMax+zMin)+(zMax-zMin)*Cos[Pi*r/(nz-1)]),{r,0,nz-1}],15];
cSpec={2}~Join~Table[1(-1)^(i),{i,nz-2}]~Join~{2*(-1)^(nz-1)};
Xspec=Transpose[Table[z,{i,nz}]];
dXspec=Xspec-Transpose[Xspec];
DSpec=Table[cSpec[[i]]*1/cSpec[[j]],{i,nz},{j,nz}]/(dXspec+IdentityMatrix[nz]);
DSpec=DSpec-IdentityMatrix[nz]*Table[Sum[Transpose[DSpec][[i,j]],{i,nz}],{j,nz}];
D2Spec=DSpec . DSpec;z);


(* ::Input::Initialization:: *)
(* initializing the time derivatives of b0[t] *)
InitTrepl[b0ini_,t_]:=Module[{},bini=b0ini;
auxb0=Prepend[Table[D[bini,{ti,i}],{i,6}],bini];
b0num=auxb0/.ti->t//Chop;];


(* ::Input::Initialization:: *)
(* setting the initial conditions for the time evolution *)
init[profileB_,a3ini_,b0ini_,tini_,npoints_,zMin_,zMax_]:=(
initGrid[npoints,zMin,zMax];
tinitial=t=tini;
nums=Table[0,{Length[funcs]},{3},{nz}];
a3=a3ini;
Bfini=Function[z,(profileB[z])(*(profileB[z]/z^4)*)(*/.replnum*)];
nums[[nB,1]]=Bfini/@z;
nums[[nB,1,-1]]=SeriesCoefficient[Series[Bfini[ztmp,tini],{ztmp,0,0}],0];
Bdtlist=Sdotlist=a3dtlist=timelist=a3list=b3list=b3dtlist=constrlist=VEVlist=VEVdtlist=Alist=Blist=Slist={};
b0numlist={b0num};
update[nB];
);


(* ::Input::Initialization:: *)
(* loop for the time evolution *)
order=4;
(* the initial time steps are smaller *)
dts=Table[10^(i/2),{i,-8,-1}]~Join~{1-Sum[10^(i/2),{i,-8,-1}]}~Join~{1};
step[change_,append_,b0ini_]:=(
InitTrepl[b0ini,t];
Do[
If[solveorder[[kk]]==2,Anum=Subscript[As, kk+30];Anum[[-1]]=Subscript[Asb, kk+30]];
If[solveorder[[kk]]>=1,Bnum=Subscript[Bs, kk+30];Bnum[[-1]]=Subscript[Bsb, kk+30]];
Cnum=Subscript[Cs, kk+30];Cnum[[-1]]=Subscript[Csb, kk+30];
Snum=Subscript[Ss, kk+30];Snum[[-1]]=Subscript[Ssb, kk+30];

eqm=Switch[solveorder[[kk]],
2,DiagonalMatrix[SparseArray[Anum]] . D2Spec+DiagonalMatrix[SparseArray[Bnum]] . DSpec+DiagonalMatrix[SparseArray[Cnum]],
1,DiagonalMatrix[SparseArray[Bnum]] . DSpec+DiagonalMatrix[SparseArray[Cnum]],
0,DiagonalMatrix[SparseArray[Cnum]]];

nums[[kk+1,1]]=LinearSolve[eqm,-Snum];
update[kk+1];
,{kk,1,Length[solvefor]}];

nums[[nSdotdt,1]]=If[Length[timelist]>5&&((timelist[[-6]]-timelist[[-5]])-(timelist[[-3]]-timelist[[-2]])<10^-12),
1/(12 dt) (3 Sdotlist[[-4]]-16 Sdotlist[[-3]]+36 Sdotlist[[-2]]-48 Sdotlist[[-1]]+25 nums[[nSdot,1]]),
Derivative[1][Interpolation[Thread[{timelist~Join~{t},Sdotlist~Join~{nums[[nSdot,1]]}}],InterpolationOrder->Min[order-1,Length[timelist]]]][t]];
update[nSdotdt];
connum=constraint[[;;-2]]~Join~{0};
mcon=Max[Abs[connum]];

If[append!= -1,
Bdtlist=Append[Drop[Bdtlist,append],dBfdtn[[;;-2]]~Join~{dBfdtbn}];
a3dtlist=Append[Drop[a3dtlist,append],a3dtn];
Sdotlist=Append[Drop[Sdotlist,append],nums[[nSdot,1]]];
Slist=Append[Drop[Slist,append],nums[[nS,1]]];
Blist=Append[Drop[Blist,append],nums[[nB,1]]];
Alist=Append[Drop[Alist,append],nums[[nA,1]]];
timelist=Append[timelist,t];
b3dtlist=Append[b3dtlist,dBfdtbn];
a3list=Append[a3list,a3];
b3list=Append[b3list,b3];
b0numlist=Append[b0numlist,b0num];
(*VEVlist=Append[VEVlist,VEVnum];*)
(*VEVdtlist=Append[VEVdtlist,VEVdtnum];*)
constrlist=Append[constrlist,mcon];
];

If[change,
t=t+dt dts[[Min[Length[dts],Length[timelist]]]];
If[Length[timelist]>5&&((timelist[[-6]]-timelist[[-5]])-(timelist[[-3]]-timelist[[-2]])<10^-12),
nums[[nB,1]]+=dt/24 (55 Bdtlist[[-1]]-59 Bdtlist[[-2]]+37 Bdtlist[[-3]]-9 Bdtlist[[-4]]);
a3+=dt/24 (55 a3dtlist[[-1]]-59 a3dtlist[[-2]]+37 a3dtlist[[-3]]-9 a3dtlist[[-4]]);,
nums[[nB,1]]+=\!\(
\*SubsuperscriptBox[\(\[Integral]\), \(timelist[\([\(-1\)]\)]\), \(t\)]\(\(Interpolation[Table[{timelist[\([\(-j\)]\)], Bdtlist[\([\(-j\)]\)]}, {j, Min[Length[Bdtlist], 6]}], InterpolationOrder -> Min[order - 1, Length[Bdtlist] - 1]]\)[T]\ \[DifferentialD]T\)\);
a3+=\!\(
\*SubsuperscriptBox[\(\[Integral]\), \(timelist[\([\(-1\)]\)]\), \(t\)]\(\(Interpolation[Table[{timelist[\([\(-j\)]\)], a3dtlist[\([\(-j\)]\)]}, {j, Min[Length[a3dtlist], 6]}], InterpolationOrder -> Min[order - 1, Length[a3dtlist] - 1]]\)[T]\ \[DifferentialD]T\)\);
];
update[nB];
];
);


(* ::Section::Initialization::Closed:: *)
(*(*(*Running*)*)*)


(* ::Input::Initialization:: *)
saverun[file_]:=Put[{dt,Bfini,b0numlist,nz,t,tinitial,Bdtlist,Sdotlist,timelist,a3list,a3dtlist,b3list,b3dtlist,constrlist,Alist,Blist,Slist},file];
loadrun[file_]:=({dt,Bfini,b0numlist,nz,t,tinitial,Bdtlist,Sdotlist,timelist,a3list,a3dtlist,b3list,b3dtlist,constrlist,Alist,Blist,Slist}=Get[file];
initGrid[nz,zMin,zMax]);


(* ::Input::Initialization:: *)
plotfunc[n_]:=Show[plj[numsstart[[n,1]]],pl[nums[[n,1]]],PlotLabel->funcsf[[n]],ImageSize->350,BaseStyle->20,AxesLabel->"z",PlotRange->{Min[numsstart[[n,1]]~Join~nums[[n,1]]],Max[numsstart[[n,1]]~Join~nums[[n,1]]]}];
run[profile_,a3ini_,b0ini_,tini_,tend_,gridpoints_,zMin_,zMax_,name_]:=(session=name;
init[profile,a3ini,b0ini,tini,gridpoints,zMin,zMax];
(*dt=1/(25 nz^2)//N;*)
(*dt=1/(2nz^2)//N;*)
dt=0.001`15;
ttotal=0;
(*Monitor[*)Do[ttotal+=First[Timing[step[True,(*If[Length[timelist]<20,0,1]*)0,b0ini]]];
If[stepnumber==4,numsstart=nums];
(*If[Mod[stepnumber,25]\[Equal]0,timeleft=ToString[(Round[(tend-tini)/dt]-stepnumber) ttotal/(stepnumber)//Round]<>" seconds left." ;];*)
(*If[Mod[stepnumber,1000]\[Equal]0,Print[{stepnumber,t,ttotal,mcon,a3,b3}];
(*saverun[name<>".dat"];*)];*)
If[Max[Abs[nums[[nB,1]]]]>10^4||(mcon[[2]]>5),(*Throw[Null];*)Goto[next];];,{stepnumber,Round[(tend-tini)/dt]}];(*,{t,mcon(*,timeleft*)(*,plots*)}];*)
saverun[name<>".dat"];
);


(* ::Section::Initialization::Closed:: *)
(*(*(*Run the simulation*)*)*)


(* ::Input::Initialization:: *)
(*SetDirectory[NotebookDirectory[]];*)
randomtable1=ToExpression[Import["randomnumbers1.dat","TSV"]]//Flatten;
randomtable2=ToExpression[Import["randomnumbers2.dat","TSV"]]//Flatten;
randomtable3=ToExpression[Import["randomnumbers3.dat","TSV"]]//Flatten;


(* ::Input::Initialization:: *)
(*SetDirectory[NotebookDirectory[]];*)
(* numerical parameters *)
ClearAll[t,z];
Btable=ToExpression[Import["Btable.dat","TSV"]];


(* ::Input::Initialization:: *)
exprToFunction[expr_,vars_]:=ToExpression[ToString[FullForm[expr]/.MapIndexed[#1->Slot@@#2&,vars]]<>"&"];


(* ::Input::Initialization:: *)
(* numerical parameters *)
ClearAll[t,z];
tstart=0; tend=3;Ngrid=60;IntOrder=8;zMin=0;zMax=1.07`15;
(* initialconditions *)
a3ini=-1;
tot=1;
b0now[t_]:=0;


(* ::Input::Initialization:: *)
Do[
(*Label[begin];*)
ClearAll[t,z];
a3now=a3ini;
fileName="ml-test-datafile-"<>ToString[i];
(*Print["profile"<>ToString[i]];*)
run[exprToFunction[Btable[[i,1]],{z}],a3now,b0now[ti],tstart,tend,Ngrid,zMin,zMax,fileName];//Quiet;
Label[next];
Clear[t];,{i,1,Length[Btable]}];
