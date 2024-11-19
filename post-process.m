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
(*(*(*(*(*(*(*Running*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
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
saverun[file_]:=Put[{dt,Bfini,b0numlist,nz,t,tinitial,Bdtlist,Sdotlist,timelist,a3list,a3dtlist,b3list,b3dtlist,constrlist,Alist,Blist,Slist},file];
loadrun[file_]:=({dt,Bfini,b0numlist,nz,t,tinitial,Bdtlist,Sdotlist,timelist,a3list,a3dtlist,b3list,b3dtlist,constrlist,Alist,Blist,Slist}=Get[file];
initGrid[nz,0,1.1`15]);


(* ::Section::Initialization::Closed:: *)
(*(*(*(*(*(*(*Extract Event and Apparent horizon data*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
kk=0;
Do[
(* loading results for IC1 *)
(*SetDirectory[NotebookDirectory[]];*)
If[FileExistsQ["ml-raw-datafile-"<>ToString[k]<>".dat"]==True,
kk =kk + 1;
loadrun["ml-raw-datafile-"<>ToString[k]<>".dat"];
a3int=Interpolation[Thread[{timelist,a3list}]];
b3int=Interpolation[Thread[{timelist,b3list}]];
energylist = Thread[{timelist[[;;]],Abs[-2a3list[[;;]]]}];
anisolist =Thread[{timelist[[;;]],Abs[6*b3list[[;;]]]}];
epsilon[t_]:=-2 a3int[t];
pl[t_]:=- a3int[t]-3b3int[t];
pt[t_]:=- a3int[t]+3b3int[t];
(* locating the apparent horizon *)
zAHList={};
Do[
SdotNow=Interpolation[Thread[{z[[;;]],Sdotlist[[i,;;]]}]];
b0dtNow=b0numlist[[i,2]];
ThetaNow[z_]:=1/2+z^3SdotNow[z]-5/16z^2b0dtNow^2;
zAHnow=zz/.FindRoot[ThetaNow[zz],{zz,1},WorkingPrecision->15]//Quiet;
AppendTo[zAHList,{timelist[[i]],zAHnow}];
,{i,1,Length[timelist]}];
b0int=Interpolation[Thread[{timelist[[;;]],b0numlist[[;;-2,1]]}]];
b0dtint=Interpolation[Thread[{timelist[[;;]],b0numlist[[;;-2,2]]}]];
SData=Flatten[Table[{{z[[i]],timelist[[j]]},Slist[[j,i]]},{i,1,Length[z]},{j,1,Length[timelist]}],1];
Sint=Interpolation[SData];
SAHlistcorrected=Table[{zAHList[[i,1]],1/(2Pi)(1/zAHList[[i,2]]+zAHList[[i,2]]^3Sint[zAHList[[i,2]],zAHList[[i,1]]]-1/8 zAHList[[i,2]] b0dtint[zAHList[[i,1]]]^2)^2},{i,1,Length[zAHList]}];
Alistfull={};Aintlist={};
Do[
Anowlist={timelist[[i]],Thread[{z[[;;]],Alist[[i,;;]]}]};
Clear[Aint];
Aint= Interpolation[Thread[{z[[;;]],Alist[[i,;;]]}]];
AppendTo[Alistfull,Anowlist];
AppendTo[Aintlist,Aint];
,{i,1,Length[timelist]}];
zEHlist={};Clear[zEHnow];AppendTo[zEHlist,{timelist[[-1]],(-a3list[[-1]])^(-1/3)}];
Do[
Clear[zEHnow];
zEHnow = zEHlist[[i,2]] - dt*(-1/2 - (zEHlist[[i,2]]^3/2 )*Aintlist[[-i]][zz]/.zz-> zEHlist[[i,2]]);
AppendTo[zEHlist,{timelist[[-i-1]],zEHnow}];,
{i,1,Length[timelist]-1}];
SEHlistcorrected=Table[{zEHlist[[i,1]],1/(2Pi)(1/zEHlist[[i,2]]+zEHlist[[i,2]]^3Sint[zEHlist[[i,2]],zEHlist[[i,1]]]-1/8 zEHlist[[i,2]] b0dtint[zEHlist[[i,1]]]^2)^2},{i,1,Length[zEHlist]}];
Put[{dt,Bfini,b0numlist,nz,t,tinitial,Bdtlist,Sdotlist,timelist,a3list,a3dtlist,b3list,b3dtlist,constrlist,Alist,Blist,Slist,Reverse[zEHlist][[;;,2]],zAHList[[;;,2]],Reverse[SEHlistcorrected][[;;,2]],SAHlistcorrected[[;;,2]]},"ml-data-pp-"<>ToString[kk]<>".dat"];
(*Print["kk="<>ToString[kk]];
Print["k="<>ToString[k]];*)],{k,1,20000}];
