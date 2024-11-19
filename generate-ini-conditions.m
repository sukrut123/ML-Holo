(* ::Package:: *)

(* ::Input::Initialization:: *)
randomtable1 =Table[RandomReal[{0,1},WorkingPrecision->30],{i,1,100}];


(* ::Input::Initialization:: *)
randomtable2 =Table[RandomReal[{0,1},WorkingPrecision->30],{i,1,100}];


(* ::Input::Initialization:: *)
randomtable3 =Table[RandomReal[{0,1},WorkingPrecision->30],{i,1,100}];


(* ::Input::Initialization:: *)
Export["randomnumbers1.dat",randomtable1 ];


(* ::Input::Initialization:: *)
Export["randomnumbers2.dat",randomtable2 ];


(* ::Input::Initialization:: *)
Export["randomnumbers3.dat",randomtable3 ];


(* ::Input::Initialization:: *)
exprToFunction[expr_,vars_]:=ToExpression[ToString[FullForm[expr]/.MapIndexed[#1->Slot@@#2&,vars]]<>"&"];


(* ::Input::Initialization:: *)
(*Create a smooth Gaussian noise function*)
smoothNoise[x_,terms_:10,stdDev_:0.01`30]:=Sum[RandomReal[NormalDistribution[0,stdDev],WorkingPrecision->30]*Sin[RandomReal[{0,2 Pi},WorkingPrecision->30]*x+RandomReal[{0,2 Pi},WorkingPrecision->30]],{i,terms}];


(* ::Input::Initialization:: *)
ClearAll[t,z];
Btable={};noisyBtable={};ClearAll[z];
Do[
ClearAll[z];
atemp=RandomChoice[randomtable1];
btemp=RandomChoice[randomtable2];
ctemp=RandomChoice[randomtable3];
temp =atemp*Exp[-200 (z-0.2`30)^2] +btemp*Exp[-200 (z-0.5`30)^2]  +ctemp*Exp[-100 (z-0.8`30)^2];
(*Add the smooth noise function to each function in the dataset*)noisytemp=temp+smoothNoise[z,50];
AppendTo[Btable,temp ];
AppendTo[noisyBtable,noisytemp ];,{j,1,20000}];
Export["Btable.dat",Btable];
Export["noisyBtable.dat",noisyBtable];
