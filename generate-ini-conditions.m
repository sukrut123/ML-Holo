(* ::Package:: *)

(* ::Input::Initialization:: *)
randomtable1 =Table[RandomReal[{0,1}],{i,1,100}];


(* ::Input::Initialization:: *)
randomtable2 =Table[RandomReal[{0,1}],{i,1,100}];


(* ::Input::Initialization:: *)
randomtable3 =Table[RandomReal[{0,1}],{i,1,100}];


(* ::Input::Initialization:: *)
Export["randomnumbers1.dat",randomtable1 ];


(* ::Input::Initialization:: *)
Export["randomnumbers2.dat",randomtable2 ];


(* ::Input::Initialization:: *)
Export["randomnumbers3.dat",randomtable3 ];


(* ::Input::Initialization:: *)
exprToFunction[expr_,vars_]:=ToExpression[ToString[FullForm[expr]/.MapIndexed[#1->Slot@@#2&,vars]]<>"&"];


(* ::Input::Initialization:: *)
ClearAll[t,z];
Btable={};ClearAll[z];
Do[
ClearAll[z];
atemp=RandomChoice[randomtable1];
btemp=RandomChoice[randomtable2];
ctemp=RandomChoice[randomtable3];
temp =atemp*Exp[-200 (z-0.2`15)^2] +btemp*Exp[-200 (z-0.5`15)^2]  +ctemp*Exp[-100 (z-0.8`15)^2];
AppendTo[Btable,temp ];,{j,1,10000}];
Export["Btable.dat",Btable];
