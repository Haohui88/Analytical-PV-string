(* ::Package:: *)

(* ::Title:: *)
(*Analytical PV strings*)


(* ::Chapter:: *)
(*Initialization*)


(* :Copyright: Haohui Liu  *)

(* :Name: Analytical PV strings *)

(* :Author: Dr. Liu Haohui *)

(* :Package Version: 1.0 *)

(* :Mathematica Version: 12.0 *)

(*:Summary:
	provides analytical models for PV strings to calculate DC output.
*)

BeginPackage["AnalyticalPVStrings`"];


If[ Not@ValueQ[StringIV::usage],
StringIV::usage = "combines individual IV curves into string IV curve assuming series connection."]

Options[StringIV]={"BypassDiode"->True,"BypassDiodeVoltage"->0.5};

If[ Not@ValueQ[CombinerIV::usage],
CombinerIV::usage = "combines individual IV curves into string IV curve assuming parallel connection."]

Options[CombinerIV]={};

If[ Not@ValueQ[MPPT::usage],
MPPT::usage = "extract maximum power point in the IV curve."]

Options[MPPT]={TrackingMethod->"Fast"};


(* ::Chapter:: *)
(*Main definitions*)


Begin["`Private`"];

Off[InterpolatingFunction::dmval];


(* ::Text:: *)
(*StringIV requires a set of IV curves in the format of {{current 1, voltage 1}, {current 2, voltage 2}...}*)


StringIV[IVset_,opt:OptionsPattern[]]:=Module[{maxJ,diodeVoltage,combinedIV,IV$fleetInterp,probe},
If[OptionValue["BypassDiode"]==True,maxJ=Max@IVset[[All,-1,1]],maxJ=Min@IVset[[All,-1,1]]];
diodeVoltage=OptionValue["BypassDiodeVoltage"];

combinedIV=With[{range=Range[0,maxJ]},
	{range,Total@Table[
		With[{interp=Interpolation[DeleteDuplicatesBy[Round[IVset[[i]],0.00001],First],InterpolationOrder->1]},
			Table[With[{interpPt=interp[x]},If[OptionValue["BypassDiode"]==True&&interpPt<-diodeVoltage,-diodeVoltage,interpPt]],{x,range}]
		]
	,{i,Length@IVset}]}\[Transpose]];

IV$fleetInterp=Interpolation[#,InterpolationOrder->1]&/@IVset;
maxJ=combinedIV[[-1,1]];
probe=maxJ*1.01;

While[combinedIV[[-1,2]]>-diodeVoltage*Length@IVset*1.1 && probe<maxJ*1.1,
	probe=combinedIV[[-1,1]]*1.01;
	AppendTo[combinedIV,{probe,Total@Table[
		With[{interpPt=IV$fleetInterp[[i]][probe]},If[OptionValue["BypassDiode"]==True&&interpPt<-diodeVoltage,-diodeVoltage,interpPt]]
		,{i,Length@IVset}]}
	];
];

Return[combinedIV];

];


MPPT[IV_,opt:OptionsPattern[]]:=Module[{method,IVP,interp,maxJ,MPP,x},
method=OptionValue[TrackingMethod];

IVP=Append[#,Times@@#]&/@IV;
MPP=First[IVP~Extract~Position[#,Max@#]&@Part[IVP\[Transpose],3]];

If[method=="Robust",
	interp=Interpolation[IV,InterpolationOrder->3];
	maxJ=IV[[-1,1]];
	MPP=NMaximize[{interp[x]*x,{0<=x<maxJ}},x,AccuracyGoal->5,PrecisionGoal->5];
	MPP={x/.#[[2,1]],interp[x/.#[[2,1]]],First@#}&@MPP;
];

Return[MPP];

];


(* ::Text:: *)
(*CombinerIV under development... *)


CombinerIV[IVset_,opt:OptionsPattern[]]:=Module[{minV,diodeVoltage,combinedIV,IV$fleetInterp,probe},

combinedIV=With[{range=Range[0,minV]},
	{range,Total@Table[
		With[{interp=Interpolation[DeleteDuplicatesBy[Round[IVset[[i]],0.00001],First],InterpolationOrder->1]},
			Table[With[{interpPt=interp[x]},If[OptionValue["BypassDiode"]==True&&interpPt<-diodeVoltage,-diodeVoltage,interpPt]],{x,range}]
		]
	,{i,Length@IVset}]}\[Transpose]];
	
];


(* ::Chapter:: *)
(*End*)


End[]; 

EndPackage[]; 
