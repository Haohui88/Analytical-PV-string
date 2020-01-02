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

If[ Not@ValueQ[PlotIV::usage],
PlotIV::usage = "plot the IV curve with key points labelled."]

If[ Not@ValueQ[ShadeAreaFraction::usage],
ShadeAreaFraction::usage = "performs simple estimation of area based inter-row shading fraction of a typical array (sheds)."]


(* ::Chapter:: *)
(*Main definitions*)


Begin["`Private`"];

Off[InterpolatingFunction::dmval];


(* ::Section:: *)
(*IV curve functions*)


(* ::Text:: *)
(*StringIV requires a set of IV curves in the format of {{current 1, voltage 1}, {current 2, voltage 2}...}. *)


StringIV[IVset_,opt:OptionsPattern[]]:=Module[{maxJ,bypass=OptionValue["BypassDiode"],diodeVoltage,combinedIV,IV$fleetInterp,probe},
If[bypass==True,maxJ=Max@IVset[[All,-1,1]],maxJ=Min@IVset[[All,-1,1]]];
diodeVoltage=OptionValue["BypassDiodeVoltage"];

combinedIV=With[{range=Range[0,maxJ]},
	{range,Total@Table[
		With[{interp=Interpolation[DeleteDuplicatesBy[Round[IVset[[i]],0.00001],First],InterpolationOrder->1]},
			Table[With[{interpPt=interp[x]},If[bypass==True&&interpPt<-diodeVoltage,-diodeVoltage,interpPt]],{x,range}]
		]
	,{i,Length@IVset}]}\[Transpose]];

IV$fleetInterp=Interpolation[#,InterpolationOrder->1]&/@IVset;
maxJ=combinedIV[[-1,1]];
probe=maxJ*1.01;

While[combinedIV[[-1,2]]>-diodeVoltage*Length@IVset*1.1 && probe<maxJ*1.1,
	probe=combinedIV[[-1,1]]*1.01;
	AppendTo[combinedIV,{probe,Total@Table[
		With[{interpPt=IV$fleetInterp[[i]][probe]},If[bypass==True&&interpPt<-diodeVoltage,-diodeVoltage,interpPt]]
		,{i,Length@IVset}]}
	];
];

Return[combinedIV];

];


(* ::Text:: *)
(*MPPT requires a set of IV curves in the format of {{current 1, voltage 1}, {current 2, voltage 2}...}. *)
(*returns {current, voltage, power} at MPP. *)


MPPT[IV_,opt:OptionsPattern[]]:=Module[{method,IVP,interp,maxJ,MPP,x},
method=OptionValue[TrackingMethod];

If[method=="Fast",
	IVP=Append[#,Times@@#]&/@IV;
	MPP=First[IVP~Extract~Position[#,Max@#]&@Part[IVP\[Transpose],3]];
];

If[method=="Robust",
	interp=Interpolation[IV,InterpolationOrder->1];
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


PlotIV[IV_]:=Module[{Isc,Voc,mpp,plot1,plot2},

Isc=Interpolation[DeleteDuplicatesBy[Reverse/@IV,First],InterpolationOrder->1][0];
Voc=Interpolation[DeleteDuplicatesBy[IV,First],InterpolationOrder->1][0];

mpp=MPPT@IV;

plot1=ListLinePlot[Reverse/@IV,AxesLabel->{"V","I"}];
plot2=ListPlot[{Callout[{Voc,0},"Voc = "<>ToString@Round[Voc,0.1]],Callout[{0,Isc},"Isc = "<>ToString@Round[Isc,0.1]],Callout[Reverse@mpp[[;;2]],"Vmpp = "<>ToString@Round[mpp[[2]],0.1]<>", Impp = "<>ToString@Round[mpp[[1]],0.1],Left]},
PlotStyle->{Red,PointSize[Large]},Epilog->{Text["FF = "<>ToString@Round[Last@mpp/(Isc*Voc),0.01],{Voc/2,Isc/2},{0,-1}],Text["Power = "<>ToString@Round[Last@mpp,0.1],{Voc/2,Isc/2},{0,1}]}];

Return[Show[plot2,plot1]];

];


(* ::Section:: *)
(*Auxiliary functions*)


(* ::Text:: *)
(*ShadeAreaFraction performs simple estimation of area based inter-row shading fraction of a typical array (sheds). *)


ShadeAreaFraction[d_,L_,w_,tilt_,orientation_,sunElev_,sunAzimuth_,n_:5,m_:10]:=Module[{\[Theta],rowFraction,widthFraction,shadeFraction},
\[Theta]=sunAzimuth-orientation;(*\[Theta] is transformed polar coordinate, + direction is clockwise, -ve direction is counter-clockwise, -90 to 0 is towards left, 0 to 90 is towards right, when facing forward*)
If[\[Theta]>270,\[Theta]=\[Theta]-360];
If[\[Theta]<-270,\[Theta]=\[Theta]+360];

rowFraction=If[Abs[\[Theta]]>=90,
0,
If[\[Theta]<0,
	Min[N[(n L-Min[d Tan[Abs[\[Theta]] \[Degree]],n L])/L],1],(*sunlight is from the left*)
	Min[N[((m-n+1) L-Min[d Tan[Abs[\[Theta]] \[Degree]],(m-n+1) L])/L],1](*sunlight is from the right*)
	]
];
widthFraction=1-Min[1/w*(d Cos[tilt \[Degree]]-d Sin[tilt \[Degree]] Tan[90\[Degree]-tilt \[Degree]-sunElev \[Degree]]),1]//N;
shadeFraction=rowFraction*widthFraction;
Return[shadeFraction];

];


(* ::Chapter:: *)
(*End*)


End[]; 

EndPackage[]; 
