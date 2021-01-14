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

If[ Not@ValueQ[MPPT::usage],
MPPT::usage = "MPPT[IV] extract maximum power point in the IV curve. Possible TrackingMethods are: 
Fast (max point in the supplied IV curve);
FindPeak (uses function FindPeak to locate maximum, avoids returning MPP at boundary of MPPT voltage range);
Absolute (uses interpolation and guarantees absolute maximum even when IV list does not have extensive coverage of points, do not choose this if there is pLimit). 
MPPT[IV,vProbe] uses a sticky method, may use previous voltage as the initial probe point, may end up with local maximum. Also uses interpolation, not confined to finding a point already in the input IV list, as is the case with Fast and FindPeak methods."]

Options[MPPT]={"TrackingMethod"->"FindPeak","VoltageRange"->{150,1500},"SearchStep"->5,PowerLimit->None};

If[ Not@ValueQ[CablingCorrection::usage],
CablingCorrection::usage = "Modifies IV to include effect of cabling series resistance. 
Inputs are [IV curve, cable length (round trip), cable cross section (optional), resistivity (optional)].
CablingCorrection[cableLength,crossSection_:6,\[Rho]_:0.023] is the operator form. "]



If[ Not@ValueQ[ShadeAreaFraction::usage],
ShadeAreaFraction::usage = "ShadeAreaFraction[pitch, table length, collector width, tilt, orientation, sun elevation, sun azimuth, table position (defulat 5), table number (default 10)] 
performs simple estimation of area based inter-row shading fraction of a typical array (sheds)."]

If[ Not@ValueQ[ElectricalShading::usage],
ElectricalShading::usage = "ElectricalShading[areaShaded, numSubString, orientation] gives estimation on electrical shading by giving."]

Options[ElectricalShading]={"ModuleDimension"->{6,12}}; (* default standard module dimension is 6x12 cells *)


If[ Not@ValueQ[IAMashrae::usage],
IAMashrae::usage = "IAM[AOI,b0_0.05] ashrae model."]


If[ Not@ValueQ[PlotIV::usage],
PlotIV::usage = "PlotIV[IV,voltageRange(optional),PlotOptions] plots the IV curve with key points labelled.
voltageRange of MPPT by default is {150,1500} (change this when dealing with cell IV). PlotOptions of ListPlot can be specified to change plot display. "]

If[ Not@ValueQ[PlotPV::usage],
PlotPV::usage = "PlotPV[IV,voltageRange(optional),PlotOptions] plots the PV curve with key points labelled.
voltageRange of MPPT by default is {150,1500} (change this when dealing with cell IV). PlotOptions of ListPlot can be specified to change plot display. "]

If[ Not@ValueQ[GetIsc::usage],
GetIsc::usage = "GetIsc[IV] returns the Isc value of an IV curve. "]

If[ Not@ValueQ[GetVoc::usage],
GetVoc::usage = "GetVoc[IV] returns the Voc value of an IV curve. "]

If[ Not@ValueQ[GetFF::usage],
GetFF::usage = "GetFF[IV] returns the FF value of an IV curve. "]


(* ::Chapter:: *)
(*Main definitions*)


Begin["`Private`"];

Off[InterpolatingFunction::dmval];


(* ::Section:: *)
(*IV curve functions*)


(* ::Text:: *)
(*Overall convention: *)
(*IV curve functions require a set of IV curves in the format of {{current 1, voltage 1}, {current 2, voltage 2}...}. *)
(*The IV curves will be sorted by current in ascending order (from small to large, towards Isc, i.e. voltage from Voc, from +ve to -ve). *)
(*Range of the input IV curves should be complete (at least covering one whole quadrant and cross the two axes). *)


(* ::Text:: *)
(*StringIV takes care of combining IV in series both on string level or sub-module level. *)
(*IV with value Null is treated as completely blocked (zero output). *)


StringIV[inputIVset_/;ContainsAny[inputIVset,{Null}],opt:OptionsPattern[]]=Null;


StringIV[inputIVset_,opt:OptionsPattern[]]:=Module[{IVset,maxJ,bypass=OptionValue["BypassDiode"],diodeVoltage,combinedIV,IV$fleetInterp,probe},
IVset=SortBy[First]/@inputIVset;
If[bypass==True,maxJ=Max@IVset[[All,-1,1]];,maxJ=Min@Cases[IVset[[All,-1,1]],_?Positive];]; (*when bypass diodes are present, scan probing current to the max available in the IV set*)
If[!Positive@maxJ||!NumericQ@maxJ,Message[StringIV::invalidJ];Return[$Failed];Continue[];];
diodeVoltage=OptionValue["BypassDiodeVoltage"];

(*for same current, add up voltage*)
IV$fleetInterp=Interpolation[DeleteDuplicatesBy[Round[#,0.0001],First],InterpolationOrder->1]&/@IVset;
combinedIV=With[{range=Range[-0.2*maxJ,-0.04*maxJ,0.04*maxJ]~Join~Range[0,maxJ,Min[maxJ/100,1]]},
	{range,Total@Table[
		With[{interp=IV$fleetInterp[[i]]},
			Table[With[{interpPt=interp[x]},If[bypass==True&&interpPt<-diodeVoltage,-diodeVoltage,interpPt]],{x,range}]
			(*when bypass diode is present, negative voltage is pinned at diode voltage beyond Isc*)
		]
	,{i,Length@IVset}]}\[Transpose]];

(* extend the IV curve towards higher current region to cover the complete voltage range *)
probe=maxJ;

While[(combinedIV[[-1,2]]>-diodeVoltage*Length@IVset*1.1&&Abs@probe<Abs@maxJ*1.1)||combinedIV[[-1,2]]>0&&Abs@probe<Abs@maxJ*1.1,
(* as long as not all bypass diodes are activated in the IVset, i.e. voltage is not negatively biased enough *)
	probe=combinedIV[[-1,1]]*1.01;
	AppendTo[combinedIV,{probe,Total@Table[
		With[{interpPt=IV$fleetInterp[[i]][probe]},If[bypass==True&&interpPt<-diodeVoltage,-diodeVoltage,interpPt]]
		,{i,Length@IVset}]}
	];
];

Return[combinedIV];
];
StringIV::invalidJ="invalid short circuit current for IVs in the input set";


(* ::Text:: *)
(*CombinerIV  takes care of combining IV in parallel both on string level or sub-module level. *)
(*In practice, IV behavior in the negative voltage range does not involve complications brought by diodes for strings/sub-strings combining in parallel (strings with bypass diodes combining in parallel, then in series with others is rarely seen in practice). *)


CombinerIV[IVset_]:=Module[{minV,maxV,maxJ,diodeVoltage,combinedIV,IV$fleetInterp,probe},
maxV=Min@IVset[[All,1,2]];
minV=Min[Max@IVset[[All,-1,2]],0]; 
(*for IVs with bypass diode combining in parallel, negative voltage beyond diode voltages is not available, this part is either not interesting to look at in practice, or can be simply interpolated*)

(*for same voltage, add up current*)
IV$fleetInterp=Interpolation[DeleteDuplicatesBy[Round[#,0.0001],First],InterpolationOrder->1]&/@Map[Reverse/@#&,IVset];
combinedIV=With[{range=Range[maxV,minV,-Min[maxV/100,0.1]]},
	{Total@Table[
		With[{interp=IV$fleetInterp[[i]]},
			Table[interp[x],{x,range}]
		]
	,{i,Length@IVset}],range}\[Transpose]];

(* extend the IV curve towards higher voltage region to cover a more complete current range *)
probe=maxV;
maxJ=combinedIV[[-1,1]];

While[combinedIV[[1,1]]>-0.2*maxJ,
(* as long as current is still positive or not negative enough *)
	probe=combinedIV[[1,2]]*1.01;
	PrependTo[combinedIV,{Total@Table[IV$fleetInterp[[i]][probe],{i,Length@IVset}],probe}
	];
];

Return[combinedIV];
];


(* ::Text:: *)
(*MPPT requires a set of IV curves in the format of {{current 1, voltage 1}, {current 2, voltage 2}...}. *)
(*returns {current, voltage, power} at MPP. *)
(*Method using NMaximize is generally very slow. *)
(*The method "FindPeak" is somewhat more realistic in the sense that it is less likely to end up at MPP at the lower boundary of MPPT voltage range. *)


MPPT[Null,opt:OptionsPattern[]]={0,0,0};


MPPT[IV_,opt:OptionsPattern[]]:=Module[{method=OptionValue["TrackingMethod"],vRange=OptionValue["VoltageRange"],searchStep=OptionValue["SearchStep"],pLimit=OptionValue[PowerLimit],IVP,IV2,interp,maxJ,MPP={0,0,0},x,MPPsearch},

If[method=="Fast",
	IVP=Append[#,Times@@#]&/@IV;
	IVP=Select[IVP,#[[3]]>=0&&vRange[[1]]<#[[2]]<vRange[[2]]&];
	If[Length@IVP!=0,
		MPP=First[IVP~Extract~Position[#,Max@#]&@Part[IVP\[Transpose],3]];
		(*MPP=First@MaximalBy[IVP,Last];*)
		
		(* power limit *)
		If[NumericQ@pLimit&&pLimit<Last@MPP,
			Sow["power limited"];
			interp=Interpolation[DeleteDuplicatesBy[Round[Reverse/@IV,0.0001],First],InterpolationOrder->1]; (* interpolation: x\[Rule]voltage, y\[Rule]current *)
			MPPsearch=Pick[MovingAverage[IV[[All,2]],2],Negative/@Times@@@Partition[Times@@@IV-pLimit,2,1]]; (* pick out the voltages near pLimit *)
			MPPsearch=Select[x/.FindRoot[interp[x]*x==pLimit,{x,MPPsearch/.{}->{vRange[[2]]}}],vRange[[1]]<=#<=vRange[[2]]&];
			If[Length@MPPsearch==0,
				MPP={0,0,0};
			,
				MPPsearch=Max@MPPsearch; (* choose the larger voltage as the operating point *)
				MPP=With[{i=interp[MPPsearch]},{i,MPPsearch,i*MPPsearch}];
			];
		];
	];
];

If[method=="FindPeak",
	vRange={#1*1.025,#2}&@@vRange; (* modify lower boundary so that inverters are less likely to find the peak close to this end *)
	IV2=Select[IV,#[[1]]>=0&&vRange[[1]]<#[[2]]<vRange[[2]]&];
	
	If[Length@IV2!=0,
		interp=Interpolation[DeleteDuplicatesBy[Round[Reverse/@IV,0.0001],First],InterpolationOrder->1]; (* interpolation: x\[Rule]voltage, y\[Rule]current *)
		IV2=Append[{interp[vRange[[2]]],vRange[[2]]}]@Append[{interp[vRange[[1]]],vRange[[1]]}]@IV2; (* add the two boundary points of MPPT working range *)
		IV2=DeleteDuplicatesBy[Last]@SortBy[Last]@IV2; (* sorted by voltage, ascending order from small voltage to large *)
		MPPsearch=FindPeaks[Times@@@IV2]; (* gives {position in list IV2, peak power} *)
		If[Length@MPPsearch>1,
			MPPsearch=DeleteCases[MPPsearch,{1,_}]; (* remove the peak exactly at lower voltage boundary, if any *)
			MPPsearch=MaximalBy[MPPsearch,Last];
		];
		MPP=Append[IV2[[MPPsearch[[1,1]]]],MPPsearch[[1,2]]];
		
		(* power limit *)
		If[NumericQ@pLimit&&pLimit<Last@MPP,
			Sow["power limited"];
			MPPsearch=Pick[MovingAverage[IV[[All,2]],2],Negative/@Times@@@Partition[Times@@@IV-pLimit,2,1]]; (* pick out the voltages near pLimit *)
			MPPsearch=Select[x/.FindRoot[interp[x]*x==pLimit,{x,MPPsearch/.{}->{vRange[[2]]}}],vRange[[1]]<=#<=vRange[[2]]&];
			If[Length@MPPsearch==0,
				MPP={0,0,0};
			,
				MPPsearch=Max@MPPsearch; (* choose the larger voltage as the operating point *)
				MPP=With[{i=interp[MPPsearch]},{i,MPPsearch,i*MPPsearch}];
			];
		];
	];
];

If[method=="Absolute", (* guarantees absolute maximum even when IV list does not have extensive coverage of points, do not choose this if there is pLimit *)
	IV2=Select[IV,#[[1]]>=0&&vRange[[1]]<#[[2]]<vRange[[2]]&];
	
	If[Length@IV2!=0,
		interp=Interpolation[DeleteDuplicatesBy[Round[IV,0.0001],First],InterpolationOrder->1]; (* interpolation: x\[Rule]current, y\[Rule]voltage *)
		maxJ=MinMax@IV2[[All,1]];
		MPP=NMaximize[{interp[x]*x,{maxJ[[1]]<=x<=maxJ[[2]]}},x,AccuracyGoal->5,PrecisionGoal->5];
		MPP={x/.#[[2,1]],interp[x/.#[[2,1]]],First@#}&@MPP;
	];
];

Return[MPP];
];


MPPT[IV_,vProbe_,opt:OptionsPattern[]]:=Module[{vRange=OptionValue["VoltageRange"],pLimit=OptionValue[PowerLimit],IVP,IV2,interp,maxJ,MPP={0,0,0},x,MPPsearch,Vmax},

(* uses interpolation, not confined to finding a point already in the input IV list, as is the case with Fast and FindPeak methods *)
(* sticky method, may use previous voltage as the initial probe point, may end up with local maximum *)

	IV2=Select[IV,#[[1]]>=0&&vRange[[1]]<#[[2]]<vRange[[2]]&];
	If[Length@IV2!=0,
		(*Voc=First[IV[[All,2]]~Extract~Position[IV[[All,1]],Min@Abs@IV[[All,1]]]];*)
		interp=Interpolation[DeleteDuplicatesBy[Round[Reverse/@IV,0.0001],First],InterpolationOrder->1]; (* interpolation: x\[Rule]voltage, y\[Rule]current *)
		
		MPPsearch=Quiet@FindMaximum[{interp[x]*x,{vRange[[1]]<=x<=vRange[[2]]}},{x,vProbe},AccuracyGoal->3,PrecisionGoal->3]; (* gives {MPP, MPP voltage} *)
		Vmax=x/.Last@MPPsearch;
		MPP={interp[Vmax],Vmax,First@MPPsearch};
		
		(* power limit *)
		If[NumericQ@pLimit&&pLimit<Last@MPP,
			Sow["power limited"];
			MPPsearch=Pick[MovingAverage[IV[[All,2]],2],Negative/@Times@@@Partition[Times@@@IV-pLimit,2,1]]; (* pick out the voltages near pLimit *)
			MPPsearch=Select[x/.FindRoot[interp[x]*x==pLimit,{x,MPPsearch/.{}->{vRange[[2]]}}],vRange[[1]]<=#<=vRange[[2]]&];
			If[Length@MPPsearch==0,
				MPP={0,0,0};
			,
				MPPsearch=Max@MPPsearch; (* choose the larger voltage as the operating point *)
				MPP=With[{i=interp[MPPsearch]},{i,MPPsearch,i*MPPsearch}];
			];
		];
	];

Return[MPP];
];


(* ::Text:: *)
(*Cable voltage drop correction. Applies to string cables on the DC side. *)
(*Default cross section is 6mm^2. Default resistivity is 0.023 ohm.mm^2/m (copper), use 0.037 for aluminum. *)


CablingCorrection[IV_List/;ArrayDepth@IV==2,cableLength_?NumericQ,crossSection_:6,\[Rho]_:0.023]:=Module[{\[Delta]V},
\[Delta]V[current_]:=current*\[Rho]*cableLength/crossSection;

(* \[Delta]V is negative when current becomes negative *)
Return[{#1,#2-\[Delta]V[#1]}&@@@IV]; 
];


CablingCorrection[cableLength_?NumericQ,crossSection_:6,\[Rho]_:0.023][IV_List/;ArrayDepth@IV==2]:=Module[{\[Delta]V},
\[Delta]V[current_]:=current*\[Rho]*cableLength/crossSection;

(* \[Delta]V is negative when current becomes negative *)
Return[{#1,#2-\[Delta]V[#1]}&@@@IV]; 
];


(* ::Section::Closed:: *)
(*Shading logics*)


(* ::Text:: *)
(*ShadeAreaFraction performs simple estimation of area based inter-row shading fraction of a typical array (sheds). *)
(*Determines shading for direct beam component only. *)
(*d - pitch of arrays in m;*)
(*L - table length in m; *)
(*w - collector width in m;*)


ShadeAreaFraction[d_,L_,w_,tilt_,orientation_,sunElev_,sunAzimuth_,n_:1,m_:10]:=Module[{\[Phi],\[Theta]z,rowFraction,widthFraction,shadeFraction},
\[Phi]=sunAzimuth-orientation;(*\[Phi] is transformed polar coordinate, + direction is clockwise, -ve direction is counter-clockwise, -90 to 0 is towards left, 0 to 90 is towards right, when facing forward*)
If[\[Phi]>270,\[Phi]=\[Phi]-360];
If[\[Phi]<-270,\[Phi]=\[Phi]+360];

(* rowFraction: percentage of row length in lateral that is covered by shadow *)
rowFraction=If[Abs[\[Phi]]>=90,
0,
If[\[Phi]<0,
Min[N[(n L-Min[d Tan[Abs[\[Phi]] \[Degree]],n L])/L],1],(*sunlight is from the left*)
Min[N[((m-n+1) L-Min[d Tan[Abs[\[Phi]] \[Degree]],(m-n+1) L])/L],1](*sunlight is from the right*)
]
];

(* widthFraction: percentage of length in transverse collector direction that is covered by shadow *)
If[Abs[\[Phi]]>=90,
widthFraction=0;
,
\[Theta]z=ArcTan[Sin[sunElev \[Degree]]/(Cos[sunElev \[Degree]] Cos[\[Phi] \[Degree]])]/Degree;
widthFraction=Min[1-Min[1/w*(d Cos[tilt \[Degree]]-d Sin[tilt \[Degree]] Tan[90\[Degree]-tilt \[Degree]-\[Theta]z \[Degree]]),1],1]//N;
];

shadeFraction=rowFraction*widthFraction;
Return[{shadeFraction,rowFraction,widthFraction}];
];


(* ::Text:: *)
(*ElectricalShading uses a set of simple logic to estimate the amount of complete and partial shading of substrings (each PV module has three substring if there are three bypass diodes). *)
(*Applicable to standard c-Si modules with bypass diodes. *)
(*Requires shaded fraction and number of substrings in both lateral and transverse dimension (table configuration). Lengths in both directions will be measured as relative number (normalized to array table dimension). *)
(*nc: number of substrings completely shaded*)
(*nf: number of substrings fractionally shaded*)
(*n0: number of substrings completely unshaded*)
(*frac: fraction of shading*)


ElectricalShading[{shadeFracLat_?(0<=#<=1&),shadeFracTrans_?(0<=#<=1&)},{numSubStrLat_,numSubStrTrans_},orientation_:"landscape",opt:OptionsPattern[]]:=Module[{subStrLength,eShadeElementFn,modDimension=OptionValue["ModuleDimension"],n0,nc,ncLat,ncTrans,nfLat,nfTrans,remain,fracLat,fracTrans,x},

(* ------------ in the transverse direction ------------------- *)
subStrLength=1/numSubStrTrans; (* sub-string length in transverse direction *)

If[orientation=="landscape",
	eShadeElementFn=Piecewise[{{1/(0.5 subStrLength)*x,0<=x<0.5*subStrLength},{1,0.5*subStrLength<=x}},Null]; (* translate fraction of substring length shaded into electrical shading proportion *)
,
(* else, portrait *)
	eShadeElementFn=Piecewise[{{1/(1/modDimension[[2]]*subStrLength)*x,0<=x<1/modDimension[[2]]*subStrLength},{1,1/modDimension[[2]]*subStrLength<=x}},Null];
];

{ncTrans,remain}=QuotientRemainder[shadeFracTrans,subStrLength];
fracTrans=(eShadeElementFn/.x->remain);

(* ------------ in the lateral direction ------------------- *)
subStrLength=1/numSubStrLat; (* sub-string length in transverse direction *)

If[orientation=="portrait",
	eShadeElementFn=Piecewise[{{1/(0.5 subStrLength)*x,0<=x<0.5*subStrLength},{1,0.5*subStrLength<=x}},Null]; (* translate fraction of substring length shaded into electrical shading proportion *)
,
(* else, landscape *)
	eShadeElementFn=Piecewise[{{1/(1/modDimension[[2]]*subStrLength)*x,0<=x<1/modDimension[[2]]*subStrLength},{1,1/modDimension[[2]]*subStrLength<=x}},Null];
];

{ncLat,remain}=QuotientRemainder[shadeFracLat,subStrLength];
fracLat=(eShadeElementFn/.x->remain);

(* determine nf as a status value indicating there is partially shaded substring in lateral/transverse direction *)
Which[fracTrans==1,
	ncTrans=ncTrans+1;
	nfTrans=0;,
fracTrans<0.05,
	nfTrans=0;,
True,
	nfTrans=1;
];

Which[fracLat==1,
	ncLat=ncLat+1;
	nfLat=0;,
fracLat<0.05,
	nfLat=0;,
True,
	nfLat=1;
];

nfTrans=nfTrans*ncLat; (* nf updated as the number of substrings partially shaded *)
nfLat=nfLat*ncTrans;
nc=ncTrans*ncLat;
(* here ignores the substring containing the corner cell which is partially shaded in both directions *)
(*If[0.05<fracLat<1&&0.05<fracTrans<1,
	nfCorner=1;
	FracCorner=fracLat*fracTrans;
];*)

n0=numSubStrTrans*numSubStrLat-nc-nfTrans-nfLat;

Return@{nc,{{nfLat,fracLat},{nfTrans,fracTrans}},n0};

];


(*depreciated*)
(*ElectricalShading[areaShaded_,numSubString_,orientation_:"landscape",opt:OptionsPattern[]]:=Module[{subStringArea,eShadeElementFn,modDimension=OptionValue["ModuleDimension"],nCompleteShade,remain,fractionalShade,x},
subStringArea=1/numSubString; (* sub-string length in transverse direction *)

If[orientation=="landscape",
	eShadeElementFn=Piecewise[{{1/(0.5 subStringArea)*x,0<=x<0.5*subStringArea},{1,0.5*subStringArea<=x}},Null];
,
(* else, portrait *)
	eShadeElementFn=Piecewise[{{1/(1/modDimension[[2]]*subStringArea)*x,0<=x<1/modDimension[[2]]*subStringArea},{1,1/modDimension[[2]]*subStringArea<=x}},Null];
];

{nCompleteShade,remain}=QuotientRemainder[areaShaded,subStringArea];
fractionalShade=(eShadeElementFn/.x->remain);

Return@{nCompleteShade,fractionalShade};

];*)


(* ::Text:: *)
(*ElectricalShading in irregular shading mode... *)


(* ::Section::Closed:: *)
(*Incidence angle modifier*)


IAMashrae[AOI_,b0_:0.05]:=Max[1-b0*(1/Cos[AOI \[Degree]]-1),0]/;0<=AOI<90;
IAMashrae[AOI_?(!0<=#<90&),b0_:0.05]:=0;


(* ::Section:: *)
(*Auxiliary functions*)


GetVoc[IV_List]:=Interpolation[DeleteDuplicatesBy[IV,First@*(Round[#,0.0001]&)],InterpolationOrder->1][0];
GetIsc[IV_List]:=Interpolation[DeleteDuplicatesBy[Reverse/@IV,First@*(Round[#,0.0001]&)],InterpolationOrder->1][0];
GetFF[IV_List]:=With[{Isc=GetIsc@IV,Voc=GetVoc@IV},Last@MPPT@IV/(Isc*Voc)];


PlotIV[IV_List,voltageRange:{_,_}:{150,1500},opt:OptionsPattern[ListPlot]]:=Module[{Isc,Voc,mpp,plot1,plot2,textLoc},

Isc=Interpolation[DeleteDuplicatesBy[Reverse/@IV,First@*(Round[#,0.0001]&)],InterpolationOrder->1][0];
Voc=Interpolation[DeleteDuplicatesBy[IV,First@*(Round[#,0.0001]&)],InterpolationOrder->1][0];

mpp=MPPT[IV,"VoltageRange"->voltageRange]; (* {Impp,Vmpp,Pmpp} *)

plot1=ListLinePlot[Reverse/@IV,AxesLabel->{"V","I"},PlotRange->All];

textLoc=If[First@mpp<Isc/2 && mpp[[2]]<Voc*0.8,
	{Voc*0.8,Isc*0.8}
, (* If shaded, place on top right corner, else *)
	{Voc/2,First@mpp/2}
];
plot2=ListPlot[{Callout[{Voc,0},"Voc = "<>ToString@Round[Voc,0.1]],Callout[{0,Isc},"Isc = "<>ToString@Round[Isc,0.1]],Callout[Reverse@mpp[[;;2]],"Vmpp = "<>ToString@Round[mpp[[2]],0.1]<>", Impp = "<>ToString@Round[mpp[[1]],0.1],Left]},
Epilog->{Text["FF = "<>ToString@Round[Last@mpp/(Isc*Voc),0.01],textLoc,{0,-1}],Text["Power = "<>ToString@Round[Last@mpp,0.1],textLoc,{0,1}]},
opt,PlotStyle->{Red,PointSize[Large]},AxesLabel->{"V","I"},ImageMargins->5,GridLines->Automatic];

Return[Show[plot2,plot1]];

];


PlotPV[IV_List,voltageRange:{_,_}:{150,1500},opt:OptionsPattern[ListPlot]]:=Module[{Isc,Voc,mpp,plot1,plot2,VP,textLoc},

Isc=Interpolation[DeleteDuplicatesBy[Reverse/@IV,First@*(Round[#,0.0001]&)],InterpolationOrder->1][0];
Voc=Interpolation[DeleteDuplicatesBy[IV,First@*(Round[#,0.0001]&)],InterpolationOrder->1][0];

mpp=MPPT[IV,"VoltageRange"->voltageRange]; (* {Impp,Vmpp,Pmpp} *)

VP={#1,#1*#2}&@@@Reverse/@IV;
plot1=ListLinePlot[VP,AxesLabel->{"V","P"},PlotRange->All];

textLoc=If[mpp[[2]] < Voc*0.7,
	{Voc/2,Last@mpp/3}
,
	{mpp[[2]]*0.8,Last@mpp/3}
];
plot2=ListPlot[{Callout[{Voc,0},"Voc = "<>ToString@Round[Voc,0.1]],Callout[mpp[[2;;]],"Vmpp = "<>ToString@Round[mpp[[2]],0.1]<>",\nPmpp = "<>ToString@Round[mpp[[3]],0.1],Top]},
Epilog->{Text["FF = "<>ToString@Round[Last@mpp/(Isc*Voc),0.01],textLoc,{0,-1}],Text["Power = "<>ToString@Round[Last@mpp,0.1],textLoc,{0,1}]},
opt,PlotStyle->{Red,PointSize[Large]},AxesLabel->{"V","P"},ImageMargins->5,GridLines->Automatic,PlotRange->{{0,Voc*1.1},{0,Max@VP[[All,2]]*1.2}}];

Return[Show[plot2,plot1]];

];


(* ::Chapter:: *)
(*End*)


End[]; 

EndPackage[]; 
