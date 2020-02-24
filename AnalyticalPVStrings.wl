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
MPPT::usage = "extract maximum power point in the IV curve."]

Options[MPPT]={"TrackingMethod"->"Fast","VoltageRange"->{150,1500},"SearchStep"->5,PowerLimit->None};

If[ Not@ValueQ[CablingCorrection::usage],
CablingCorrection::usage = "Modifies IV to include effect of cabling series resistance. 
Inputs are [IV curve, cable length (round trip), cable cross section (optional), resistivity (optional)]."]




If[ Not@ValueQ[ShadeAreaFraction::usage],
ShadeAreaFraction::usage = "ShadeAreaFraction[pitch, table length, collector width, tilt, orientation, sun elevation, sun azimuth, table position (defulat 5), table number (default 10)] 
performs simple estimation of area based inter-row shading fraction of a typical array (sheds)."]

If[ Not@ValueQ[ElectricalShading::usage],
ElectricalShading::usage = "ElectricalShading[areaShaded, numSubString, orientation] gives estimation on electrical shading by giving."]

Options[ElectricalShading]={"ModuleDimension"->{6,12}}; (* default standard module dimension is 6x12 cells *)


If[ Not@ValueQ[IAMashrae::usage],
IAMashrae::usage = "IAM[AOI,b0_0.05] ashrae model."]


If[ Not@ValueQ[PlotIV::usage],
PlotIV::usage = "plot the IV curve with key points labelled."]



(* ::Chapter:: *)
(*Main definitions*)


Begin["`Private`"];

Off[InterpolatingFunction::dmval];


(* ::Section:: *)
(*IV curve functions*)


(* ::Text:: *)
(*Overall convention: *)
(*IV curve functions require a set of IV curves in the format of {{current 1, voltage 1}, {current 2, voltage 2}...}. *)
(*The IV curves must be sorted by current in ascending order (from small to large, towards Isc, i.e. voltage from Voc, from +ve to -ve). *)
(*Range of the input IV curves should be complete (at least covering one whole quadrant and cross the two axes). *)


(* ::Text:: *)
(*StringIV takes care of combining IV in series both on string level or sub-module level. *)


StringIV[IVset_,opt:OptionsPattern[]]:=Module[{maxJ,bypass=OptionValue["BypassDiode"],diodeVoltage,combinedIV,IV$fleetInterp,probe},
If[bypass==True,maxJ=Max@IVset[[All,-1,1]],maxJ=Min@IVset[[All,-1,1]]]; (*when bypass diodes are present, scan probing current to the max available in the IV set*)
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

While[(combinedIV[[-1,2]]>-diodeVoltage*Length@IVset*1.1 && probe<maxJ*1.1)||combinedIV[[-1,2]]>0,
(* as long as not all bypass diodes are activated in the IVset, i.e. voltage is not negatively biased enough *)
	probe=combinedIV[[-1,1]]*1.01;
	AppendTo[combinedIV,{probe,Total@Table[
		With[{interpPt=IV$fleetInterp[[i]][probe]},If[bypass==True&&interpPt<-diodeVoltage,-diodeVoltage,interpPt]]
		,{i,Length@IVset}]}
	];
];

Return[combinedIV];
];


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


MPPT[IV_,opt:OptionsPattern[]]:=Module[{method=OptionValue["TrackingMethod"],vRange=OptionValue["VoltageRange"],searchStep=OptionValue["SearchStep"],pLimit=OptionValue[PowerLimit],IVP,IV2,interp,maxJ,MPP={0,0,0},x,MPPsearch},

If[method=="Fast",
	IVP=Append[#,Times@@#]&/@IV;
	IVP=Select[IVP,#[[3]]>=0&&vRange[[1]]<#[[2]]<vRange[[2]]&];
	If[Length@IVP!=0,
		MPP=First[IVP~Extract~Position[#,Max@#]&@Part[IVP\[Transpose],3]];
		(*MPP=First@MaximalBy[IVP,Last];*)
		
		(* power limit *)
		If[NumericQ@pLimit&&pLimit<Last@MPP,
			interp=Interpolation[DeleteDuplicatesBy[Round[Reverse/@IV,0.0001],First],InterpolationOrder->1]; (* interpolation: x\[Rule]voltage, y\[Rule]current *)
			MPPsearch=Pick[MovingAverage[IV[[All,2]],2],Negative/@Times@@@Partition[Times@@@IV-pLimit,2,1]]; (* pick out the voltages near pLimit *)
			If[Length@MPPsearch==0,
				MPP={0,0,0};
			,
				MPPsearch=Select[x/.FindRoot[interp[x]*x==pLimit,{x,MPPsearch}],vRange[[1]]<=#<=vRange[[2]]&];
				If[Length@MPPsearch==0,
				MPP={0,0,0};
				,
				MPPsearch=Max@MPPsearch; (* choose the larger voltage as the operating point *)
				MPP=With[{i=interp[MPPsearch]},{i,MPPsearch,i*MPPsearch}];
				];
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


CablingCorrection[IV_,cableLength_,crossSection_:6,\[Rho]_:0.023]:=Module[{\[Delta]V},
\[Delta]V[current_]:=current*\[Rho]*cableLength/crossSection;

(* \[Delta]V is negative when current becomes negative *)
Return[{#1,#2-\[Delta]V[#1]}&@@@IV]; 
];


(* ::Section::Closed:: *)
(*Shading logics*)


(* ::Text:: *)
(*ShadeAreaFraction performs simple estimation of area based inter-row shading fraction of a typical array (sheds). *)
(*Determines shading for direct beam component only. *)


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
Return[{shadeFraction,rowFraction,widthFraction}];

];


(* ::Text:: *)
(*ElectricalShading uses a set of simple logic to estimate the amount of complete and partial shading of substrings in terms of rows of substrings for landscape (each PV module has three substring if there are three bypass diodes) or rows of cells for portrait. *)
(*Applicable to standard c-Si modules with bypass diodes, in the case of horizontally symmetric shading mode. *)
(*Requires shaded area fraction of the table, and number of substrings in the transverse direction of the table. *)
(*In case there are more than one string per table, i.e. upper and lower sub-row belongs to different PV strings, additional translation is needed separately for each upper and lower PV string. *)


ElectricalShading[areaShaded_,numSubString_,orientation_:"landscape",opt:OptionsPattern[]]:=Module[{subStringArea,eShadeElementFn,modDimension=OptionValue["ModuleDimension"],nCompleteShade,remain,fractionalShade,x},
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

];


(* ::Text:: *)
(*ElectricalShading in irregular shading mode... *)


(* ::Section::Closed:: *)
(*Incidence angle modifier*)


IAMashrae[AOI_,b0_:0.05]:=Max[1-b0*(1/Cos[AOI \[Degree]]-1),0]/;0<=AOI<90;
IAMashrae[AOI_?(!0<=#<90&),b0_:0.05]:=0;


(* ::Section:: *)
(*Auxiliary functions*)


PlotIV[IV_List,opt:OptionsPattern[]]:=Module[{Isc,Voc,mpp,plot1,plot2},

Isc=Interpolation[DeleteDuplicatesBy[Reverse/@IV,First],InterpolationOrder->1][0];
Voc=Interpolation[DeleteDuplicatesBy[IV,First],InterpolationOrder->1][0];

mpp=MPPT@IV;

plot1=ListLinePlot[Reverse/@IV,AxesLabel->{"V","I"},PlotRange->All];
plot2=ListPlot[{Callout[{Voc,0},"Voc = "<>ToString@Round[Voc,0.1]],Callout[{0,Isc},"Isc = "<>ToString@Round[Isc,0.1]],Callout[Reverse@mpp[[;;2]],"Vmpp = "<>ToString@Round[mpp[[2]],0.1]<>", Impp = "<>ToString@Round[mpp[[1]],0.1],Left]},
PlotStyle->{Red,PointSize[Large]},Epilog->{Text["FF = "<>ToString@Round[Last@mpp/(Isc*Voc),0.01],{Voc/2,Isc/2},{0,-1}],Text["Power = "<>ToString@Round[Last@mpp,0.1],{Voc/2,Isc/2},{0,1}]},opt];

Return[Show[plot2,plot1]];

];


(* ::Chapter:: *)
(*End*)


End[]; 

EndPackage[]; 
