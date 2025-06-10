(* ::Package:: *)

(* MeanCovMissingMLE - Maximum likelihood estimate of the mean and covariance with missing data. *)

(* FQS Capital Management, LLP (c) 2009, 2010, 2011, 2015 *)

(* Robert J Frey - 2009-11-15 *)
(*     01 - RJF - 2010-06-11 - Correct E-Step to ensure cross product matrix is symmetric. *)
(*                              Improve feedback on Verbose mode to show log likelihood at start of iteration. *)
(*     02 - RJF - 2010-10-30 - Extend initial option to include a user provided mean and covariance. *)
(*     03 - RJF - 2011-01-10 - Increase stability forming augmented covariance to ensure symmetry *)
(*     04 - RJF - 2011-08-13 - Add a preliminary E and M step before computing an initial log likelihood.  This *)
(*                             avoids certain numerical instabilities was computing the log likelihood from an *)
(*                             intial mean vector and covariance matrix. *)
(*     05 - RJF - 2015-07-23 - Add an X-Step (i.e., extrapolation) prior to each E-Step to accelerate convergence. *)
(*                             Fix the iteration counter so that it includes the starter EM-iteration that *)
(*                             occurs prior to the main loop. Add checkpointing for restart. *)
(*     06 - RJF - 2015-08-03 - Add log file controlled by "LogFile" boolean option." *)
(*     07 - RJF - 2015-08-04 - Improvements and standardization of verbose and logging operations. Added Run ID to *)
(*                             function return. *)
(*     08 - RJF - 2015-08-08 - Incorporate Giri's suggestions to make filenaming more standardized and simpler. *)

BeginPackage["MeanCovMissingMLE`"]

(* Define usage attributes for public functions *)

xSweep::usage = 
"Swept = xSweep[AugCov, SrcIndices] - Dempster's form of Beaton's sweep operator.\n\n"<>
"AugCov - augmented cov matrix. e.g., with means joined to bottom row and right column and -1 in the lower right corner.\n"<>
"SrcIndices - Integer vector indicating source rows and columns of AugCov.\n\n"<>
"Swept\[LeftDoubleBracket]Prepend[SrcIndices,=1], TgtIndices\[RightDoubleBracket] - Matrix of regression coefficients with intercept leading and TgtIndices the complement of the SrcIndices.\n"<>
"Swept\[LeftDoubleBracket]TgtIndices, TgtIndices\[RightDoubleBracket] - The cov matrix of the errors.";

xReverseSweep::usage =
"Unswept = xReverseSweep[Swept, SrcIndices] - Inverse sweep operator; see xSweep.";

xCovMLE::usage =
"CovMatrix = xCovMLE[DataMatrix] - MLE covariance matrix given the DataMatrix in which each row is an observation.\n"<>
"Covariance = xCovMLE[VectorX, VectorY] - MLE covariance estimate given the data vectors.";

xMeanCovMissingMLE::usage = 
"{{Mean, Cov}, {Completed, Cross}, LLHistory, RunID} = xMeanCovMissingMLE[dbReferenceReturns, opts] - MLE estimate of mean and cov with missing data.\n"<>
"    Uses XEM-algorithm (extrapolation, expectation, minimization).\n\n" <>
"dbReferenceReturns - db of reference returns; missing values indicated by Missing[].\n"<>
"Option names are strings:\n"<>
"   \"InitialEstimate\" - \"Diagonal\" (default), \"CompleteObs\", and \"ImputedMean\" or {Mean, Cov}.\n"<>
"      If a string, then defines how the intial estimate is produced; otherwise, intial starting values for the mean vector and coveariance matrix.\n"<>
"   \"ToleranceGoal\" - Rate of change in log likelihood, criterion for EM termination, default 10^-6.\n"<>
"   \"MaxIterations\" - Max number of iterations, criterion for EM termination, default 100.\n"<>
"   \"Extrapolation\" - True to use the X-Step (default); False to use standard EM algorithm.\n"<>
"   \"CheckPoint\" - Frequency at which state sufficient for restart saved; 0 (default) for no check-pointing.\n" <>
"   \"Verbose\" - True to print info on each iteration, default False.\n\n"<>
"Mean - Mean vector\n"<>
"Cov - Covariance matrix\n"<>
"Completed - Data with missing values replaced by expectations based on regression estimates.\n"<>
"Cross - Matrix for second order correction of cross products involving missing data.\n"<>
"LLHistory - History of log likelihood at each EM iteration.\n" <>
"RunID - Unique string representing the environment and time stamp of the run.";

(* HistoryLength is (for now) an undocumented feature; it defaults to 2 (linear extrapolation). *)
Options[xMeanCovMissingMLE]={
     "InitialEstimate" -> "Diagonal",
     "ToleranceGoal" -> 10.^-6,
     "MaxIterations" -> 400,
     "Verbose" -> False,
     "LogFile" -> False,
     "Extrapolation" -> True,
     "CheckPoint" -> 0,
     "HistoryLength" -> 2
};


(* Define all functions in a private context within package *)
Begin["`Private`"]


xSweep[mnM_,iK_]:=Array[
Which[
iK==#1==#2,-1/mnM[[iK,iK]],
#1!=iK&&#2==iK,mnM[[#1,iK]]/mnM[[iK,iK]],
#1==iK&&#2!=iK,mnM[[#2,iK]]/mnM[[iK,iK]],
#1!=iK&&#2!=iK,mnM[[#1,#2]]-mnM[[#1,iK]]mnM[[iK,#2]]/mnM[[iK,iK]]
]&,
Dimensions[mnM]
]/;1<=iK<=Length[mnM];


xReverseSweep[mnM_,iK_]:=Array[
Which[
iK==#1==#2,-1/mnM[[iK,iK]],
#1!=iK&&#2==iK,-mnM[[#1,iK]]/mnM[[iK,iK]],
#1==iK&&#2!=iK,-mnM[[#2,iK]]/mnM[[iK,iK]],
#1!=iK&&#2!=iK,mnM[[#1,#2]]-mnM[[#1,iK]]mnM[[iK,#2]]/mnM[[iK,iK]]
]&,
Dimensions[mnM]
]/;1<=iK<=Length[mnM];


xMeanCovMissingMLE[dbReferenceReturns_,OptionsPattern[]]:=Module[
{vqAvail, iCheckPoint, mnCov, mnData, mnErrorCov, qExtrapolation, iHistoryLength, xInitialEstimate, iIndexPad, iIter, stmLog, sLogDelimit, sLogExt, qLogFile, vnLogLikHistory, miMis, iMaxIter, 
      vnMean, vxMeanCovHistory, vsMessage, sModelDate, iN, miObs, vnPriorMean, mnPriorCov, mnPriorWorking, mnPriorErrorCov, vnPriorLogLikHistory, sRunID, nTestLogLik, nTol, qVerbose,
      sVerboseDelimit, sVerboseMargin, mnWorking, qXSuccess},
(* Extract data matrix *)
mnData = dbReferenceReturns[[3]];
(* Set the run ID based on starting time and environment *)
sModelDate = DateString[Last@First[dbReferenceReturns],{"Year","Month","Day"}];
sRunID = xSetRunID[sModelDate];
(* Off message on possible inverse stabiltiy. *)
Off[Inverse::luc];
(* Resolve options and set up run parameters *)
xInitialEstimate=OptionValue["InitialEstimate"];
iMaxIter=OptionValue["MaxIterations"];
iIndexPad = iIndexPad =Length@IntegerDigits[iMaxIter];
nTol=OptionValue["ToleranceGoal"];
iHistoryLength = Max[OptionValue["HistoryLength"], 2];
qExtrapolation = OptionValue["Extrapolation"];
iCheckPoint = OptionValue["CheckPoint"];
qVerbose=OptionValue["Verbose"];
If[qVerbose, sVerboseDelimit = " - "; sVerboseMargin = "   "];
qLogFile = OptionValue["LogFile"];
If[qLogFile, sLogDelimit=","; sLogExt = "csv"; stmLog = xOpenLogFile[sRunID, sLogExt]];
vsMessage = {"Run ID", sRunID};
If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
vsMessage = {DateString[],"Commencing mean and covariance estimation"};
If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
vsMessage = {"Logging on", xHMS[]};
If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
iN = Length[mnData];
mnWorking=mnData;
{vnMean, mnCov, vqAvail, miObs, miMis} = xInit[mnWorking, xInitialEstimate];
vxMeanCovHistory = {{vnMean, mnCov}};
vnLogLikHistory={xLogLik[mnWorking,miObs,vnMean,mnCov]};
vsMessage = {xIndexString[0, iIndexPad], ToString[Last[vnLogLikHistory]], xHMS[]};
If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
(* E-Step *)
{mnWorking,mnErrorCov}=EStep[mnWorking,vqAvail,miObs,miMis,vnMean,mnCov];
{{vnPriorMean, mnPriorCov}, {mnPriorWorking, mnPriorErrorCov}, vnPriorLogLikHistory} = {{vnMean, mnCov}, {mnWorking, mnErrorCov}, vnLogLikHistory};
vsMessage = { "E", xHMS[]};
If[qVerbose, Print[sVerboseMargin<>StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
(* M-Step *)
{vnMean, mnCov} = MStep[mnWorking, mnErrorCov];
vsMessage = {"M", xHMS[]};
If[qVerbose, Print[sVerboseMargin<>StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
(* Log likelihood *)
vnLogLikHistory = Append[vnLogLikHistory, xLogLik[mnWorking,miObs,vnMean,mnCov]];
(* Check for early termination *)
If[Not[vnLogLikHistory[[-1]] \[Element] Reals] || (vnLogLikHistory[[-1]] <= vnLogLikHistory[[-2]]) || (vnLogLikHistory[[-1]] / vnLogLikHistory[[-2]] - 1 <= nTol),
     vsMessage = {"Likelihood exception - restoring prior iterates", ToString[vnLogLikHistory[[-1]]], xHMS[]};
     If[qVerbose, Print[sVerboseMargin<>StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
     If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
     {{vnMean, mnCov}, {mnWorking, mnErrorCov}, vnLogLikHistory} = {{vnPriorMean, mnPriorCov}, {mnPriorWorking, mnPriorErrorCov}, vnPriorLogLikHistory};
     vnLogLikHistory = Append[vnLogLikHistory, vnLogLikHistory[[-1]]];  
     vsMessage = {xIndexString[1, iIndexPad], ToString[Last[vnLogLikHistory]], xHMS[]};
     If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
     If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
     vsMessage = {DateString[] ,"Finished mean and covariance estimation", sRunID};
     If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
     If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
     If[qLogFile, Close[stmLog]];
     On[Inverse::luc];
     Return[{{vnMean, mnCov}, {mnWorking, mnErrorCov}, vnLogLikHistory, sRunID}]
];
vsMessage = {xIndexString[1, iIndexPad], ToString[Last[vnLogLikHistory]], xHMS[]};
If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
If[iMaxIter == 1, 
     vsMessage = {DateString[] ,"Finished mean and covariance estimation (MaxIter = 1)", sRunID};
     If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
     If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
     If[qLogFile, Close[stmLog]]; 
     On[Inverse::luc]; 
     Return[{{vnMean, mnCov}, {mnWorking, mnErrorCov}, vnLogLikHistory, sRunID}]
];
For[iIter=2 , iIter <= iMaxIter, iIter++,
     {{vnPriorMean, mnPriorCov}, {mnPriorWorking, mnPriorErrorCov}, vnPriorLogLikHistory} = {{vnMean, mnCov}, {mnWorking, mnErrorCov}, vnLogLikHistory};
     (* X-Step - Attempt if sufficient history available. *)
     If[qExtrapolation && Length[vxMeanCovHistory]>= iHistoryLength,
          {vnMean, mnCov, qXSuccess, nTestLogLik} = XStep[mnWorking, miObs, vnMean, mnCov, vxMeanCovHistory, vnLogLikHistory];
          vsMessage = {"X", If[qXSuccess, "Accepted", "Rejected"], If[NumericQ[nTestLogLik], ToString[nTestLogLik], "NaN"], xHMS[]};
          If[qVerbose, Print[sVerboseMargin<>StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
          If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
     ];
	(* E-Step *)
	{mnWorking, mnErrorCov} = EStep[mnWorking, vqAvail, miObs, miMis, vnMean, mnCov];
    vsMessage = {"E", xHMS[]};
    If[qVerbose, Print[sVerboseMargin<>StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
    If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
	(* M-Step *)
    {vnMean, mnCov} = MStep[mnWorking, mnErrorCov];
	vsMessage = {"M", xHMS[]};
    If[qVerbose, Print[sVerboseMargin<>StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
    If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
	(* Log likelihood *)
	vnLogLikHistory = Append[vnLogLikHistory, xLogLik[mnWorking, miObs, vnMean, mnCov]];
    If[Not[vnLogLikHistory[[-1]] \[Element] Reals]||vnLogLikHistory[[-1]] < vnLogLikHistory[[-2]],
        vsMessage = {"Likelihood exception - restoring prior iterates", ToString[vnLogLikHistory[[-1]]]. xHMS[]};
        If[qVerbose, Print[sVerboseMargin<>StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
        If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
         {{vnMean, mnCov}, {mnWorking, mnErrorCov}, vnLogLikHistory} = {{vnPriorMean, mnPriorCov}, {mnPriorWorking, mnPriorErrorCov}, vnPriorLogLikHistory};
         vnLogLikHistory = Append[vnLogLikHistory, vnLogLikHistory[[-1]]];
     ];
    (* Update mean-cov history. *)
    vxMeanCovHistory = If[Length[vxMeanCovHistory] < iHistoryLength,
        Append[vxMeanCovHistory , {vnMean, mnCov}],
        Append[Rest@vxMeanCovHistory, {vnMean, mnCov}]
    ];
    (* Check Pointing *)
    If[iCheckPoint > 0 && Mod[iIter, iCheckPoint] == 0,
        vsMessage = {"C", xHMS[]};
        If[qVerbose, Print[sVerboseMargin<>StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
        If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
        xCheckPoint[iIter, sRunID];
        xCheckPoint[vnMean, sRunID];
        xCheckPoint[mnCov, sRunID];
        xCheckPoint[mnWorking, sRunID];
        xCheckPoint[vnLogLikHistory, sRunID]
    ];        
    vsMessage = {xIndexString[iIter, iIndexPad], ToString[Last[vnLogLikHistory]], xHMS[]};
    If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
    If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
	(* Termination test *)
	If[(vnLogLikHistory[[-1]] <= vnLogLikHistory[[-2]]) || (vnLogLikHistory[[-1]] / vnLogLikHistory[[-2]] - 1 <= nTol), Break[]];
];
vsMessage = {DateString[] ,"Finished mean and covariance estimation", sRunID};
If[qVerbose, Print[StringJoin@Riffle[vsMessage,sVerboseDelimit]]];
If[qLogFile, WriteLine[stmLog, StringJoin@Riffle[vsMessage,sLogDelimit]]];
If[qLogFile, Close[stmLog]];
(* On message on possible inverse stabiltiy. *)
On[Inverse::luc];
{{vnMean, mnCov}, {mnWorking, mnErrorCov}, vnLogLikHistory, sRunID}
];


xSweep[mnM_?SymmetricMatrixQ,viK_?(VectorQ[#,IntegerQ]&)]:=Fold[
xSweep[#1,#2]&,mnM,viK
]/;And@@(0<=#<=Length[mnM]&/@viK);


xReverseSweep[mnM_?SymmetricMatrixQ,viK_?(VectorQ[#,IntegerQ]&)]:=Fold[
xReverseSweep[#1,#2]&,mnM,viK
]/;And@@(0<=#<=Length[mnM]&/@viK);


xCovMLE[x_]:=x\[Transpose].x/Length[x]-KroneckerProduct[#,#]&[Mean[x]];
xCovMLE[x_,y_]:=x.y/Length[x]-Mean[x]Mean[y];


xFormAugCov[vnM_,mnC_]:=(0.5 (#+Transpose[#]))&[Join[Append[mnC,vnM],{Append[vnM,-1]}\[Transpose],2]];


xHMS ={}\[Function]DateString[{"Hour24",":","Minute",":","Second"}, "TimeZone" -> 0];


xIndexString[i_, n_] := StringTake[ToString@NumberForm[i, n, NumberPadding -> {"0",""}], -n];


xSetRunID[sModelDate_String]:=
    $UserName<>"_"<>$SystemID <>"_"<>ToString[$KernelID]<>"_"<> sModelDate <> "_" <>
    DateString[{"Year","Month","Day","Hour24","Minute","Second"}, "TimeZone" -> 0];


SetAttributes[xCheckPoint, HoldFirst];

xCheckPoint[xVar_Symbol, sID_String] := Module[
{sDir},
sDir=FileNameJoin[{NotebookDirectory[], "CheckPoint_" <> sID}];
If[!FileExistsQ[sDir], CreateDirectory[sDir]];
Export[
     FileNameJoin[{
          sDir,
          (* Strip off contexts and post "$" sequence; e.g., Context`name$123 \[Rule] name *)
          First[StringSplit[Last[StringSplit[ToString[HoldForm[xVar]], "`"]], "$"]] <> ".m"
     }],
     xVar
]
];


xOpenLogFile[sID_String, sLogExt_] := Module[
{sFile},
sFile = FileNameJoin[{NotebookDirectory[], "MeanCovMissingMLE_Log_" <> sID <> "." <> sLogExt}];
OpenWrite[sFile]
];


xInit[mnX_,xIE_]:=Module[
{mqA,vqA,mnC,iD,iI,iJ,miM,vnM,iN,miO,vqS},
(* drop all-missing obs, get dims and mark missing data in remaining *)
mqA=Map[(#/.{_Missing->False,_->True})&,mnX,{2}];
{iN,iD}=Dimensions[mnX];
(* apply appropriate InitalEstimate *)
If[
	StringQ[xIE],
	(* resolve string option *)
	Which[
	(* use available data for each mean and covariance element *)
	xIE=="AvailableObs",
	vnM=Array[Missing[]&,iD];
	mnC=Array[Missing[]&,{iD,iD}];
	For[iI=1,iI<=iD,iI++,
	If[Count[mqA[[All,iI]],True]>0,
	vnM[[iI]]=Mean[Pick[mnX[[All,iI]],mqA[[All,iI]]]];
	];
	For[iJ =1,iJ<=iI,iJ++,
	vqS=And@@#&/@mqA[[All,{iI,iJ}]];
	mnC[[iJ,iI]]=mnC[[iI,iJ]]=xCovMLE@@Transpose[Pick[mnX[[All,{iI,iJ}]],vqS]];
	];
	],
	(* use only the subset of complete observations *)
	xIE=="CompleteObs",
	{vnM,mnC}=Through[{Mean,xCovMLE}[Pick[mnX,And@@#&/@mqA]]],
	(* use available data on means and variances with a diagonal covariance matrix *)
    xIE=="Diagonal",
	{vnM,mnC}=Through[{Mean/@#&,DiagonalMatrix[Mean[#^2]-Mean[#]^2&/@#]&}[Pick[mnX\[Transpose],mqA\[Transpose]]]],
	(* replace missing data with the mean *)
	xIE=="ImputedMean",
	{vnM,mnC}=Through[{Mean,xCovMLE}[
	MapThread[#1/.{_Missing->Mean[Pick[#1,#2]]}&,{mnX\[Transpose],mqA\[Transpose]}]\[Transpose]
	]],
	(* abort if no match on InitialEstimate string *)
	True,
	Abort[]
	],
	(* else use user provided estimate *)
	{vnM, mnC} =xIE
];
vqA=And@@#&/@mqA;
miO=Pick[Range[iD],#]&/@mqA;
miM=Pick[Range[iD],Not/@#]&/@mqA;
(* return mean vector and covariance matrix *)
{vnM,mnC,vqA,miO,miM}
];


xLogLik[mnX_,miO_,vnM_,mnC_]:=Module[
{mnCo,iD,iI,viI,nK,nL,vnMo,iN},
{iN,iD}=Dimensions[mnX];
viI=Range[iD];
nL=0.;
For[iI=1,iI<=iN,iI++,
vnMo=vnM[[miO[[iI]]]];
mnCo=mnC[[miO[[iI]],miO[[iI]]]];
nK=Log[2.\[Pi]];
nL+=-0.5(Length[iN]nK+Log[Det[mnCo]]+(#.Inverse[mnCo].#&[mnX[[iI,miO[[iI]]]]-vnMo]));
];
nL
];


XStep[mnWorking_, miObs_, vnMean_, mnCov_, vxMeanCovHistory_, vnLogLikHistory_] := Module[
{iInterpolationOrder, vnNewMean, mnNewCov, vxNewHistory, nTestLogLik, qXSuccess},
(* Extrapolation *)
vnNewMean = (Last[#]-First[#])+Last[#]&[vxMeanCovHistory[[-2;;,1]]];
mnNewCov = (Last[#]-First[#])+Last[#]&[vxMeanCovHistory[[-2;;,2]]];
(* Compute the log likelihood of the candidate extrapolation *)
nTestLogLik = xLogLik[mnWorking, miObs, vnNewMean, mnNewCov];
(* Test if extrapolation has resulted in a valid increase in likelihood *)
If[Not[nTestLogLik \[Element] Reals && nTestLogLik > vnLogLikHistory[[-1]]],
     (* If likelihood is not both valid and an improvement, then undo the extrapolation *)
     vnNewMean = vnMean;
     mnNewCov = mnCov;
     qXSuccess = False,
     (* else retain extrapolation *)
      qXSuccess = True
];
{vnNewMean, mnNewCov, qXSuccess, nTestLogLik}
];


MStep[mnX_,mnV_]:={Mean[mnX],xCovMLE[mnX]+mnV};


EStep[mnX_,vqA_,miO_,miM_,vnM_,mnC_]:=Module[
{mnAC,iD,iI,viI,viM,iN,viO,mnSw,xSwp,mnVsum,mnXnew},
{iN,iD}=Dimensions[mnX];
viI=Range[iD];
mnAC=xFormAugCov[vnM,mnC];
mnXnew=ConstantArray[0.,Dimensions[mnX]];
mnVsum=ConstantArray[0.,Dimensions[mnC]];
(* A hash is set-up for a given signature so that it doesn't need to be recomputed if it recurs. *)
xSwp[x_]:=xSwp[x]=xSweep[mnAC,x];
For[iI=1,iI<=iN,iI++,
If[vqA[[iI]],mnXnew[[iI]]=mnX[[iI]],
mnSw=xSwp[miO[[iI]]];
If[Not[MatrixQ[mnSw,NumericQ]], Abort[]];
mnVsum[[miM[[iI]],miM[[iI]]]]+=mnSw[[miM[[iI]],miM[[iI]]]];
mnXnew[[iI,miO[[iI]]]]=mnX[[iI,miO[[iI]]]];
mnXnew[[iI,miM[[iI]]]]=mnSw[[Prepend[miO[[iI]],-1],miM[[iI]]]]\[Transpose].Prepend[mnX[[iI,miO[[iI]]]],1]];
];
Remove[xSwp];
{mnXnew,(mnVsum+Transpose[mnVsum])/(2. iN)}
];


(* End Private` context *)
End[]


(* End package *)
EndPackage[]
