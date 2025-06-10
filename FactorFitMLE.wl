(* ::Package:: *)

BeginPackage["FactorFitMLE`"]

(* Define usage attributes for public functions *)

xFactorFitMLE::usage = 
  "{{FactorLoadings, ErrorMatrix}, LLHistory, BIC} = xFactorModelFitMLE[NumbObs, CovMatrix, {InitLoadings, InitErrMatrix}, opts} - Fit an orthonormal factor model\n\n"<>
  "NumbObs - Number of observations in dataset\n"<>
  "CovMatrix - Covariance matrix\n"<>
  "InitLoadings - Initialized factor loadings\n"<>
  "InitErrMatrix - Initialized diagonal matrix of error variances\n\n"<>
  "Option names are strings:\n"<>
  "  \"ToleranceGoal\" - Rate of change in log likelihood, criterion for EM termination (default 10^-8)\n"<>
  "  \"MaxIterations\" - Max number of iterations, criterion for EM termination (default 400)\n\n"<>
  "FactorLoadings - Factor loading matrix\n"<>
  "ErrorMatrix - Diagonal matrix of error variances\n"<>
  "LLHistory - History of log likelihood at each EM iteration\n"<>
  "BIC - Bayesian information criterion for the model";

Options[xFactorFitMLE] = {"ToleranceGoal" -> 10.^-8, 
   "MaxIterations" -> 400};

xInitializeFactorModel::usage =
  "{InitLoadings, InitErrMatrix} = xInitializeFactorModel[CovMatrix, Order] - Produce intial estimates for factor model\n\n"<>
  "CovMatrix - Covariance matrix\n"<>
  "Order - Number of factors in model\n\n"<>
  "InitLoadings - Initialized factor loadings\n"<>
  "InitErrMatrix - Initialized diagonal matrix of error variances";

(* Define all functions in the private context *)

Begin["`Private`"]

xInverse[mnB_, mnD_] := Module[
   {mnID},
   mnID = DiagonalMatrix[1/Tr[mnD, List]];
   mnID - 
    mnID.mnB.Inverse[
      IdentityMatrix[Last[Dimensions[mnB]]] + 
       mnB\[Transpose].mnID.mnB].mnB\[Transpose].mnID
   ];

xInitializeFactorModel[mnObsCov_, iOrder_] := Module[
   {mnErrDiag, mnFactor},
   mnErrDiag = DiagonalMatrix[Diagonal[mnObsCov]/2];
   mnFactor = (First[#].Sqrt[#[[2]]]) &[
     SingularValueDecomposition[mnObsCov - mnErrDiag, iOrder]
     ];
   {mnFactor, mnErrDiag}
   ];
xLogLik[iN_, mnOC_, mnIQ_, 
  nDQ_] := -0.5 iN (Log[nDQ] + Total[Diagonal[mnOC.mnIQ]])

xFactorFitMLE[iN_, mnObsCov_, {mnInitFactors_, mnInitErrDiag_}, 
   OptionsPattern[]] := Module[
   {nBIC, nDetCov, iF, mnErrDiag, mnFactors, mnInvCov, vnLogLikHistory, mnM,
     iMaxIter, mnProj, nTol},
   (* Resolve options and set up run parameters *)
   
   iMaxIter = OptionValue["MaxIterations"];
   nTol = OptionValue["ToleranceGoal"];
   iF = Last[Dimensions[mnInitFactors]];
   mnFactors = mnInitFactors;
   mnErrDiag = mnInitErrDiag;
   (* Pre-loop *)
   mnInvCov = xInverse[mnFactors, mnErrDiag];
   nDetCov = Det[mnFactors.mnFactors\[Transpose] + mnErrDiag];
   vnLogLikHistory = {xLogLik[iN, mnObsCov, mnInvCov, nDetCov]};
   For[iIter = 1, iIter <= iMaxIter, iIter++,
    (* E-Step *)
    mnProj = mnFactors\[Transpose].mnInvCov;
    mnM = mnObsCov.Transpose[mnProj];
    (* M-Step *)
    mnFactors = 
     mnM.Inverse[IdentityMatrix[iF] - mnProj.mnFactors + mnProj.mnM];
    mnErrDiag = 
     DiagonalMatrix[Diagonal[mnObsCov - mnFactors.Transpose[mnM]]];
    (* Log likelihood *)
    mnInvCov = xInverse[mnFactors, mnErrDiag];
    nDetCov = Det[mnFactors.mnFactors\[Transpose] + mnErrDiag];
    vnLogLikHistory = 
     Append[vnLogLikHistory, xLogLik[iN, mnObsCov, mnInvCov, nDetCov]];
    (* Termination test *)
    If[vnLogLikHistory[[-1]]/vnLogLikHistory[[-2]] - 1 <= nTol, 
     Break[]];
   ];
   nBIC = -2 Last[vnLogLikHistory] + (Times@@Dimensions[mnFactors] + Length[mnErrDiag]) Log[iN];
   {{mnFactors, mnErrDiag}, vnLogLikHistory, nBIC}
  ];

End[ ]

EndPackage[ ]






