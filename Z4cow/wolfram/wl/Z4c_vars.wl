(* ::Package:: *)

(* Z4c_vars.wl *)

(* (c) Liwei Ji, 07/2024 *)

(***********************)

(* Evolution Variables *)

(***********************)

dtEvolVarlist =
  GridTensors[
    {dtW[], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)W"},
    {dtgamt[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"},
    {dtexKh[], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)\!\(\*OverscriptBox[\(K\), \(^\)]\)"},
    {dtexAt[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)\!\(\*OverscriptBox[\(A\), \(~\)]\)"},
    {dttrGt[i], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"},
    {dtTheta[], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)\[CapitalTheta]"},
    {dtalpha[], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)\[Alpha]"},
    {dtbeta[i], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)\[Beta]"}
  ];

EvolVarlist =
  GridTensors[
    {W[], PrintAs -> "W"},
    {gamt[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"},
    {exKh[], PrintAs -> "\!\(\*OverscriptBox[\(K\), \(^\)]\)"},
    {exAt[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(A\), \(~\)]\)"},
    {trGt[i], PrintAs -> "\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"},
    {Theta[], PrintAs -> "\[CapitalTheta]"}, {alpha[], PrintAs -> "\[Alpha]"},
    {beta[i], PrintAs -> "\[Beta]"}
  ];

(*************************)

(* Other Input Variables *)

(*************************)

dEvolVarlist =
  GridTensors[
    {dW[-k], PrintAs -> "\[PartialD]W"},
    {dgamt[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"},
    {dexKh[-k], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(K\), \(^\)]\)"},
    {dexAt[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(A\), \(~\)]\)"},
    {dtrGt[-k, i], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"},
    {dTheta[-k], PrintAs -> "\[PartialD]\[CapitalTheta]"},
    {dalpha[-k], PrintAs -> "\[PartialD]\[Alpha]"},
    {dbeta[-k, i], PrintAs -> "\[PartialD]\[Beta]"}
  ];

ddEvolVarlist =
  GridTensors[
    {ddW[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[PartialD]W"},
    {ddgamt[-l, -m, -i, -j], GenSet[Cycles[{1, 2}], Cycles[{3, 4}]], PrintAs -> "\[PartialD]\[PartialD]\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"},
    {ddalpha[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[PartialD]\[Alpha]"},
    {ddbeta[-i, -j, k], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[PartialD]\[Beta]"}
  ];

TmunuVarlist =
  GridTensors[
    {eTtt[], PrintAs -> "\!\(\*SubscriptBox[\(T\), \(tt\)]\)"},
    {eTt[-i], PrintAs -> "\!\(\*SubscriptBox[\(T\), \(t\)]\)"},
    {eT[-i, -j], Symmetric[{-i, -j}], PrintAs -> "T"}
  ];

(**************************)

(* Intermediate Variables *)

(**************************)

IntermediateVarlist =
  TempTensors[
    {dlnW[-k], PrintAs -> "\[PartialD]lnW"},
    (* {detinvgamt[], PrintAs -> "1/\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"}, *)
    {invgamt[i, j], Symmetric[{i, j}], PrintAs -> "\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"},
    {invgam[i, j], Symmetric[{i, j}], PrintAs -> "\[Gamma]"},
    {gam[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[Gamma]"},
    {GtDDD[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"},
    {GtDDU[-k, -i, j], PrintAs -> "\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"},
    {Gt[k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"},
    {trGtd[i], PrintAs -> "(\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalGamma]\), \(~\)], \(d\)]\))"},
    {dgam[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[Gamma]"},
    {GamDDD[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[CapitalGamma]"},
    {Gam[k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[CapitalGamma]"},
    {exAtUU[i, j], Symmetric[{i, j}], PrintAs -> "\!\(\*OverscriptBox[\(A\), \(~\)]\)"}
  ];

DDVarlist =
  TempTensors[
    {tDtDW[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(D\), \(~\)]\)\!\(\*OverscriptBox[\(D\), \(~\)]\)W"},
    {DDalpha[-i, -j], Symmetric[{-i, -j}], PrintAs -> "DD\[Alpha]"}
  ];

RVarlist =
  TempTensors[
    {RtW[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*SuperscriptBox[OverscriptBox[\(R\), \(~\)], \(W\)]\)"},
    {Rt[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(R\), \(~\)]\)"},
    {R[-i, -j], Symmetric[{-i, -j}]},
    {trR[], PrintAs -> "R"}
  ];

(* Matter *)

MatterVarlist =
  TempTensors[
    {rho[], PrintAs -> "\[Rho]"},
    {Sm[-i], PrintAs -> "S"},
    {Ss[-i, -j], Symmetric[{-i, -j}], PrintAs -> "S"},
    {trSs[], PrintAs -> "S"}
  ];

(* Constraints *)

dAtUUVarlist =
  TempTensors[
    {trdexAtUU[i], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(A\), \(~\)]\)"}
  ];

ConstraintVarlist =
  GridTensors[
    {ZtC[i], PrintAs -> "\!\(\*OverscriptBox[\(Z\), \(~\)]\)"},
    {HC[], PrintAs -> "H"},
    {MtC[i], PrintAs -> "\!\(\*OverscriptBox[\(M\), \(~\)]\)"}
  ];

(************)

(* Constant *)

(************)

DefConstantSymbol[cpi, PrintAs -> "\[Pi]"];

DefConstantSymbol[ckappa1, PrintAs -> "\!\(\*SubscriptBox[\(\[Kappa]\), \(1\)]\)"];

DefConstantSymbol[ckappa2, PrintAs -> "\!\(\*SubscriptBox[\(\[Kappa]\), \(2\)]\)"];

DefConstantSymbol[cmuL, PrintAs -> "\!\(\*SubscriptBox[\(\[Mu]\), \(L\)]\)"];

DefConstantSymbol[cmuS, PrintAs -> "\!\(\*SubscriptBox[\(\[Mu]\), \(S\)]\)"];

DefConstantSymbol[ceta, PrintAs -> "\[Eta]"];
