(* ::Package:: *)

(* ADM_vars.wl *)

(* (c) Liwei Ji, 07/2024 *)

(* ADM variables *)

(*******************)

(* Input Variables *)

(*******************)

ADMVarlist =
  GridTensors[
    {ADMgam[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*SuperscriptBox[\(\[Gamma]\), \((ADM)\)]\)"},
    {ADMK[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*SuperscriptBox[\(K\), \((ADM)\)]\)"},
    {ADMalpha[], PrintAs -> "\!\(\*SuperscriptBox[\(\[Alpha]\), \((ADM)\)]\)"},
    {ADMbeta[i], PrintAs -> "\!\(\*SuperscriptBox[\(\[Beta]\), \((ADM)\)]\)"}
  ];

dADMVarlist =
  GridTensors[
    {dADMgam[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\!\(\*SuperscriptBox[\(\[Gamma]\), \((ADM)\)]\)"},
    {dADMK[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\!\(\*SuperscriptBox[\(K\), \((ADM)\)]\)"}
  ];

ddADMVarlist =
  GridTensors[
    {ddADMgam[-k, -l, -i, -j], GenSet[Cycles[{1, 2}], Cycles[{3, 4}]], PrintAs -> "\[PartialD]\!\(\*SuperscriptBox[\(\[Gamma]\), \((ADM)\)]\)"}
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
    {detinvgam[], PrintAs -> "1/\[Gamma]"},
    {invgam[i, j], Symmetric[{i, j}], PrintAs -> "\[Gamma]"},
    {GamDDD [-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[CapitalGamma]"},
    {Gam     [k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[CapitalGamma]"},
    {dGam[-m, k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[CapitalGamma]"},
    {R[-i, -j], Symmetric[{-i, -j}]},
    {trK[], PrintAs -> "K"},
    {DADMK[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[Del]K"}
  ];

(* Matter *)

MatterVarlist =
  TempTensors[
    {rho[], PrintAs -> "\[Rho]"},
    {Sm[-i], PrintAs -> "S"}
  ];

(***************)

(* Constraints *)

(***************)

ConstraintVarlist =
  GridTensors[
    {HC[], PrintAs -> "H"},
    {MtC[i], PrintAs -> "M"}
  ];

(************)

(* Constant *)

(************)

DefConstantSymbol[cpi, PrintAs -> "\[Pi]"];
