(* ::Package:: *)

(* ADM_vars.wl *)

(* (c) Liwei Ji, 07/2024 *)

(* ADM variables *)

(*******************)

(* Input Variables *)

(*******************)

ADMVarlist =
  GridTensors[
    {gam[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[Gamma]"},
    {exK[-i, -j], Symmetric[{-i, -j}], PrintAs -> "K"},
    {alpha[], PrintAs -> "\[Alpha]"},
    {beta[i], PrintAs -> "\[Beta]"}
  ];

dADMVarlist =
  GridTensors[
    {dgam[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[Gamma]"},
    {dexK[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]K"}
  ];

ddADMVarlist =
  GridTensors[
    {ddgam[-k, -l, -i, -j], GenSet[Cycles[{1, 2}], Cycles[{3, 4}]], PrintAs -> "\[PartialD]\[PartialD]\[Gamma]"}
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
    {tr1dGam[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((tr1)\)]\)"},
    {tr2dGam[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((tr2)\)]\)"},
    {R[-i, -j], Symmetric[{-i, -j}], PrintAs -> "R"},
    {trexK[], PrintAs -> "K"},
    {DexK[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "DK"}
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
    {MC[i], PrintAs -> "M"}
  ];

(************)

(* Constant *)

(************)

DefConstantSymbol[cpi, PrintAs -> "\[Pi]"];
