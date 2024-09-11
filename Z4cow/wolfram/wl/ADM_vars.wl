(* ::Package:: *)

(* ADM_vars.wl *)

(* (c) Liwei Ji, 07/2024 *)

(* ADM variables *)

ADMVarlist =
  GridTensors[
    {ADMgam[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*SuperscriptBox[\(\[Gamma]\), \((ADM)\)]\)"},
    {ADMK[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*SuperscriptBox[\(K\), \((ADM)\)]\)"},
    {ADMalpha[], PrintAs -> "\!\(\*SuperscriptBox[\(\[Alpha]\), \((ADM)\)]\)"},
    {ADMbeta[i], PrintAs -> "\!\(\*SuperscriptBox[\(\[Beta]\), \((ADM)\)]\)"}
    (*,
    {ADMdtalpha[], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \\(t\)]\)\!\(\*SuperscriptBox[\(\[Alpha]\), \((ADM)\)]\)"},
    {ADMdtbeta[i], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \\(t\)]\)\!\(\*SuperscriptBox[\(\[Beta]\), \((ADM)\)]\)"}
    *)
  ];
