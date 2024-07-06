(* ::Package:: *)

(* Z4c_set_rhs.wl *)

(* (c) Liwei Ji, 07/2024 *)

Needs["xAct`xCoba`", FileNameJoin[{Environment["GENERATO"], "src/Generato.wl"
  }]]

SetPVerbose[False];

SetPrintDate[False];

SetGridPointIndex[""];

SetTempVariableType["vreal"];

DefManifold[M3, 3, IndexRange[a, z]];

DefChart[cart, M3, {1, 2, 3}, {X[], Y[], Z[]}, ChartColor -> Blue];

(* Define Variables *)

<<wl/Z4c_vars.wl

<<wl/ADM_vars.wl

<<wl/ADM_rhs.wl

SetOutputFile[FileNameJoin[{Directory[], "Z4co_set_ADM.hxx"}]];

$MainPrint[] :=
  Module[{},
    (*PrintInitializations[{Mode -> "GF3D2Out"}, ADMVarlist];*)
    PrintInitializations[{Mode -> "GF3D2In"}, Drop[EvolVarlist, {5}]];
    pr[];
    PrintEquations[{Mode -> "MainCarpetX"}, ADMVarlist];
  ];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetX.wl"}]];
