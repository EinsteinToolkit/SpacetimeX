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

SetOutputFile[FileNameJoin[{Directory[], "Z4c_set_rhs.hxx"}]];

$MainPrint[] :=
  Module[{},
    PrintInitializations[{Mode -> "MainIn"}, EvolVarlist];
    PrintInitializations[{Mode -> "MainIn"}, dEvolVarlist];
    PrintInitializations[{Mode -> "MainIn"}, ddEvolVarlist];
  ];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetX.wl"}]];
