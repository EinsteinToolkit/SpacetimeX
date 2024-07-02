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

<<wl/Z4c_rhs.wl

Module[{Mat, invMat},
  Mat = Table[gamt[{ii, -cart}, {jj, -cart}] // ToValues, {ii, 1, 3}, {jj,
     1, 3}];
  invMat = Inverse[Mat] /. {1 / Det[Mat] -> (detinvgamt[] // ToValues)};
  SetEQNDelayed[detinvgamt[], 1 / Det[Mat] // Simplify];
  SetEQNDelayed[invgamt[i_, j_], invMat[[i[[1]], j[[1]]]] // Simplify]
];

SetOutputFile[FileNameJoin[{Directory[], "Z4co_set_rhs.hxx"}]];

$MainPrint[] :=
  Module[{},
    PrintInitializations[{Mode -> "GF3D5"}, dtEvolVarlist];
    PrintInitializations[{Mode -> "GF3D5"}, ConstraintVarlist];
    PrintInitializations[{Mode -> "GF3D5"}, TmunuVarlist];
    PrintInitializations[{Mode -> "GF3D2"}, EvolVarlist];
    PrintInitializations[{Mode -> "VecGF3D2"}, dEvolVarlist];
    PrintInitializations[{Mode -> "SmatGF3D2"}, ddEvolVarlist];
    pr[];
    PrintEquations[{Mode -> "Temp"}, IntermediateVarlist];
    PrintEquations[{Mode -> "Temp"}, DDVarlist];
    PrintEquations[{Mode -> "Temp"}, RVarlist];
    PrintEquations[{Mode -> "Temp"}, MatterVarlist];
    PrintEquations[{Mode -> "Temp"}, dAtUUVarlist];
    pr[];
    PrintEquations[{Mode -> "Temp"}, ConstraintVarlist];
    PrintEquations[{Mode -> "Main"}, dtEvolVarlist];
    pr[];
  ];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetX.wl"}]];
