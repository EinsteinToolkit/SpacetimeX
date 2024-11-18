(* ::Package:: *)

(* ADM_set_constraint.wl *)

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

<<wl/ADM_vars.wl

<<wl/ADM_rhs.wl

Module[{Mat, invMat},
  Mat = Table[gam[{ii, -cart}, {jj, -cart}] // ToValues, {ii, 1, 3}, {jj, 1, 3}];
  invMat = Inverse[Mat] /. {1 / Det[Mat] -> detinvgam};
  SetEQNDelayed[detinvgam[], 1 / Det[Mat] // Simplify];
  SetEQNDelayed[invgam[i_, j_], invMat[[i[[1]], j[[1]]]] // Simplify]
];

SetOutputFile[FileNameJoin[{Directory[], "ADM_set_constraint.hxx"}]];

$MainPrint[] :=
  Module[{},
    PrintInitializations[{Mode -> "MainOut"}, ConstraintVarlist];
    pr[];

    pr["noinline([&]() __attribute__((__flatten__, __hot__)) {"];
    pr["  grid.loop_int_device<0, 0, 0, vsize>("];
    pr["    grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {"];
    pr["    const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);"];
    pr["    const GF3D2index index2(layout2, p.I);"];
    pr["    const GF3D5index index5(layout5, p.I);"];
    pr[];

    PrintListInitializations[TmunuVarlist, "gf_", "index2"];
    PrintListInitializations[ADMVarlist, "tl_", "index5"];
    PrintListInitializations[dADMVarlist, "tl_", "index5"];
    PrintListInitializations[ddADMVarlist, "tl_", "index5"];
    pr[];

    PrintInitializations[{Mode -> "MainIn"}, TmunuVarlist];
    PrintInitializations[{Mode -> "MainIn", StorageType -> "Tile"},
                         ADMVarlist];
    PrintInitializations[{Mode -> "MainIn", StorageType -> "Tile", TensorType -> "Vect"},
                         dADMVarlist];
    PrintInitializations[{Mode -> "MainIn", StorageType -> "Tile", TensorType -> "Smat"},
                         ddADMVarlist];
    pr[];
    PrintEquations[{Mode -> "Temp"}, IntermediateVarlist];
    PrintEquations[{Mode -> "Temp"}, MatterVarlist];
    pr[];
    PrintEquations[{Mode -> "Main"}, ConstraintVarlist];
    pr[];
    pr["  });"];
    pr["});"];
  ];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetX.wl"}]];
