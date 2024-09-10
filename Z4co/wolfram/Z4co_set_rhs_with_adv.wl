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

<<wl/Z4c_rhs_with_adv.wl

Module[{Mat, invMat},
  Mat = Table[gamt[{ii, -cart}, {jj, -cart}] // ToValues, {ii, 1, 3}, {jj, 1, 3}];
  invMat = Inverse[Mat] /. {1 / Det[Mat] -> 1}; (* since we enforced that det(gamt) = 1 *)
  (* SetEQNDelayed[detinvgamt[], 1 / Det[Mat] // Simplify]; *)
  SetEQNDelayed[invgamt[i_, j_], invMat[[i[[1]], j[[1]]]] // Simplify]
];

SetOutputFile[FileNameJoin[{Directory[], "Z4co_set_rhs_with_adv.hxx"}]];

$MainPrint[] :=
  Module[{},
    PrintInitializations[{Mode -> "MainOut"}, dtEvolVarlist];
    pr[];

    pr["noinline([&]() __attribute__((__flatten__, __hot__)) {"];
    pr["  grid.loop_int_device<0, 0, 0, vsize>("];
    pr["    grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {"];
    pr["    const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);"];
    pr["    const GF3D2index index2(layout2, p.I);"];
    pr["    const GF3D5index index5(layout5, p.I);"];
    pr[];

    PrintListInitializations[TmunuVarlist, "gf_", "index2"];
    PrintListInitializations[EvolVarlist, "tl_", "index5"];
    PrintListInitializations[dEvolVarlist, "tl_", "index5"];
    PrintListInitializations[ddEvolVarlist, "tl_", "index5"];
    pr[];

    PrintInitializations[{Mode -> "MainIn"}, TmunuVarlist];
    PrintInitializations[{Mode -> "MainIn", StorageType -> "Tile"},
                         EvolVarlist];
    PrintInitializations[{Mode -> "MainIn", StorageType -> "Tile", TensorType -> "Vect"},
                         dEvolVarlist];
    PrintInitializations[{Mode -> "MainIn", StorageType -> "Tile", TensorType -> "Smat"},
                         ddEvolVarlist];
    pr[];
    PrintEquations[{Mode -> "Temp"}, IntermediateVarlist];
    PrintEquations[{Mode -> "Temp"}, DDVarlist];
    PrintEquations[{Mode -> "Temp"}, RVarlist];
    PrintEquations[{Mode -> "Temp"}, MatterVarlist];
    pr[];
    PrintEquations[{Mode -> "Main"}, dtEvolVarlist];
    pr[];
    pr["  });"];
    pr["});"];
  ];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetX.wl"}]];
