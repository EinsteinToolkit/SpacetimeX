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

SetOutputFile[FileNameJoin[{Directory[], "Z4cow_set_ADM.hxx"}]];

$MainPrint[] :=
  Module[{},
    PrintInitializations[{Mode -> "MainOut"}, ADMVarlist];
    pr[];

    pr["grid.loop_all_device<0, 0, 0, vsize>("];
    pr["  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {"];
    pr["  const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);"];
    pr["  const GF3D2index index2(layout2, p.I);"];
    pr[];

    PrintListInitializations[Drop[EvolVarlist, {5}], "gf_", "index2"];
    pr[];

    PrintInitializations[{Mode -> "MainIn"}, Drop[EvolVarlist, {5}]];
    pr[];
    PrintEquations[{Mode -> "Main"}, ADMVarlist];
    pr[];
    pr["});"];
  ];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetX.wl"}]];
