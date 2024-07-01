(* ::Package:: *)

(* Z4c_vars.wl *)

(* (c) Liwei Ji, 07/2024 *)

(***********************)

(* Evolution Variables *)

(***********************)

DefTensors[{chi[], PrintAs -> "\[Chi]"}, {gamt[-i, -j], Symmetric[{-i, 
  -j}], PrintAs -> "\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"}, {exKh[],
   PrintAs -> "\!\(\*OverscriptBox[\(K\), \(^\)]\)"}, {exAt[-i, -j], Symmetric[
  {-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(A\), \(~\)]\)"}, {trGt[i],
   PrintAs -> "\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"}, {Theta[
  ], PrintAs -> "\[CapitalTheta]"}, {alpha[], PrintAs -> "\[Alpha]"}, {beta[
  i], PrintAs -> "\[Beta]"}];

(*************************)

(* Other Input Variables *)

(*************************)

DefTensors[{dchi[-k], PrintAs -> "\[PartialD]\[Chi]"}, {dgamt[-k, -i, -
  j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"
  }, {dexKh[-k], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(K\), \(^\)]\)"
  }, {dexAt[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(A\), \(~\)]\)"
  }, {dtrGt[-k, i], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"
  }, {dTheta[-k], PrintAs -> "\[PartialD]\[CapitalTheta]"}, {dalpha[-k], 
  PrintAs -> "\[PartialD]\[Alpha]"}, {dbeta[-k, i], PrintAs -> "\[PartialD]\[Beta]"
  }];

DefTensors[{ddchi[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[PartialD]\[Chi]"
  }, {ddgamt[-l, -m, -i, -j], GenSet[Cycles[{1, 2}], Cycles[{3, 4}]], PrintAs
   -> "\[PartialD]\[PartialD]\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"},
   {ddalpha[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[PartialD]\[Alpha]"
  }, {ddbeta[-i, -j, k], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[PartialD]\[Beta]"
  }];

DefTensors[{eTtt[], PrintAs -> "\!\(\*SubscriptBox[\(T\), \(tt\)]\)"}, 
  {eTt[-i], PrintAs -> "\!\(\*SubscriptBox[\(T\), \(t\)]\)"}, {eT[-i, -j],
   PrintAs -> "T"}];

(**************************)

(* Intermediate Variables *)

(**************************)

DefTensors[{invgamt[i, j], Symmetric[{i, j}], PrintAs -> "\!\(\*OverscriptBox[\(\[Gamma]\), \(~\)]\)"
  }, {invgam[i, j], Symmetric[{i, j}], PrintAs -> "\[Gamma]"}, {gam[-i, -
  j], Symmetric[{-i, -j}], PrintAs -> "\[Gamma]"}, {GtDDD[-k, -i, -j], Symmetric[
  {-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"
  }, {GtDDU[-k, -i, j], PrintAs -> "\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"
  }, {Gt[k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(~\)]\)"
  }, {trGtd[i], PrintAs -> "(\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalGamma]\), \(~\)], \(d\)]\))"
  }, {exAtUU[i, j], Symmetric[{i, j}], PrintAs -> "\!\(\*OverscriptBox[\(A\), \(~\)]\)"
  }];

DefTensors[{tDtDchi[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(D\), \(~\)]\)\!\(\*OverscriptBox[\(D\), \(~\)]\)\[Chi]"
  }, {DDalpha[-i, -j], Symmetric[{-i, -j}], PrintAs -> "DD\[Alpha]"}];

DefTensors[{Rtchi[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*SuperscriptBox[OverscriptBox[\(R\), \(~\)], \(\[Chi]\)]\)"
  }, {Rt[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\!\(\*OverscriptBox[\(R\), \(~\)]\)"
  }, {R[-i, -j], Symmetric[{-i, -j}]}, {trR[], PrintAs -> "R"}];

(* Matter *)

DefTensors[{rho[], PrintAs -> "\[Rho]"}, {Sm[-i], PrintAs -> "S"}, {Ss[
  -i, -j], Symmetric[{-i, -j}], PrintAs -> "S"}, {trSs[], PrintAs -> "S"}
  ];

(* Constraints *)

DefTensors[{dexAtUU[-k, i, j], Symmetric[{i, j}], PrintAs -> "\[PartialD]\!\(\*OverscriptBox[\(A\), \(~\)]\)"
  }];

DefTensors[{ZtC[i], PrintAs -> "\!\(\*OverscriptBox[\(Z\), \(~\)]\)"}, 
  {HC[], PrintAs -> "H"}, {MtC[i], PrintAs -> "\!\(\*OverscriptBox[\(M\), \(~\)]\)"
  }];

(************)

(* Constant *)

(************)

DefConstantSymbol[cpi, PrintAs -> "\[Pi]"];

DefConstantSymbol[ckappa1, PrintAs -> "\!\(\*SubscriptBox[\(\[Kappa]\), \(1\)]\)"
  ];

DefConstantSymbol[ckappa2, PrintAs -> "\!\(\*SubscriptBox[\(\[Kappa]\), \(2\)]\)"
  ];

DefConstantSymbol[cmuL, PrintAs -> "\!\(\*SubscriptBox[\(\[Mu]\), \(L\)]\)"
  ];

DefConstantSymbol[cmuS, PrintAs -> "\!\(\*SubscriptBox[\(\[Mu]\), \(S\)]\)"
  ];

DefConstantSymbol[ceta, PrintAs -> "\[Eta]"];
