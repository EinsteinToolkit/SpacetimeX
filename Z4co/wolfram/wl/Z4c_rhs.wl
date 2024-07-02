(* ::Package:: *)

(* Z4c_rhs.wl *)

(* (c) Liwei Ji, 07/2024 *)

(* See arXiv:1212.2901 [gr-qc] *)

(****************)

(* Intermediate *)

(****************)

SetEQN[invgam[i_, j_], chi[] invgamt[i, j]];

SetEQN[gam[i_, j_], chi[] ^ -1 gamt[i, j]];

SetEQN[GtDDD[k_, i_, j_], 1/2 (dgamt[i, j, k] + dgamt[j, k, i] - dgamt[
  k, i, j])];

SetEQN[GtDDU[i_, j_, k_], invgamt[k, l] GtDDD[i, j, -l]];

SetEQN[Gt[k_, i_, j_], invgamt[k, l] GtDDD[-l, i, j]];

SetEQN[trGtd[i_], invgamt[k, l] Gt[i, -k, -l]];

SetEQN[exAtUU[i_, j_], invgamt[i, k] invgamt[j, l] exAt[-k, -l]];

SetEQN[tDtDchi[i_, j_], ddchi[i, j] - Gt[k, i, j] dchi[-k]];

SetEQN[DDalpha[i_, j_], ddalpha[i, j] - Gt[k, i, j] dalpha[-k] + 1/2 chi[
  ] ^ -1 (dchi[i] dalpha[j] + dchi[j] dalpha[i]) - 1/2 gamt[i, j] invgamt[
  k, l] chi[] ^ -1 dchi[-l] dalpha[-k]];

(* (8) *)

SetEQN[Rtchi[i_, j_], 1 / (2 chi[]) tDtDchi[i, j] + 1 / (2 chi[]) gamt[
  i, j] invgamt[k, l] tDtDchi[-k, -l] - 1 / (4 chi[] ^ 2) dchi[i] dchi[j]
   - 3 / (4 chi[] ^ 2) gamt[i, j] invgamt[k, l] dchi[-k] dchi[-l]];

(* (9) *)

SetEQN[Rt[i_, j_], -(1/2) invgamt[l, m] ddgamt[-l, -m, i, j] + 1/2 (gamt[
  -k, i] dtrGt[j, k] + gamt[-k, j] dtrGt[i, k]) + 1/2 trGtd[k] (GtDDD[i, 
  j, -k] + GtDDD[j, i, -k]) + ((Gt[k, -l, i] GtDDU[j, -k, l] + Gt[k, -l, 
  j] GtDDU[i, -k, l]) + Gt[k, i, -m] GtDDU[-k, j, m])];

(* (10) *)

SetEQN[R[i_, j_], Rtchi[i, j] + Rt[i, j]];

SetEQN[trR[], invgam[k, l] R[-k, -l]];

SetEQN[rho[], alpha[] ^ -2 (eTtt[] - 2 beta[j] eTt[-j] + beta[i] beta[j
  ] eT[-i, -j])];

SetEQN[Sm[i_], -alpha[] ^ -1 (eTt[i] - beta[k] eT[-k, i])];

SetEQN[Ss[i_, j_], eT[i, j]];

SetEQN[trSs[], invgam[k, l] Ss[-k, -l]];

SetEQN[trdexAtUU[i_], -2 invgamt[i, l] exAtUU[j, m] dgamt[-j, -l, -m] +
   invgamt[i, l] invgamt[j, m] dexAt[-j, -l, -m]];

(* (13) *)

SetEQN[ZtC[i_], trGt[i] - trGtd[i]];

(* (14) *)

SetEQN[HC[], trR[] + exAt[-k, -l] exAtUU[k, l] - 2/3 (exKh[] + 2 Theta[
  ]) ^ 2 - 16 cpi rho[]];

(* (15) *)

SetEQN[MtC[i_], trdexAtUU[i] + Gt[i, -j, -k] exAtUU[j, k] - 2/3 invgamt[
  i, j] (dexKh[-j] + 2 dTheta[-j]) - 2/3 exAtUU[i, j] chi[] ^ -1 dchi[-j]
   - 8 cpi invgamt[i, j] Sm[-j]];

(*******)

(* EOM *)

(*******)

(* (1) *)

SetEQN[dtchi[], 2/3 chi[] (alpha[] (exKh[] + 2 Theta[]) - dbeta[-i, i])
  ];

(* (2) *)

SetEQN[dtgamt[i_, j_], -2 alpha[] exAt[i, j] + (gamt[-k, i] dbeta[j, k]
   + gamt[-k, j] dbeta[i, k]) - 2/3 gamt[i, j] dbeta[-k, k]];

(* (3) *)

SetEQN[dtexKh[], -invgam[k, l] DDalpha[-k, -l] + alpha[] (exAt[-k, -l] 
  exAtUU[k, l] + 1/3 (exKh[] + 2 Theta[]) ^ 2) + 4 cpi alpha[] (trSs[] + 
  rho[]) + alpha[] ckappa1 (1 - ckappa2) Theta[]];

(* (4) *)

SetEQN[dtexAt[i_, j_], chi[] ((-DDalpha[i, j] + alpha[] (R[i, j] - 8 cpi
   Ss[i, j])) - 1/3 gam[i, j] invgam[k, l] (-DDalpha[-k, -l] + alpha[] (R[
  -k, -l] - 8 cpi Ss[-k, -l]))) + alpha[] ((exKh[] + 2 Theta[]) exAt[i, j
  ] - 2 invgamt[k, l] exAt[-k, i] exAt[-l, j]) + (exAt[-k, i] dbeta[j, k]
   + exAt[-k, j] dbeta[i, k]) - 2/3 exAt[i, j] dbeta[-k, k]];

(* (5) *)

SetEQN[dttrGt[i_], -2 exAtUU[i, j] dalpha[-j] + 2 alpha[] (Gt[i, -j, -k
  ] exAtUU[j, k] - 3/2 exAtUU[i, j] dchi[-j] / chi[] - 1/3 invgamt[i, j] 
  (2 dexKh[-j] + dTheta[-j]) - 8 cpi invgam[i, j] Sm[-j]) + invgamt[j, k]
   ddbeta[-j, -k, i] + 1/3 invgamt[i, j] ddbeta[-j, -k, k] - 2/3 trGtd[i]
   dbeta[-j, j] - 2 alpha[] ckappa1 (trGt[i] - trGtd[i])];

(* (6) *)

SetEQN[dtTheta[], 1/2 alpha[] (trR[] - exAt[-k, -l] exAtUU[k, l] + 2/3 
  (exKh[] + 2 Theta[]) ^ 2) - alpha[] (8 cpi rho[] + ckappa1 (2 + ckappa2
  ) Theta[])];

(* (11) *)

SetEQN[dtalpha[], -alpha[] ^ 2 cmuL exKh[]];

(* (12) *)

SetEQN[dtbeta[i_], alpha[] ^ 2 cmuS trGt[i] - ceta beta[i]];
