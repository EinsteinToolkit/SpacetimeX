(* ::Package:: *)

(* Z4c_rhs.wl *)

(* (c) Liwei Ji, 07/2024 *)

(* See arXiv:1212.2901 [gr-qc] *)

(****************)

(* Intermediate *)

(****************)

SetEQN[dlnW[i_], W[] ^ -1 dW[i]];

SetEQN[invgam[i_, j_], W[] ^ 2 invgamt[i, j]];

SetEQN[gam[i_, j_], W[] ^ -2 gamt[i, j]];

SetEQN[GtDDD[k_, i_, j_], 1/2 (dgamt[i, j, k] + dgamt[j, k, i] - dgamt[k, i, j])];

SetEQN[GtDDU[i_, j_, k_], invgamt[k, l] GtDDD[i, j, -l]];

SetEQN[Gt[k_, i_, j_], invgamt[k, l] GtDDD[-l, i, j]];

SetEQN[trGtd[i_], invgamt[k, l] Gt[i, -k, -l]];

SetEQN[dgam[k_, i_, j_], W[] ^ -2 (-2 dlnW[k] gamt[i, j] + dgamt[k, i, j])];

SetEQN[GamDDD[k_, i_, j_], 1/2 (dgam[i, j, k] + dgam[j, k, i] - dgam[k, i, j])];

SetEQN[Gam[k_, i_, j_], invgam[k, l] GamDDD[-l, i, j]];

SetEQN[exAtUU[i_, j_], invgamt[i, k] invgamt[j, l] exAt[-k, -l]];

SetEQN[tDtDW[i_, j_], ddW[i, j] - Gt[k, i, j] dW[-k]];

SetEQN[DDalpha[i_, j_], ddalpha[i, j] - Gam[k, i, j] dalpha[-k]];

(* (8) *)

SetEQN[RtW[i_, j_], 1 / W[] tDtDW[i, j] + 1 / W[] gamt[i, j] invgamt[k, l] tDtDW[-k, -l] - 2 gamt[i, j] invgamt[k, l] dlnW[-k] dlnW[-l]];

(* (9) *)

SetEQN[Rt[i_, j_], -(1/2) invgamt[l, m] ddgamt[-l, -m, i, j] + 1/2 (gamt[-k, i] dtrGt[j, k] + gamt[-k, j] dtrGt[i, k]) + 1/2 trGtd[k] (GtDDD[i, j, -k] + GtDDD[j, i, -k]) + ((Gt[k, -l, i] GtDDU[j, -k, l] + Gt[k, -l, j] GtDDU[i, -k, l]) + Gt[k, i, -m] GtDDU[-k, j, m])];

(* (10) *)

SetEQN[R[i_, j_], RtW[i, j] + Rt[i, j]];

SetEQN[trR[], invgam[k, l] R[-k, -l]];

SetEQN[rho[], alpha[] ^ -2 (eTtt[] - 2 beta[j] eTt[-j] + beta[i] beta[j] eT[-i, -j])];

SetEQN[Sm[i_], -alpha[] ^ -1 (eTt[i] - beta[k] eT[-k, i])];

SetEQN[Ss[i_, j_], eT[i, j]];

SetEQN[trSs[], invgam[k, l] Ss[-k, -l]];

SetEQN[trdexAtUU[i_], -invgamt[i, l] exAtUU[j, m] dgamt[-j, -l, -m] - invgamt[j, l] exAtUU[i, m] dgamt[-j, -l, -m] + invgamt[i, l] invgamt[j, m] dexAt[-j, -l, -m]];

(* (13) *)

SetEQN[ZtC[i_], (trGt[i] - trGtd[i]) / 2];

(* (14) *)

SetEQN[HC[], trR[] - exAt[-k, -l] exAtUU[k, l] + 2/3 (exKh[] + 2 Theta[]) ^ 2 - 16 cpi rho[]];

(* (15) *)

SetEQN[MtC[i_], trdexAtUU[i] + Gt[i, -j, -k] exAtUU[j, k] - 2/3 invgamt[i, j] (dexKh[-j] + 2 dTheta[-j]) - 3 exAtUU[i, j] dlnW[-j] - 8 cpi invgamt[i, j] Sm[-j]];

(*******)

(* EOM *)

(*******)

(* (1) *)

SetEQN[dtW[], beta[k] dW[-k] + 1/3 W[] (alpha[] (exKh[] + 2 Theta[]) - dbeta[-i, i])];

(* (2) *)

SetEQN[dtgamt[i_, j_], beta[k] dgamt[-k, i, j] - 2 alpha[] exAt[i, j] + (gamt[-k, i] dbeta[j, k] + gamt[-k, j] dbeta[i, k]) - 2/3 gamt[i, j] dbeta[-k, k]];

(* (3) *)

SetEQN[dtexKh[], beta[k] dexKh[-k] - invgam[k, l] DDalpha[-k, -l] + alpha[] (exAt[-k, -l] exAtUU[k, l] + 1/3 (exKh[] + 2 Theta[]) ^ 2) + 4 cpi alpha[] (trSs[] + rho[]) + alpha[] ckappa1 (1 - ckappa2) Theta[]];

(* (4) *)

SetEQN[dtexAt[i_, j_], beta[k] dexAt[-k, i, j] + W[] ^ 2 ((-DDalpha[i, j] + alpha[] (R[i, j] - 8 cpi Ss[i, j])) - 1/3 gam[i, j] invgam[k, l] (-DDalpha[-k, -l] + alpha[] (R[-k, -l] - 8 cpi Ss[-k, -l]))) + alpha[] ((exKh[] + 2 Theta[]) exAt[i, j] - 2 invgamt[k, l] exAt[-k, i] exAt[-l, j]) + (exAt[-k, i] dbeta[j, k] + exAt[-k, j] dbeta[i, k]) - 2/3 exAt[i, j] dbeta[-k, k]];

(* (5) *)

SetEQN[dttrGt[i_], beta[k] dtrGt[-k, i] - 2 exAtUU[i, j] dalpha[-j] + 2 alpha[] (Gt[i, -j, -k] exAtUU[j, k] - 3 exAtUU[i, j] dlnW[-j] - 1/3 invgamt[i, j] (2 dexKh[-j] + dTheta[-j]) - 8 cpi invgamt[i, j] Sm[-j]) + invgamt[j, k] ddbeta[-j, -k, i] + 1/3 invgamt[i, j] ddbeta[-j, -k, k] - trGtd[j] dbeta[-j, i] + 2/3 trGtd[i] dbeta[-j, j] - 2 alpha[] ckappa1 (trGt[i] - trGtd[i])];

(* (6) *)

SetEQN[dtTheta[], beta[k] dTheta[-k] + 1/2 alpha[] (trR[] - exAt[-k, -l] exAtUU[k, l] + 2/3 (exKh[] + 2 Theta[]) ^ 2) - alpha[] (8 cpi rho[] + ckappa1 (2 + ckappa2) Theta[])];

(* (11) *)

SetEQN[dtalpha[], beta[k] dalpha[-k] - alpha[] cmuL exKh[]];

(* (12) *)

SetEQN[dtbeta[i_], beta[k] dbeta[-k, i] + cmuS trGt[i] - ceta beta[i]];
