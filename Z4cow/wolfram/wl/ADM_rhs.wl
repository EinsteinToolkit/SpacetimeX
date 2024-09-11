(* ::Package:: *)

(* ADM_rhs.wl *)

(* (c) Liwei Ji, 07/2024 *)

(* See arXiv:1212.2901 [gr-qc] *)

SetEQN[ADMgam[i_, j_], W[] ^ -2 gamt[i, j]];

SetEQN[ADMK[i_, j_], W[] ^ -2 (exAt[i, j] + 1/3 (exKh[] + 2 Theta[]) gamt[i, j])];

SetEQN[ADMalpha[], alpha[]];

SetEQN[ADMbeta[i_], beta[i]];

SetEQN[ADMdtalpha[], dtalpha[]];

SetEQN[ADMdtbeta[i_], dtbeta[i]];
