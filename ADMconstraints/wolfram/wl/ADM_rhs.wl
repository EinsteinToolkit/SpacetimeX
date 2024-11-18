(* ::Package:: *)

(* ADM_rhs.wl *)

(* (c) Liwei Ji, 11/2024 *)

(****************)

(* Intermediate *)

(****************)

SetEQN[GamDDD[k_, i_, j_], 1/2 (dADMgam[i, j, k] + dADMgam[j, k, i] - dADMgam[k, i, j])];

SetEQN[Gam[k_, i_, j_], invgam[k, l] GamDDD[-l, i, j]];

SetEQN[tr1dGam[i_, j_], -invgam[k, p] invgam[l, q] dADMgam[-k, -p, -q] GamDDD[-l, i, j] + 1/2 invgam[k, l] (ddADMgam[-k, i, j, -l] + ddADMgam[-k, j, -l, i] - ddADMgam[-k, -l, i, j])];

SetEQN[tr2dGam[m_, i_], 1/2 (invgam[k, l] ddADMgam[m, i, -k, -l] - invgam[k, p] invgam[l, q] dADMgam[m, -p, -q] dADMgam[i, -k, -l])];

SetEQN[R[j_, k_], tr1dGam[j, k] - tr2dGam[j, k] + Gam[i, -i, -p] Gam[p, j, k] - Gam[i, j, -p] Gam[p, -i, k]];

SetEQN[trK[], invgam[k, l] ADMK[-k, -l]];

SetEQN[DADMK[k_, i_, j_], dADMK[k, i, j] - Gam[l, k, i] ADMK[-l, j] - Gam[l, k, j] ADMK[-l, i]];

(* matter *)

SetEQN[rho[], ADMalpha[]^-2 (eTtt[] - 2 ADMbeta[j] eTt[-j] + ADMbeta[i] ADMbeta[j] eT[-i, -j])];

SetEQN[Sm[i_], -ADMalpha[]^-1 (eTt[i] - ADMbeta[k] eT[-k, i])];

(***************)

(* Constraints *)

(***************)

SetEQN[HC[], invgam[k, l] R[-k, -l] + trK[]^2 - invgam[i, k] invgam[j, l] ADMK[-i, -j] ADMK[-k, -l] - 16 cpi rho[]];

SetEQN[MtC[i_], (invgam[i, k] invgam[j, l] - invgam[i, j] invgam[k, l]) DADMK[-j, -k, -l] - 8 cpi invgam[i, j] Sm[-j]];
