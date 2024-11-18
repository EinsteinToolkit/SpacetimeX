(* ::Package:: *)

(* ADM_rhs.wl *)

(* (c) Liwei Ji, 11/2024 *)

(****************)

(* Intermediate *)

(****************)

SetEQN[GamDDD[k_, i_, j_], 1/2 (dgam[i, j, k] + dgam[j, k, i] - dgam[k, i, j])];

SetEQN[Gam[k_, i_, j_], invgam[k, l] GamDDD[-l, i, j]];

SetEQN[tr1dGam[i_, j_], -invgam[k, p] invgam[l, q] dgam[-k, -p, -q] GamDDD[-l, i, j] + 1/2 invgam[k, l] (ddgam[-k, i, j, -l] + ddgam[-k, j, -l, i] - ddgam[-k, -l, i, j])];

SetEQN[tr2dGam[m_, i_], 1/2 (invgam[k, l] ddgam[m, i, -k, -l] - invgam[k, p] invgam[l, q] dgam[m, -p, -q] dgam[i, -k, -l])];

SetEQN[R[j_, k_], tr1dGam[j, k] - tr2dGam[j, k] + Gam[i, -i, -p] Gam[p, j, k] - Gam[i, j, -p] Gam[p, -i, k]];

SetEQN[trexK[], invgam[k, l] exK[-k, -l]];

SetEQN[DexK[k_, i_, j_], dexK[k, i, j] - Gam[l, k, i] exK[-l, j] - Gam[l, k, j] exK[-l, i]];

(* matter *)

SetEQN[rho[], alpha[]^-2 (eTtt[] - 2 beta[j] eTt[-j] + beta[i] beta[j] eT[-i, -j])];

SetEQN[Sm[i_], -alpha[]^-1 (eTt[i] - beta[k] eT[-k, i])];

(***************)

(* Constraints *)

(***************)

SetEQN[HC[], invgam[k, l] R[-k, -l] + trexK[]^2 - invgam[i, k] invgam[j, l] exK[-i, -j] exK[-k, -l] - 16 cpi rho[]];

SetEQN[MC[i_], (invgam[i, k] invgam[j, l] - invgam[i, j] invgam[k, l]) DexK[-j, -k, -l] - 8 cpi invgam[i, j] Sm[-j]];
