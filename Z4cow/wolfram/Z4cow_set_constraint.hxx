/* Z4cow_set_constraint.hxx */
/* Produced with Mathematica */

#ifndef Z4COW_SET_CONSTRAINT_HXX
#define Z4COW_SET_CONSTRAINT_HXX

const GF3D2<CCTK_REAL> &local_ZtC1 = gf_ZtC(0);
const GF3D2<CCTK_REAL> &local_ZtC2 = gf_ZtC(1);
const GF3D2<CCTK_REAL> &local_ZtC3 = gf_ZtC(2);
const GF3D2<CCTK_REAL> &local_HC = gf_HC;
const GF3D2<CCTK_REAL> &local_MtC1 = gf_MtC(0);
const GF3D2<CCTK_REAL> &local_MtC2 = gf_MtC(1);
const GF3D2<CCTK_REAL> &local_MtC3 = gf_MtC(2);

noinline([&]() __attribute__((__flatten__, __hot__)) {
  grid.loop_int_device<0, 0, 0, vsize>(
    grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
    const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
    const GF3D2index index2(layout2, p.I);
    const GF3D5index index5(layout5, p.I);

const auto &tmp_eTtt = gf_eTtt(mask, index2);
const auto &tmp_eTt = gf_eTt(mask, index2);
const auto &tmp_eT = gf_eT(mask, index2);
const auto &tmp_W = tl_W(mask, index5);
const auto &tmp_gamt = tl_gamt(mask, index5);
const auto &tmp_exKh = tl_exKh(mask, index5);
const auto &tmp_exAt = tl_exAt(mask, index5);
const auto &tmp_trGt = tl_trGt(mask, index5);
const auto &tmp_Theta = tl_Theta(mask, index5);
const auto &tmp_alpha = tl_alpha(mask, index5);
const auto &tmp_beta = tl_beta(mask, index5);
const auto &tmp_dW = tl_dW(mask, index5);
const auto &tmp_dgamt = tl_dgamt(mask, index5);
const auto &tmp_dexKh = tl_dexKh(mask, index5);
const auto &tmp_dexAt = tl_dexAt(mask, index5);
const auto &tmp_dtrGt = tl_dtrGt(mask, index5);
const auto &tmp_dTheta = tl_dTheta(mask, index5);
const auto &tmp_ddW = tl_ddW(mask, index5);
const auto &tmp_ddgamt = tl_ddgamt(mask, index5);

const vreal eTtt = tmp_eTtt;
const vreal eTt1 = tmp_eTt(0);
const vreal eTt2 = tmp_eTt(1);
const vreal eTt3 = tmp_eTt(2);
const vreal eT11 = tmp_eT(0,0);
const vreal eT12 = tmp_eT(0,1);
const vreal eT13 = tmp_eT(0,2);
const vreal eT22 = tmp_eT(1,1);
const vreal eT23 = tmp_eT(1,2);
const vreal eT33 = tmp_eT(2,2);
const vreal W = tmp_W;
const vreal gamt11 = tmp_gamt(0,0);
const vreal gamt12 = tmp_gamt(0,1);
const vreal gamt13 = tmp_gamt(0,2);
const vreal gamt22 = tmp_gamt(1,1);
const vreal gamt23 = tmp_gamt(1,2);
const vreal gamt33 = tmp_gamt(2,2);
const vreal exKh = tmp_exKh;
const vreal exAt11 = tmp_exAt(0,0);
const vreal exAt12 = tmp_exAt(0,1);
const vreal exAt13 = tmp_exAt(0,2);
const vreal exAt22 = tmp_exAt(1,1);
const vreal exAt23 = tmp_exAt(1,2);
const vreal exAt33 = tmp_exAt(2,2);
const vreal trGt1 = tmp_trGt(0);
const vreal trGt2 = tmp_trGt(1);
const vreal trGt3 = tmp_trGt(2);
const vreal Theta = tmp_Theta;
const vreal alpha = tmp_alpha;
const vreal beta1 = tmp_beta(0);
const vreal beta2 = tmp_beta(1);
const vreal beta3 = tmp_beta(2);
const vreal dW1 = tmp_dW(0);
const vreal dW2 = tmp_dW(1);
const vreal dW3 = tmp_dW(2);
const vreal dgamt111 = tmp_dgamt(0,0)(0);
const vreal dgamt112 = tmp_dgamt(0,1)(0);
const vreal dgamt113 = tmp_dgamt(0,2)(0);
const vreal dgamt122 = tmp_dgamt(1,1)(0);
const vreal dgamt123 = tmp_dgamt(1,2)(0);
const vreal dgamt133 = tmp_dgamt(2,2)(0);
const vreal dgamt211 = tmp_dgamt(0,0)(1);
const vreal dgamt212 = tmp_dgamt(0,1)(1);
const vreal dgamt213 = tmp_dgamt(0,2)(1);
const vreal dgamt222 = tmp_dgamt(1,1)(1);
const vreal dgamt223 = tmp_dgamt(1,2)(1);
const vreal dgamt233 = tmp_dgamt(2,2)(1);
const vreal dgamt311 = tmp_dgamt(0,0)(2);
const vreal dgamt312 = tmp_dgamt(0,1)(2);
const vreal dgamt313 = tmp_dgamt(0,2)(2);
const vreal dgamt322 = tmp_dgamt(1,1)(2);
const vreal dgamt323 = tmp_dgamt(1,2)(2);
const vreal dgamt333 = tmp_dgamt(2,2)(2);
const vreal dexKh1 = tmp_dexKh(0);
const vreal dexKh2 = tmp_dexKh(1);
const vreal dexKh3 = tmp_dexKh(2);
const vreal dexAt111 = tmp_dexAt(0,0)(0);
const vreal dexAt112 = tmp_dexAt(0,1)(0);
const vreal dexAt113 = tmp_dexAt(0,2)(0);
const vreal dexAt122 = tmp_dexAt(1,1)(0);
const vreal dexAt123 = tmp_dexAt(1,2)(0);
const vreal dexAt133 = tmp_dexAt(2,2)(0);
const vreal dexAt211 = tmp_dexAt(0,0)(1);
const vreal dexAt212 = tmp_dexAt(0,1)(1);
const vreal dexAt213 = tmp_dexAt(0,2)(1);
const vreal dexAt222 = tmp_dexAt(1,1)(1);
const vreal dexAt223 = tmp_dexAt(1,2)(1);
const vreal dexAt233 = tmp_dexAt(2,2)(1);
const vreal dexAt311 = tmp_dexAt(0,0)(2);
const vreal dexAt312 = tmp_dexAt(0,1)(2);
const vreal dexAt313 = tmp_dexAt(0,2)(2);
const vreal dexAt322 = tmp_dexAt(1,1)(2);
const vreal dexAt323 = tmp_dexAt(1,2)(2);
const vreal dexAt333 = tmp_dexAt(2,2)(2);
const vreal dtrGt11 = tmp_dtrGt(0)(0);
const vreal dtrGt12 = tmp_dtrGt(1)(0);
const vreal dtrGt13 = tmp_dtrGt(2)(0);
const vreal dtrGt21 = tmp_dtrGt(0)(1);
const vreal dtrGt22 = tmp_dtrGt(1)(1);
const vreal dtrGt23 = tmp_dtrGt(2)(1);
const vreal dtrGt31 = tmp_dtrGt(0)(2);
const vreal dtrGt32 = tmp_dtrGt(1)(2);
const vreal dtrGt33 = tmp_dtrGt(2)(2);
const vreal dTheta1 = tmp_dTheta(0);
const vreal dTheta2 = tmp_dTheta(1);
const vreal dTheta3 = tmp_dTheta(2);
const vreal ddW11 = tmp_ddW(0,0);
const vreal ddW12 = tmp_ddW(0,1);
const vreal ddW13 = tmp_ddW(0,2);
const vreal ddW22 = tmp_ddW(1,1);
const vreal ddW23 = tmp_ddW(1,2);
const vreal ddW33 = tmp_ddW(2,2);
const vreal ddgamt1111 = tmp_ddgamt(0,0)(0,0);
const vreal ddgamt1112 = tmp_ddgamt(0,1)(0,0);
const vreal ddgamt1113 = tmp_ddgamt(0,2)(0,0);
const vreal ddgamt1122 = tmp_ddgamt(1,1)(0,0);
const vreal ddgamt1123 = tmp_ddgamt(1,2)(0,0);
const vreal ddgamt1133 = tmp_ddgamt(2,2)(0,0);
const vreal ddgamt1211 = tmp_ddgamt(0,0)(0,1);
const vreal ddgamt1212 = tmp_ddgamt(0,1)(0,1);
const vreal ddgamt1213 = tmp_ddgamt(0,2)(0,1);
const vreal ddgamt1222 = tmp_ddgamt(1,1)(0,1);
const vreal ddgamt1223 = tmp_ddgamt(1,2)(0,1);
const vreal ddgamt1233 = tmp_ddgamt(2,2)(0,1);
const vreal ddgamt1311 = tmp_ddgamt(0,0)(0,2);
const vreal ddgamt1312 = tmp_ddgamt(0,1)(0,2);
const vreal ddgamt1313 = tmp_ddgamt(0,2)(0,2);
const vreal ddgamt1322 = tmp_ddgamt(1,1)(0,2);
const vreal ddgamt1323 = tmp_ddgamt(1,2)(0,2);
const vreal ddgamt1333 = tmp_ddgamt(2,2)(0,2);
const vreal ddgamt2211 = tmp_ddgamt(0,0)(1,1);
const vreal ddgamt2212 = tmp_ddgamt(0,1)(1,1);
const vreal ddgamt2213 = tmp_ddgamt(0,2)(1,1);
const vreal ddgamt2222 = tmp_ddgamt(1,1)(1,1);
const vreal ddgamt2223 = tmp_ddgamt(1,2)(1,1);
const vreal ddgamt2233 = tmp_ddgamt(2,2)(1,1);
const vreal ddgamt2311 = tmp_ddgamt(0,0)(1,2);
const vreal ddgamt2312 = tmp_ddgamt(0,1)(1,2);
const vreal ddgamt2313 = tmp_ddgamt(0,2)(1,2);
const vreal ddgamt2322 = tmp_ddgamt(1,1)(1,2);
const vreal ddgamt2323 = tmp_ddgamt(1,2)(1,2);
const vreal ddgamt2333 = tmp_ddgamt(2,2)(1,2);
const vreal ddgamt3311 = tmp_ddgamt(0,0)(2,2);
const vreal ddgamt3312 = tmp_ddgamt(0,1)(2,2);
const vreal ddgamt3313 = tmp_ddgamt(0,2)(2,2);
const vreal ddgamt3322 = tmp_ddgamt(1,1)(2,2);
const vreal ddgamt3323 = tmp_ddgamt(1,2)(2,2);
const vreal ddgamt3333 = tmp_ddgamt(2,2)(2,2);

vreal dlnW1
=
dW1/W
;

vreal dlnW2
=
dW2/W
;

vreal dlnW3
=
dW3/W
;

vreal invgamt11
=
-Power(gamt23,2) + gamt22*gamt33
;

vreal invgamt12
=
gamt13*gamt23 - gamt12*gamt33
;

vreal invgamt13
=
-(gamt13*gamt22) + gamt12*gamt23
;

vreal invgamt22
=
-Power(gamt13,2) + gamt11*gamt33
;

vreal invgamt23
=
gamt12*gamt13 - gamt11*gamt23
;

vreal invgamt33
=
-Power(gamt12,2) + gamt11*gamt22
;

vreal gam11
=
gamt11/Power(W,2)
;

vreal gam12
=
gamt12/Power(W,2)
;

vreal gam13
=
gamt13/Power(W,2)
;

vreal gam22
=
gamt22/Power(W,2)
;

vreal gam23
=
gamt23/Power(W,2)
;

vreal gam33
=
gamt33/Power(W,2)
;

vreal GtDDD111
=
dgamt111/2.
;

vreal GtDDD112
=
dgamt211/2.
;

vreal GtDDD113
=
dgamt311/2.
;

vreal GtDDD122
=
-0.5*dgamt122 + dgamt212
;

vreal GtDDD123
=
(-dgamt123 + dgamt213 + dgamt312)/2.
;

vreal GtDDD133
=
-0.5*dgamt133 + dgamt313
;

vreal GtDDD211
=
dgamt112 - dgamt211/2.
;

vreal GtDDD212
=
dgamt122/2.
;

vreal GtDDD213
=
(dgamt123 - dgamt213 + dgamt312)/2.
;

vreal GtDDD222
=
dgamt222/2.
;

vreal GtDDD223
=
dgamt322/2.
;

vreal GtDDD233
=
-0.5*dgamt233 + dgamt323
;

vreal GtDDD311
=
dgamt113 - dgamt311/2.
;

vreal GtDDD312
=
(dgamt123 + dgamt213 - dgamt312)/2.
;

vreal GtDDD313
=
dgamt133/2.
;

vreal GtDDD322
=
dgamt223 - dgamt322/2.
;

vreal GtDDD323
=
dgamt233/2.
;

vreal GtDDD333
=
dgamt333/2.
;

vreal GtDDU111
=
GtDDD111*invgamt11 + GtDDD112*invgamt12 + GtDDD113*invgamt13
;

vreal GtDDU112
=
GtDDD111*invgamt12 + GtDDD112*invgamt22 + GtDDD113*invgamt23
;

vreal GtDDU113
=
GtDDD111*invgamt13 + GtDDD112*invgamt23 + GtDDD113*invgamt33
;

vreal GtDDU121
=
GtDDD112*invgamt11 + GtDDD122*invgamt12 + GtDDD123*invgamt13
;

vreal GtDDU122
=
GtDDD112*invgamt12 + GtDDD122*invgamt22 + GtDDD123*invgamt23
;

vreal GtDDU123
=
GtDDD112*invgamt13 + GtDDD122*invgamt23 + GtDDD123*invgamt33
;

vreal GtDDU131
=
GtDDD113*invgamt11 + GtDDD123*invgamt12 + GtDDD133*invgamt13
;

vreal GtDDU132
=
GtDDD113*invgamt12 + GtDDD123*invgamt22 + GtDDD133*invgamt23
;

vreal GtDDU133
=
GtDDD113*invgamt13 + GtDDD123*invgamt23 + GtDDD133*invgamt33
;

vreal GtDDU211
=
GtDDD211*invgamt11 + GtDDD212*invgamt12 + GtDDD213*invgamt13
;

vreal GtDDU212
=
GtDDD211*invgamt12 + GtDDD212*invgamt22 + GtDDD213*invgamt23
;

vreal GtDDU213
=
GtDDD211*invgamt13 + GtDDD212*invgamt23 + GtDDD213*invgamt33
;

vreal GtDDU221
=
GtDDD212*invgamt11 + GtDDD222*invgamt12 + GtDDD223*invgamt13
;

vreal GtDDU222
=
GtDDD212*invgamt12 + GtDDD222*invgamt22 + GtDDD223*invgamt23
;

vreal GtDDU223
=
GtDDD212*invgamt13 + GtDDD222*invgamt23 + GtDDD223*invgamt33
;

vreal GtDDU231
=
GtDDD213*invgamt11 + GtDDD223*invgamt12 + GtDDD233*invgamt13
;

vreal GtDDU232
=
GtDDD213*invgamt12 + GtDDD223*invgamt22 + GtDDD233*invgamt23
;

vreal GtDDU233
=
GtDDD213*invgamt13 + GtDDD223*invgamt23 + GtDDD233*invgamt33
;

vreal GtDDU311
=
GtDDD311*invgamt11 + GtDDD312*invgamt12 + GtDDD313*invgamt13
;

vreal GtDDU312
=
GtDDD311*invgamt12 + GtDDD312*invgamt22 + GtDDD313*invgamt23
;

vreal GtDDU313
=
GtDDD311*invgamt13 + GtDDD312*invgamt23 + GtDDD313*invgamt33
;

vreal GtDDU321
=
GtDDD312*invgamt11 + GtDDD322*invgamt12 + GtDDD323*invgamt13
;

vreal GtDDU322
=
GtDDD312*invgamt12 + GtDDD322*invgamt22 + GtDDD323*invgamt23
;

vreal GtDDU323
=
GtDDD312*invgamt13 + GtDDD322*invgamt23 + GtDDD323*invgamt33
;

vreal GtDDU331
=
GtDDD313*invgamt11 + GtDDD323*invgamt12 + GtDDD333*invgamt13
;

vreal GtDDU332
=
GtDDD313*invgamt12 + GtDDD323*invgamt22 + GtDDD333*invgamt23
;

vreal GtDDU333
=
GtDDD313*invgamt13 + GtDDD323*invgamt23 + GtDDD333*invgamt33
;

vreal Gt111
=
GtDDD111*invgamt11 + GtDDD211*invgamt12 + GtDDD311*invgamt13
;

vreal Gt112
=
GtDDD112*invgamt11 + GtDDD212*invgamt12 + GtDDD312*invgamt13
;

vreal Gt113
=
GtDDD113*invgamt11 + GtDDD213*invgamt12 + GtDDD313*invgamt13
;

vreal Gt122
=
GtDDD122*invgamt11 + GtDDD222*invgamt12 + GtDDD322*invgamt13
;

vreal Gt123
=
GtDDD123*invgamt11 + GtDDD223*invgamt12 + GtDDD323*invgamt13
;

vreal Gt133
=
GtDDD133*invgamt11 + GtDDD233*invgamt12 + GtDDD333*invgamt13
;

vreal Gt211
=
GtDDD111*invgamt12 + GtDDD211*invgamt22 + GtDDD311*invgamt23
;

vreal Gt212
=
GtDDD112*invgamt12 + GtDDD212*invgamt22 + GtDDD312*invgamt23
;

vreal Gt213
=
GtDDD113*invgamt12 + GtDDD213*invgamt22 + GtDDD313*invgamt23
;

vreal Gt222
=
GtDDD122*invgamt12 + GtDDD222*invgamt22 + GtDDD322*invgamt23
;

vreal Gt223
=
GtDDD123*invgamt12 + GtDDD223*invgamt22 + GtDDD323*invgamt23
;

vreal Gt233
=
GtDDD133*invgamt12 + GtDDD233*invgamt22 + GtDDD333*invgamt23
;

vreal Gt311
=
GtDDD111*invgamt13 + GtDDD211*invgamt23 + GtDDD311*invgamt33
;

vreal Gt312
=
GtDDD112*invgamt13 + GtDDD212*invgamt23 + GtDDD312*invgamt33
;

vreal Gt313
=
GtDDD113*invgamt13 + GtDDD213*invgamt23 + GtDDD313*invgamt33
;

vreal Gt322
=
GtDDD122*invgamt13 + GtDDD222*invgamt23 + GtDDD322*invgamt33
;

vreal Gt323
=
GtDDD123*invgamt13 + GtDDD223*invgamt23 + GtDDD323*invgamt33
;

vreal Gt333
=
GtDDD133*invgamt13 + GtDDD233*invgamt23 + GtDDD333*invgamt33
;

vreal trGtd1
=
Gt111*invgamt11 + 2*Gt112*invgamt12 + 2*Gt113*invgamt13 + Gt122*invgamt22 + 
  2*Gt123*invgamt23 + Gt133*invgamt33
;

vreal trGtd2
=
Gt211*invgamt11 + 2*Gt212*invgamt12 + 2*Gt213*invgamt13 + Gt222*invgamt22 + 
  2*Gt223*invgamt23 + Gt233*invgamt33
;

vreal trGtd3
=
Gt311*invgamt11 + 2*Gt312*invgamt12 + 2*Gt313*invgamt13 + Gt322*invgamt22 + 
  2*Gt323*invgamt23 + Gt333*invgamt33
;

vreal exAtUU11
=
exAt11*Power(invgamt11,2) + 2*exAt12*invgamt11*invgamt12 + 
  exAt22*Power(invgamt12,2) + 2*exAt13*invgamt11*invgamt13 + 
  2*exAt23*invgamt12*invgamt13 + exAt33*Power(invgamt13,2)
;

vreal exAtUU12
=
exAt11*invgamt11*invgamt12 + exAt13*invgamt12*invgamt13 + 
  exAt22*invgamt12*invgamt22 + exAt23*invgamt13*invgamt22 + 
  exAt12*(Power(invgamt12,2) + invgamt11*invgamt22) + 
  exAt13*invgamt11*invgamt23 + exAt23*invgamt12*invgamt23 + 
  exAt33*invgamt13*invgamt23
;

vreal exAtUU13
=
exAt11*invgamt11*invgamt13 + exAt12*invgamt12*invgamt13 + 
  exAt13*Power(invgamt13,2) + exAt12*invgamt11*invgamt23 + 
  exAt22*invgamt12*invgamt23 + exAt23*invgamt13*invgamt23 + 
  exAt13*invgamt11*invgamt33 + exAt23*invgamt12*invgamt33 + 
  exAt33*invgamt13*invgamt33
;

vreal exAtUU22
=
exAt11*Power(invgamt12,2) + 2*exAt12*invgamt12*invgamt22 + 
  exAt22*Power(invgamt22,2) + 2*exAt13*invgamt12*invgamt23 + 
  2*exAt23*invgamt22*invgamt23 + exAt33*Power(invgamt23,2)
;

vreal exAtUU23
=
exAt11*invgamt12*invgamt13 + exAt12*invgamt13*invgamt22 + 
  exAt12*invgamt12*invgamt23 + exAt13*invgamt13*invgamt23 + 
  exAt22*invgamt22*invgamt23 + exAt23*Power(invgamt23,2) + 
  exAt13*invgamt12*invgamt33 + exAt23*invgamt22*invgamt33 + 
  exAt33*invgamt23*invgamt33
;

vreal exAtUU33
=
exAt11*Power(invgamt13,2) + 2*exAt12*invgamt13*invgamt23 + 
  exAt22*Power(invgamt23,2) + 2*exAt13*invgamt13*invgamt33 + 
  2*exAt23*invgamt23*invgamt33 + exAt33*Power(invgamt33,2)
;

vreal tDtDW11
=
ddW11 - dW1*Gt111 - dW2*Gt211 - dW3*Gt311
;

vreal tDtDW12
=
ddW12 - dW1*Gt112 - dW2*Gt212 - dW3*Gt312
;

vreal tDtDW13
=
ddW13 - dW1*Gt113 - dW2*Gt213 - dW3*Gt313
;

vreal tDtDW22
=
ddW22 - dW1*Gt122 - dW2*Gt222 - dW3*Gt322
;

vreal tDtDW23
=
ddW23 - dW1*Gt123 - dW2*Gt223 - dW3*Gt323
;

vreal tDtDW33
=
ddW33 - dW1*Gt133 - dW2*Gt233 - dW3*Gt333
;

vreal RtW11
=
(tDtDW11 + gamt11*invgamt11*tDtDW11 + 
    gamt11*(2*invgamt13*tDtDW13 + invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 + 
       invgamt33*tDtDW33 - 2*Power(dlnW1,2)*invgamt11*W - 
       4*dlnW1*dlnW3*invgamt13*W - 2*Power(dlnW2,2)*invgamt22*W - 
       4*dlnW2*dlnW3*invgamt23*W - 2*Power(dlnW3,2)*invgamt33*W + 
       2*invgamt12*(tDtDW12 - 2*dlnW1*dlnW2*W)))/W
;

vreal RtW12
=
(tDtDW12 + gamt12*(2*invgamt12*tDtDW12 + 2*invgamt13*tDtDW13 + 
       invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 + invgamt33*tDtDW33 - 
       4*dlnW1*dlnW2*invgamt12*W - 4*dlnW1*dlnW3*invgamt13*W - 
       2*Power(dlnW2,2)*invgamt22*W - 4*dlnW2*dlnW3*invgamt23*W - 
       2*Power(dlnW3,2)*invgamt33*W + 
       invgamt11*(tDtDW11 - 2*Power(dlnW1,2)*W)))/W
;

vreal RtW13
=
(tDtDW13 + gamt13*(2*invgamt12*tDtDW12 + 2*invgamt13*tDtDW13 + 
       invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 + invgamt33*tDtDW33 - 
       4*dlnW1*dlnW2*invgamt12*W - 4*dlnW1*dlnW3*invgamt13*W - 
       2*Power(dlnW2,2)*invgamt22*W - 4*dlnW2*dlnW3*invgamt23*W - 
       2*Power(dlnW3,2)*invgamt33*W + 
       invgamt11*(tDtDW11 - 2*Power(dlnW1,2)*W)))/W
;

vreal RtW22
=
(tDtDW22 + gamt22*(2*invgamt12*tDtDW12 + 2*invgamt13*tDtDW13 + 
       invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 + invgamt33*tDtDW33 - 
       4*dlnW1*dlnW2*invgamt12*W - 4*dlnW1*dlnW3*invgamt13*W - 
       2*Power(dlnW2,2)*invgamt22*W - 4*dlnW2*dlnW3*invgamt23*W - 
       2*Power(dlnW3,2)*invgamt33*W + 
       invgamt11*(tDtDW11 - 2*Power(dlnW1,2)*W)))/W
;

vreal RtW23
=
(tDtDW23 + gamt23*(2*invgamt12*tDtDW12 + 2*invgamt13*tDtDW13 + 
       invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 + invgamt33*tDtDW33 - 
       4*dlnW1*dlnW2*invgamt12*W - 4*dlnW1*dlnW3*invgamt13*W - 
       2*Power(dlnW2,2)*invgamt22*W - 4*dlnW2*dlnW3*invgamt23*W - 
       2*Power(dlnW3,2)*invgamt33*W + 
       invgamt11*(tDtDW11 - 2*Power(dlnW1,2)*W)))/W
;

vreal RtW33
=
(tDtDW33 + gamt33*(2*invgamt12*tDtDW12 + 2*invgamt13*tDtDW13 + 
       invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 + invgamt33*tDtDW33 - 
       4*dlnW1*dlnW2*invgamt12*W - 4*dlnW1*dlnW3*invgamt13*W - 
       2*Power(dlnW2,2)*invgamt22*W - 4*dlnW2*dlnW3*invgamt23*W - 
       2*Power(dlnW3,2)*invgamt33*W + 
       invgamt11*(tDtDW11 - 2*Power(dlnW1,2)*W)))/W
;

vreal Rt11
=
dtrGt11*gamt11 + dtrGt12*gamt12 + dtrGt13*gamt13 + 3*Gt111*GtDDU111 + 
  3*Gt112*GtDDU112 + 3*Gt113*GtDDU113 + 2*Gt211*GtDDU121 + 
  2*Gt212*GtDDU122 + 2*Gt213*GtDDU123 + 2*Gt311*GtDDU131 + 
  2*Gt312*GtDDU132 + 2*Gt313*GtDDU133 + Gt211*GtDDU211 + Gt212*GtDDU212 + 
  Gt213*GtDDU213 + Gt311*GtDDU311 + Gt312*GtDDU312 + Gt313*GtDDU313 - 
  (ddgamt1111*invgamt11)/2. - ddgamt1211*invgamt12 - ddgamt1311*invgamt13 - 
  (ddgamt2211*invgamt22)/2. - ddgamt2311*invgamt23 - 
  (ddgamt3311*invgamt33)/2. + GtDDD111*trGtd1 + GtDDD112*trGtd2 + 
  GtDDD113*trGtd3
;

vreal Rt12
=
(dtrGt21*gamt11 + dtrGt11*gamt12 + dtrGt22*gamt12 + dtrGt23*gamt13 + 
    dtrGt12*gamt22 + dtrGt13*gamt23 + 2*Gt112*GtDDU111 + 2*Gt122*GtDDU112 + 
    2*Gt123*GtDDU113 + 2*Gt111*GtDDU121 + 2*Gt212*GtDDU121 + 
    2*Gt112*GtDDU122 + 2*Gt222*GtDDU122 + 2*Gt113*GtDDU123 + 
    2*Gt223*GtDDU123 + 2*Gt312*GtDDU131 + 2*Gt322*GtDDU132 + 
    2*Gt323*GtDDU133 + 2*Gt111*GtDDU211 + 2*Gt112*GtDDU212 + 
    2*Gt113*GtDDU213 + 4*Gt211*GtDDU221 + 4*Gt212*GtDDU222 + 
    4*Gt213*GtDDU223 + 2*Gt311*GtDDU231 + 2*Gt312*GtDDU232 + 
    2*Gt313*GtDDU233 + 2*Gt311*GtDDU321 + 2*Gt312*GtDDU322 + 
    2*Gt313*GtDDU323 - ddgamt1112*invgamt11 - 2*ddgamt1212*invgamt12 - 
    2*ddgamt1312*invgamt13 - ddgamt2212*invgamt22 - 
    2*ddgamt2312*invgamt23 - ddgamt3312*invgamt33 + GtDDD112*trGtd1 + 
    GtDDD211*trGtd1 + GtDDD122*trGtd2 + GtDDD212*trGtd2 + GtDDD123*trGtd3 + 
    GtDDD213*trGtd3)/2.
;

vreal Rt13
=
(dtrGt31*gamt11 + dtrGt32*gamt12 + dtrGt11*gamt13 + dtrGt33*gamt13 + 
    dtrGt12*gamt23 + dtrGt13*gamt33 + 2*Gt113*GtDDU111 + 2*Gt123*GtDDU112 + 
    2*Gt133*GtDDU113 + 2*Gt213*GtDDU121 + 2*Gt223*GtDDU122 + 
    2*Gt233*GtDDU123 + 2*Gt111*GtDDU131 + 2*Gt313*GtDDU131 + 
    2*Gt112*GtDDU132 + 2*Gt323*GtDDU132 + 2*Gt113*GtDDU133 + 
    2*Gt333*GtDDU133 + 2*Gt211*GtDDU231 + 2*Gt212*GtDDU232 + 
    2*Gt213*GtDDU233 + 2*Gt111*GtDDU311 + 2*Gt112*GtDDU312 + 
    2*Gt113*GtDDU313 + 2*Gt211*GtDDU321 + 2*Gt212*GtDDU322 + 
    2*Gt213*GtDDU323 + 4*Gt311*GtDDU331 + 4*Gt312*GtDDU332 + 
    4*Gt313*GtDDU333 - ddgamt1113*invgamt11 - 2*ddgamt1213*invgamt12 - 
    2*ddgamt1313*invgamt13 - ddgamt2213*invgamt22 - 
    2*ddgamt2313*invgamt23 - ddgamt3313*invgamt33 + GtDDD113*trGtd1 + 
    GtDDD311*trGtd1 + GtDDD123*trGtd2 + GtDDD312*trGtd2 + GtDDD133*trGtd3 + 
    GtDDD313*trGtd3)/2.
;

vreal Rt22
=
dtrGt21*gamt12 + dtrGt22*gamt22 + dtrGt23*gamt23 + Gt112*GtDDU121 + 
  Gt122*GtDDU122 + Gt123*GtDDU123 + 2*Gt112*GtDDU211 + 2*Gt122*GtDDU212 + 
  2*Gt123*GtDDU213 + 3*Gt212*GtDDU221 + 3*Gt222*GtDDU222 + 
  3*Gt223*GtDDU223 + 2*Gt312*GtDDU231 + 2*Gt322*GtDDU232 + 
  2*Gt323*GtDDU233 + Gt312*GtDDU321 + Gt322*GtDDU322 + Gt323*GtDDU323 - 
  (ddgamt1122*invgamt11)/2. - ddgamt1222*invgamt12 - ddgamt1322*invgamt13 - 
  (ddgamt2222*invgamt22)/2. - ddgamt2322*invgamt23 - 
  (ddgamt3322*invgamt33)/2. + GtDDD212*trGtd1 + GtDDD222*trGtd2 + 
  GtDDD223*trGtd3
;

vreal Rt23
=
(dtrGt31*gamt12 + dtrGt21*gamt13 + dtrGt32*gamt22 + dtrGt22*gamt23 + 
    dtrGt33*gamt23 + dtrGt23*gamt33 + 2*Gt112*GtDDU131 + 2*Gt122*GtDDU132 + 
    2*Gt123*GtDDU133 + 2*Gt113*GtDDU211 + 2*Gt123*GtDDU212 + 
    2*Gt133*GtDDU213 + 2*Gt213*GtDDU221 + 2*Gt223*GtDDU222 + 
    2*Gt233*GtDDU223 + 2*Gt212*GtDDU231 + 2*Gt313*GtDDU231 + 
    2*Gt222*GtDDU232 + 2*Gt323*GtDDU232 + 2*Gt223*GtDDU233 + 
    2*Gt333*GtDDU233 + 2*Gt112*GtDDU311 + 2*Gt122*GtDDU312 + 
    2*Gt123*GtDDU313 + 2*Gt212*GtDDU321 + 2*Gt222*GtDDU322 + 
    2*Gt223*GtDDU323 + 4*Gt312*GtDDU331 + 4*Gt322*GtDDU332 + 
    4*Gt323*GtDDU333 - ddgamt1123*invgamt11 - 2*ddgamt1223*invgamt12 - 
    2*ddgamt1323*invgamt13 - ddgamt2223*invgamt22 - 
    2*ddgamt2323*invgamt23 - ddgamt3323*invgamt33 + GtDDD213*trGtd1 + 
    GtDDD312*trGtd1 + GtDDD223*trGtd2 + GtDDD322*trGtd2 + GtDDD233*trGtd3 + 
    GtDDD323*trGtd3)/2.
;

vreal Rt33
=
dtrGt31*gamt13 + dtrGt32*gamt23 + dtrGt33*gamt33 + Gt113*GtDDU131 + 
  Gt123*GtDDU132 + Gt133*GtDDU133 + Gt213*GtDDU231 + Gt223*GtDDU232 + 
  Gt233*GtDDU233 + 2*Gt113*GtDDU311 + 2*Gt123*GtDDU312 + 2*Gt133*GtDDU313 + 
  2*Gt213*GtDDU321 + 2*Gt223*GtDDU322 + 2*Gt233*GtDDU323 + 
  3*Gt313*GtDDU331 + 3*Gt323*GtDDU332 + 3*Gt333*GtDDU333 - 
  (ddgamt1133*invgamt11)/2. - ddgamt1233*invgamt12 - ddgamt1333*invgamt13 - 
  (ddgamt2233*invgamt22)/2. - ddgamt2333*invgamt23 - 
  (ddgamt3333*invgamt33)/2. + GtDDD313*trGtd1 + GtDDD323*trGtd2 + 
  GtDDD333*trGtd3
;

vreal R11
=
Rt11 + RtW11
;

vreal R12
=
Rt12 + RtW12
;

vreal R13
=
Rt13 + RtW13
;

vreal R22
=
Rt22 + RtW22
;

vreal R23
=
Rt23 + RtW23
;

vreal R33
=
Rt33 + RtW33
;

vreal trR
=
invgam11*R11 + 2*invgam12*R12 + 2*invgam13*R13 + invgam22*R22 + 
  2*invgam23*R23 + invgam33*R33
;

vreal rho
=
(Power(beta1,2)*eT11 + Power(beta2,2)*eT22 + 2*beta2*beta3*eT23 + 
    Power(beta3,2)*eT33 + 2*beta1*(beta2*eT12 + beta3*eT13 - eTt1) - 
    2*beta2*eTt2 - 2*beta3*eTt3 + eTtt)/Power(alpha,2)
;

vreal Sm1
=
(beta1*eT11 + beta2*eT12 + beta3*eT13 - eTt1)/alpha
;

vreal Sm2
=
(beta1*eT12 + beta2*eT22 + beta3*eT23 - eTt2)/alpha
;

vreal Sm3
=
(beta1*eT13 + beta2*eT23 + beta3*eT33 - eTt3)/alpha
;

vreal trdexAtUU1
=
-2*dgamt111*exAtUU11*invgamt11 - dgamt211*exAtUU12*invgamt11 - 
  2*dgamt113*exAtUU13*invgamt11 - dgamt311*exAtUU13*invgamt11 - 
  dgamt212*exAtUU22*invgamt11 - dgamt213*exAtUU23*invgamt11 - 
  dgamt312*exAtUU23*invgamt11 - dgamt313*exAtUU33*invgamt11 + 
  dexAt111*Power(invgamt11,2) - dgamt211*exAtUU11*invgamt12 - 
  2*dgamt122*exAtUU12*invgamt12 - 2*dgamt212*exAtUU12*invgamt12 - 
  2*dgamt123*exAtUU13*invgamt12 - dgamt213*exAtUU13*invgamt12 - 
  dgamt312*exAtUU13*invgamt12 - dgamt222*exAtUU22*invgamt12 - 
  dgamt223*exAtUU23*invgamt12 - dgamt322*exAtUU23*invgamt12 - 
  dgamt323*exAtUU33*invgamt12 + 2*dexAt112*invgamt11*invgamt12 + 
  dexAt211*invgamt11*invgamt12 + dexAt122*Power(invgamt12,2) + 
  dexAt212*Power(invgamt12,2) - 
  2*dgamt112*(exAtUU12*invgamt11 + exAtUU11*invgamt12) - 
  2*dgamt113*exAtUU11*invgamt13 - dgamt311*exAtUU11*invgamt13 - 
  2*dgamt123*exAtUU12*invgamt13 - dgamt213*exAtUU12*invgamt13 - 
  dgamt312*exAtUU12*invgamt13 - 2*dgamt133*exAtUU13*invgamt13 - 
  2*dgamt313*exAtUU13*invgamt13 - dgamt223*exAtUU22*invgamt13 - 
  dgamt233*exAtUU23*invgamt13 - dgamt323*exAtUU23*invgamt13 - 
  dgamt333*exAtUU33*invgamt13 + 2*dexAt113*invgamt11*invgamt13 + 
  dexAt311*invgamt11*invgamt13 + 2*dexAt123*invgamt12*invgamt13 + 
  dexAt213*invgamt12*invgamt13 + dexAt312*invgamt12*invgamt13 + 
  dexAt133*Power(invgamt13,2) + dexAt313*Power(invgamt13,2) - 
  dgamt212*exAtUU11*invgamt22 - dgamt222*exAtUU12*invgamt22 - 
  dgamt223*exAtUU13*invgamt22 + dexAt212*invgamt11*invgamt22 + 
  dexAt222*invgamt12*invgamt22 + dexAt223*invgamt13*invgamt22 - 
  dgamt213*exAtUU11*invgamt23 - dgamt312*exAtUU11*invgamt23 - 
  dgamt223*exAtUU12*invgamt23 - dgamt322*exAtUU12*invgamt23 - 
  dgamt233*exAtUU13*invgamt23 - dgamt323*exAtUU13*invgamt23 + 
  dexAt213*invgamt11*invgamt23 + dexAt312*invgamt11*invgamt23 + 
  dexAt223*invgamt12*invgamt23 + dexAt322*invgamt12*invgamt23 + 
  dexAt233*invgamt13*invgamt23 + dexAt323*invgamt13*invgamt23 - 
  dgamt313*exAtUU11*invgamt33 - dgamt323*exAtUU12*invgamt33 - 
  dgamt333*exAtUU13*invgamt33 + dexAt313*invgamt11*invgamt33 + 
  dexAt323*invgamt12*invgamt33 + dexAt333*invgamt13*invgamt33
;

vreal trdexAtUU2
=
-(dgamt113*exAtUU23*invgamt11) - 2*dgamt211*exAtUU12*invgamt12 - 
  dgamt113*exAtUU13*invgamt12 - dgamt311*exAtUU13*invgamt12 - 
  dgamt122*exAtUU22*invgamt12 - 2*dgamt212*exAtUU22*invgamt12 - 
  dgamt123*exAtUU23*invgamt12 - 2*dgamt213*exAtUU23*invgamt12 - 
  dgamt312*exAtUU23*invgamt12 - dgamt313*exAtUU33*invgamt12 + 
  dexAt111*invgamt11*invgamt12 + dexAt112*Power(invgamt12,2) + 
  dexAt211*Power(invgamt12,2) - 
  dgamt111*(exAtUU12*invgamt11 + exAtUU11*invgamt12) - 
  dgamt113*exAtUU12*invgamt13 - dgamt311*exAtUU12*invgamt13 - 
  dgamt123*exAtUU22*invgamt13 - dgamt312*exAtUU22*invgamt13 - 
  dgamt133*exAtUU23*invgamt13 - dgamt313*exAtUU23*invgamt13 + 
  dexAt113*invgamt12*invgamt13 + dexAt311*invgamt12*invgamt13 - 
  dgamt122*exAtUU12*invgamt22 - 2*dgamt212*exAtUU12*invgamt22 - 
  dgamt123*exAtUU13*invgamt22 - dgamt312*exAtUU13*invgamt22 - 
  2*dgamt222*exAtUU22*invgamt22 - 2*dgamt223*exAtUU23*invgamt22 - 
  dgamt322*exAtUU23*invgamt22 - dgamt323*exAtUU33*invgamt22 + 
  dexAt112*invgamt11*invgamt22 + dexAt122*invgamt12*invgamt22 + 
  2*dexAt212*invgamt12*invgamt22 + dexAt123*invgamt13*invgamt22 + 
  dexAt312*invgamt13*invgamt22 + dexAt222*Power(invgamt22,2) - 
  dgamt112*(exAtUU22*invgamt11 + 2*exAtUU12*invgamt12 + 
     exAtUU11*invgamt22) - dgamt113*exAtUU11*invgamt23 - 
  dgamt123*exAtUU12*invgamt23 - 2*dgamt213*exAtUU12*invgamt23 - 
  dgamt312*exAtUU12*invgamt23 - dgamt133*exAtUU13*invgamt23 - 
  dgamt313*exAtUU13*invgamt23 - 2*dgamt223*exAtUU22*invgamt23 - 
  dgamt322*exAtUU22*invgamt23 - 2*dgamt233*exAtUU23*invgamt23 - 
  2*dgamt323*exAtUU23*invgamt23 - dgamt333*exAtUU33*invgamt23 + 
  dexAt113*invgamt11*invgamt23 + dexAt123*invgamt12*invgamt23 + 
  2*dexAt213*invgamt12*invgamt23 + dexAt312*invgamt12*invgamt23 + 
  dexAt133*invgamt13*invgamt23 + dexAt313*invgamt13*invgamt23 + 
  2*dexAt223*invgamt22*invgamt23 + dexAt322*invgamt22*invgamt23 + 
  dexAt233*Power(invgamt23,2) + dexAt323*Power(invgamt23,2) - 
  dgamt313*exAtUU12*invgamt33 - dgamt323*exAtUU22*invgamt33 - 
  dgamt333*exAtUU23*invgamt33 + dexAt313*invgamt12*invgamt33 + 
  dexAt323*invgamt22*invgamt33 + dexAt333*invgamt23*invgamt33
;

vreal trdexAtUU3
=
-(dgamt113*exAtUU33*invgamt11) - dgamt211*exAtUU13*invgamt12 - 
  dgamt122*exAtUU23*invgamt12 - dgamt212*exAtUU23*invgamt12 - 
  dgamt123*exAtUU33*invgamt12 - dgamt213*exAtUU33*invgamt12 - 
  dgamt211*exAtUU12*invgamt13 - 2*dgamt113*exAtUU13*invgamt13 - 
  2*dgamt311*exAtUU13*invgamt13 - dgamt212*exAtUU22*invgamt13 - 
  dgamt123*exAtUU23*invgamt13 - dgamt213*exAtUU23*invgamt13 - 
  2*dgamt312*exAtUU23*invgamt13 - dgamt133*exAtUU33*invgamt13 - 
  2*dgamt313*exAtUU33*invgamt13 + dexAt111*invgamt11*invgamt13 + 
  dexAt112*invgamt12*invgamt13 + dexAt211*invgamt12*invgamt13 + 
  dexAt113*Power(invgamt13,2) + dexAt311*Power(invgamt13,2) - 
  dgamt111*(exAtUU13*invgamt11 + exAtUU11*invgamt13) - 
  dgamt212*exAtUU13*invgamt22 - dgamt222*exAtUU23*invgamt22 - 
  dgamt223*exAtUU33*invgamt22 + dexAt212*invgamt13*invgamt22 - 
  dgamt122*exAtUU12*invgamt23 - dgamt212*exAtUU12*invgamt23 - 
  dgamt123*exAtUU13*invgamt23 - dgamt213*exAtUU13*invgamt23 - 
  2*dgamt312*exAtUU13*invgamt23 - dgamt222*exAtUU22*invgamt23 - 
  2*dgamt223*exAtUU23*invgamt23 - 2*dgamt322*exAtUU23*invgamt23 - 
  dgamt233*exAtUU33*invgamt23 - 2*dgamt323*exAtUU33*invgamt23 + 
  dexAt112*invgamt11*invgamt23 + dexAt122*invgamt12*invgamt23 + 
  dexAt212*invgamt12*invgamt23 + dexAt123*invgamt13*invgamt23 + 
  dexAt213*invgamt13*invgamt23 + 2*dexAt312*invgamt13*invgamt23 + 
  dexAt222*invgamt22*invgamt23 + dexAt223*Power(invgamt23,2) + 
  dexAt322*Power(invgamt23,2) - 
  dgamt112*(exAtUU23*invgamt11 + exAtUU13*invgamt12 + exAtUU12*invgamt13 + 
     exAtUU11*invgamt23) - dgamt113*exAtUU11*invgamt33 - 
  dgamt123*exAtUU12*invgamt33 - dgamt213*exAtUU12*invgamt33 - 
  dgamt133*exAtUU13*invgamt33 - 2*dgamt313*exAtUU13*invgamt33 - 
  dgamt223*exAtUU22*invgamt33 - dgamt233*exAtUU23*invgamt33 - 
  2*dgamt323*exAtUU23*invgamt33 - 2*dgamt333*exAtUU33*invgamt33 + 
  dexAt113*invgamt11*invgamt33 + dexAt123*invgamt12*invgamt33 + 
  dexAt213*invgamt12*invgamt33 + dexAt133*invgamt13*invgamt33 + 
  2*dexAt313*invgamt13*invgamt33 + dexAt223*invgamt22*invgamt33 + 
  dexAt233*invgamt23*invgamt33 + 2*dexAt323*invgamt23*invgamt33 + 
  dexAt333*Power(invgamt33,2)
;


local_ZtC1.store(mask, index2, 
(trGt1 - trGtd1)/2.
);

local_ZtC2.store(mask, index2, 
(trGt2 - trGtd2)/2.
);

local_ZtC3.store(mask, index2, 
(trGt3 - trGtd3)/2.
);

local_HC.store(mask, index2, 
exAt11*exAtUU11 + 2*exAt12*exAtUU12 + 2*exAt13*exAtUU13 + exAt22*exAtUU22 + 
  2*exAt23*exAtUU23 + exAt33*exAtUU33 - (2*Power(exKh,2))/3. - 16*cpi*rho - 
  (8*exKh*Theta)/3. - (8*Power(Theta,2))/3. + trR
);

local_MtC1.store(mask, index2, 
(-4*dlnW1*exAtUU11 - 4*dlnW2*exAtUU12 - 4*dlnW3*exAtUU13 + 
    3*exAtUU11*Gt111 + 6*exAtUU12*Gt112 + 6*exAtUU13*Gt113 + 
    3*exAtUU22*Gt122 + 6*exAtUU23*Gt123 + 3*exAtUU33*Gt133 - 
    2*dexKh1*invgamt11 - 4*dTheta1*invgamt11 - 2*dexKh2*invgamt12 - 
    4*dTheta2*invgamt12 - 2*dexKh3*invgamt13 - 4*dTheta3*invgamt13 - 
    24*cpi*invgamt11*Sm1 - 24*cpi*invgamt12*Sm2 - 24*cpi*invgamt13*Sm3 + 
    3*trdexAtUU1)/3.
);

local_MtC2.store(mask, index2, 
(-4*dlnW1*exAtUU12 - 4*dlnW2*exAtUU22 - 4*dlnW3*exAtUU23 + 
    3*exAtUU11*Gt211 + 6*exAtUU12*Gt212 + 6*exAtUU13*Gt213 + 
    3*exAtUU22*Gt222 + 6*exAtUU23*Gt223 + 3*exAtUU33*Gt233 - 
    2*dexKh1*invgamt12 - 4*dTheta1*invgamt12 - 2*dexKh2*invgamt22 - 
    4*dTheta2*invgamt22 - 2*dexKh3*invgamt23 - 4*dTheta3*invgamt23 - 
    24*cpi*invgamt12*Sm1 - 24*cpi*invgamt22*Sm2 - 24*cpi*invgamt23*Sm3 + 
    3*trdexAtUU2)/3.
);

local_MtC3.store(mask, index2, 
(-4*dlnW1*exAtUU13 - 4*dlnW2*exAtUU23 - 4*dlnW3*exAtUU33 + 
    3*exAtUU11*Gt311 + 6*exAtUU12*Gt312 + 6*exAtUU13*Gt313 + 
    3*exAtUU22*Gt322 + 6*exAtUU23*Gt323 + 3*exAtUU33*Gt333 - 
    2*dexKh1*invgamt13 - 4*dTheta1*invgamt13 - 2*dexKh2*invgamt23 - 
    4*dTheta2*invgamt23 - 2*dexKh3*invgamt33 - 4*dTheta3*invgamt33 - 
    24*cpi*invgamt13*Sm1 - 24*cpi*invgamt23*Sm2 - 24*cpi*invgamt33*Sm3 + 
    3*trdexAtUU3)/3.
);


  });
});

#endif // #ifndef Z4COW_SET_CONSTRAINT_HXX

/* Z4cow_set_constraint.hxx */
