/* Z4co_set_rhs.hxx */
/* Produced with Mathematica */

const GF3D2<CCTK_REAL> &local_dtchi = gf_dtchi;
const GF3D2<CCTK_REAL> &local_dtgamt11 = gf_dtgamt(0,0);
const GF3D2<CCTK_REAL> &local_dtgamt12 = gf_dtgamt(0,1);
const GF3D2<CCTK_REAL> &local_dtgamt13 = gf_dtgamt(0,2);
const GF3D2<CCTK_REAL> &local_dtgamt22 = gf_dtgamt(1,1);
const GF3D2<CCTK_REAL> &local_dtgamt23 = gf_dtgamt(1,2);
const GF3D2<CCTK_REAL> &local_dtgamt33 = gf_dtgamt(2,2);
const GF3D2<CCTK_REAL> &local_dtexKh = gf_dtexKh;
const GF3D2<CCTK_REAL> &local_dtexAt11 = gf_dtexAt(0,0);
const GF3D2<CCTK_REAL> &local_dtexAt12 = gf_dtexAt(0,1);
const GF3D2<CCTK_REAL> &local_dtexAt13 = gf_dtexAt(0,2);
const GF3D2<CCTK_REAL> &local_dtexAt22 = gf_dtexAt(1,1);
const GF3D2<CCTK_REAL> &local_dtexAt23 = gf_dtexAt(1,2);
const GF3D2<CCTK_REAL> &local_dtexAt33 = gf_dtexAt(2,2);
const GF3D2<CCTK_REAL> &local_dttrGt1 = gf_dttrGt(0);
const GF3D2<CCTK_REAL> &local_dttrGt2 = gf_dttrGt(1);
const GF3D2<CCTK_REAL> &local_dttrGt3 = gf_dttrGt(2);
const GF3D2<CCTK_REAL> &local_dtTheta = gf_dtTheta;
const GF3D2<CCTK_REAL> &local_dtalpha = gf_dtalpha;
const GF3D2<CCTK_REAL> &local_dtbeta1 = gf_dtbeta(0);
const GF3D2<CCTK_REAL> &local_dtbeta2 = gf_dtbeta(1);
const GF3D2<CCTK_REAL> &local_dtbeta3 = gf_dtbeta(2);

noinline([&]() __attribute__((__flatten__, __hot__)) {
  grid.loop_int_device<0, 0, 0, vsize>(
    grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
    const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
    const GF3D2index index2(layout2, p.I);
    const GF3D5index index5(layout5, p.I);

const auto &tmp_eTtt = gf_eTtt(mask, index2);
const auto &tmp_eTt = gf_eTt(mask, index2);
const auto &tmp_eT = gf_eT(mask, index2);
const auto &tmp_chi = tl_chi(mask, index5);
const auto &tmp_gamt = tl_gamt(mask, index5);
const auto &tmp_exKh = tl_exKh(mask, index5);
const auto &tmp_exAt = tl_exAt(mask, index5);
const auto &tmp_trGt = tl_trGt(mask, index5);
const auto &tmp_Theta = tl_Theta(mask, index5);
const auto &tmp_alpha = tl_alpha(mask, index5);
const auto &tmp_beta = tl_beta(mask, index5);
const auto &tmp_dchi = tl_dchi(mask, index5);
const auto &tmp_dgamt = tl_dgamt(mask, index5);
const auto &tmp_dexKh = tl_dexKh(mask, index5);
const auto &tmp_dtrGt = tl_dtrGt(mask, index5);
const auto &tmp_dTheta = tl_dTheta(mask, index5);
const auto &tmp_dalpha = tl_dalpha(mask, index5);
const auto &tmp_dbeta = tl_dbeta(mask, index5);
const auto &tmp_ddchi = tl_ddchi(mask, index5);
const auto &tmp_ddgamt = tl_ddgamt(mask, index5);
const auto &tmp_ddalpha = tl_ddalpha(mask, index5);
const auto &tmp_ddbeta = tl_ddbeta(mask, index5);

const vreal eTtt = tmp_eTtt;
const vreal eTt1 = tmp_eTt(0);
const vreal eTt2 = tmp_eTt(1);
const vreal eTt3 = tmp_eTt(2);
const vreal eT11 = tmp_eT(0,0);
const vreal eT12 = tmp_eT(0,1);
const vreal eT13 = tmp_eT(0,2);
const vreal eT21 = tmp_eT(1,0);
const vreal eT22 = tmp_eT(1,1);
const vreal eT23 = tmp_eT(1,2);
const vreal eT31 = tmp_eT(2,0);
const vreal eT32 = tmp_eT(2,1);
const vreal eT33 = tmp_eT(2,2);
const vreal chi = tmp_chi;
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
const vreal dchi1 = tmp_dchi(0);
const vreal dchi2 = tmp_dchi(1);
const vreal dchi3 = tmp_dchi(2);
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
const vreal dalpha1 = tmp_dalpha(0);
const vreal dalpha2 = tmp_dalpha(1);
const vreal dalpha3 = tmp_dalpha(2);
const vreal dbeta11 = tmp_dbeta(0)(0);
const vreal dbeta12 = tmp_dbeta(1)(0);
const vreal dbeta13 = tmp_dbeta(2)(0);
const vreal dbeta21 = tmp_dbeta(0)(1);
const vreal dbeta22 = tmp_dbeta(1)(1);
const vreal dbeta23 = tmp_dbeta(2)(1);
const vreal dbeta31 = tmp_dbeta(0)(2);
const vreal dbeta32 = tmp_dbeta(1)(2);
const vreal dbeta33 = tmp_dbeta(2)(2);
const vreal ddchi11 = tmp_ddchi(0,0);
const vreal ddchi12 = tmp_ddchi(0,1);
const vreal ddchi13 = tmp_ddchi(0,2);
const vreal ddchi22 = tmp_ddchi(1,1);
const vreal ddchi23 = tmp_ddchi(1,2);
const vreal ddchi33 = tmp_ddchi(2,2);
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
const vreal ddalpha11 = tmp_ddalpha(0,0);
const vreal ddalpha12 = tmp_ddalpha(0,1);
const vreal ddalpha13 = tmp_ddalpha(0,2);
const vreal ddalpha22 = tmp_ddalpha(1,1);
const vreal ddalpha23 = tmp_ddalpha(1,2);
const vreal ddalpha33 = tmp_ddalpha(2,2);
const vreal ddbeta111 = tmp_ddbeta(0)(0,0);
const vreal ddbeta112 = tmp_ddbeta(1)(0,0);
const vreal ddbeta113 = tmp_ddbeta(2)(0,0);
const vreal ddbeta121 = tmp_ddbeta(0)(0,1);
const vreal ddbeta122 = tmp_ddbeta(1)(0,1);
const vreal ddbeta123 = tmp_ddbeta(2)(0,1);
const vreal ddbeta131 = tmp_ddbeta(0)(0,2);
const vreal ddbeta132 = tmp_ddbeta(1)(0,2);
const vreal ddbeta133 = tmp_ddbeta(2)(0,2);
const vreal ddbeta221 = tmp_ddbeta(0)(1,1);
const vreal ddbeta222 = tmp_ddbeta(1)(1,1);
const vreal ddbeta223 = tmp_ddbeta(2)(1,1);
const vreal ddbeta231 = tmp_ddbeta(0)(1,2);
const vreal ddbeta232 = tmp_ddbeta(1)(1,2);
const vreal ddbeta233 = tmp_ddbeta(2)(1,2);
const vreal ddbeta331 = tmp_ddbeta(0)(2,2);
const vreal ddbeta332 = tmp_ddbeta(1)(2,2);
const vreal ddbeta333 = tmp_ddbeta(2)(2,2);

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

vreal invgam11
=
chi*invgamt11
;

vreal invgam12
=
chi*invgamt12
;

vreal invgam13
=
chi*invgamt13
;

vreal invgam22
=
chi*invgamt22
;

vreal invgam23
=
chi*invgamt23
;

vreal invgam33
=
chi*invgamt33
;

vreal gam11
=
gamt11/chi
;

vreal gam12
=
gamt12/chi
;

vreal gam13
=
gamt13/chi
;

vreal gam22
=
gamt22/chi
;

vreal gam23
=
gamt23/chi
;

vreal gam33
=
gamt33/chi
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

vreal dgam111
=
(chi*dgamt111 - dchi1*gamt11)/Power(chi,2)
;

vreal dgam112
=
(chi*dgamt112 - dchi1*gamt12)/Power(chi,2)
;

vreal dgam113
=
(chi*dgamt113 - dchi1*gamt13)/Power(chi,2)
;

vreal dgam122
=
(chi*dgamt122 - dchi1*gamt22)/Power(chi,2)
;

vreal dgam123
=
(chi*dgamt123 - dchi1*gamt23)/Power(chi,2)
;

vreal dgam133
=
(chi*dgamt133 - dchi1*gamt33)/Power(chi,2)
;

vreal dgam211
=
(chi*dgamt211 - dchi2*gamt11)/Power(chi,2)
;

vreal dgam212
=
(chi*dgamt212 - dchi2*gamt12)/Power(chi,2)
;

vreal dgam213
=
(chi*dgamt213 - dchi2*gamt13)/Power(chi,2)
;

vreal dgam222
=
(chi*dgamt222 - dchi2*gamt22)/Power(chi,2)
;

vreal dgam223
=
(chi*dgamt223 - dchi2*gamt23)/Power(chi,2)
;

vreal dgam233
=
(chi*dgamt233 - dchi2*gamt33)/Power(chi,2)
;

vreal dgam311
=
(chi*dgamt311 - dchi3*gamt11)/Power(chi,2)
;

vreal dgam312
=
(chi*dgamt312 - dchi3*gamt12)/Power(chi,2)
;

vreal dgam313
=
(chi*dgamt313 - dchi3*gamt13)/Power(chi,2)
;

vreal dgam322
=
(chi*dgamt322 - dchi3*gamt22)/Power(chi,2)
;

vreal dgam323
=
(chi*dgamt323 - dchi3*gamt23)/Power(chi,2)
;

vreal dgam333
=
(chi*dgamt333 - dchi3*gamt33)/Power(chi,2)
;

vreal GamDDD111
=
dgam111/2.
;

vreal GamDDD112
=
dgam211/2.
;

vreal GamDDD113
=
dgam311/2.
;

vreal GamDDD122
=
-0.5*dgam122 + dgam212
;

vreal GamDDD123
=
(-dgam123 + dgam213 + dgam312)/2.
;

vreal GamDDD133
=
-0.5*dgam133 + dgam313
;

vreal GamDDD211
=
dgam112 - dgam211/2.
;

vreal GamDDD212
=
dgam122/2.
;

vreal GamDDD213
=
(dgam123 - dgam213 + dgam312)/2.
;

vreal GamDDD222
=
dgam222/2.
;

vreal GamDDD223
=
dgam322/2.
;

vreal GamDDD233
=
-0.5*dgam233 + dgam323
;

vreal GamDDD311
=
dgam113 - dgam311/2.
;

vreal GamDDD312
=
(dgam123 + dgam213 - dgam312)/2.
;

vreal GamDDD313
=
dgam133/2.
;

vreal GamDDD322
=
dgam223 - dgam322/2.
;

vreal GamDDD323
=
dgam233/2.
;

vreal GamDDD333
=
dgam333/2.
;

vreal Gam111
=
GamDDD111*invgam11 + GamDDD211*invgam12 + GamDDD311*invgam13
;

vreal Gam112
=
GamDDD112*invgam11 + GamDDD212*invgam12 + GamDDD312*invgam13
;

vreal Gam113
=
GamDDD113*invgam11 + GamDDD213*invgam12 + GamDDD313*invgam13
;

vreal Gam122
=
GamDDD122*invgam11 + GamDDD222*invgam12 + GamDDD322*invgam13
;

vreal Gam123
=
GamDDD123*invgam11 + GamDDD223*invgam12 + GamDDD323*invgam13
;

vreal Gam133
=
GamDDD133*invgam11 + GamDDD233*invgam12 + GamDDD333*invgam13
;

vreal Gam211
=
GamDDD111*invgam12 + GamDDD211*invgam22 + GamDDD311*invgam23
;

vreal Gam212
=
GamDDD112*invgam12 + GamDDD212*invgam22 + GamDDD312*invgam23
;

vreal Gam213
=
GamDDD113*invgam12 + GamDDD213*invgam22 + GamDDD313*invgam23
;

vreal Gam222
=
GamDDD122*invgam12 + GamDDD222*invgam22 + GamDDD322*invgam23
;

vreal Gam223
=
GamDDD123*invgam12 + GamDDD223*invgam22 + GamDDD323*invgam23
;

vreal Gam233
=
GamDDD133*invgam12 + GamDDD233*invgam22 + GamDDD333*invgam23
;

vreal Gam311
=
GamDDD111*invgam13 + GamDDD211*invgam23 + GamDDD311*invgam33
;

vreal Gam312
=
GamDDD112*invgam13 + GamDDD212*invgam23 + GamDDD312*invgam33
;

vreal Gam313
=
GamDDD113*invgam13 + GamDDD213*invgam23 + GamDDD313*invgam33
;

vreal Gam322
=
GamDDD122*invgam13 + GamDDD222*invgam23 + GamDDD322*invgam33
;

vreal Gam323
=
GamDDD123*invgam13 + GamDDD223*invgam23 + GamDDD323*invgam33
;

vreal Gam333
=
GamDDD133*invgam13 + GamDDD233*invgam23 + GamDDD333*invgam33
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

vreal tDtDchi11
=
ddchi11 - dchi1*Gt111 - dchi2*Gt211 - dchi3*Gt311
;

vreal tDtDchi12
=
ddchi12 - dchi1*Gt112 - dchi2*Gt212 - dchi3*Gt312
;

vreal tDtDchi13
=
ddchi13 - dchi1*Gt113 - dchi2*Gt213 - dchi3*Gt313
;

vreal tDtDchi22
=
ddchi22 - dchi1*Gt122 - dchi2*Gt222 - dchi3*Gt322
;

vreal tDtDchi23
=
ddchi23 - dchi1*Gt123 - dchi2*Gt223 - dchi3*Gt323
;

vreal tDtDchi33
=
ddchi33 - dchi1*Gt133 - dchi2*Gt233 - dchi3*Gt333
;

vreal DDalpha11
=
ddalpha11 - dalpha1*Gam111 - dalpha2*Gam211 - dalpha3*Gam311
;

vreal DDalpha12
=
ddalpha12 - dalpha1*Gam112 - dalpha2*Gam212 - dalpha3*Gam312
;

vreal DDalpha13
=
ddalpha13 - dalpha1*Gam113 - dalpha2*Gam213 - dalpha3*Gam313
;

vreal DDalpha22
=
ddalpha22 - dalpha1*Gam122 - dalpha2*Gam222 - dalpha3*Gam322
;

vreal DDalpha23
=
ddalpha23 - dalpha1*Gam123 - dalpha2*Gam223 - dalpha3*Gam323
;

vreal DDalpha33
=
ddalpha33 - dalpha1*Gam133 - dalpha2*Gam233 - dalpha3*Gam333
;

vreal Rtchi11
=
(-(Power(dchi1,2)*(1 + 3*gamt11*invgamt11)) - 
    6*dchi1*gamt11*(dchi2*invgamt12 + dchi3*invgamt13) - 
    3*Power(dchi2,2)*gamt11*invgamt22 - 6*dchi2*dchi3*gamt11*invgamt23 - 
    3*Power(dchi3,2)*gamt11*invgamt33 + 2*chi*tDtDchi11 + 
    2*chi*gamt11*invgamt11*tDtDchi11 + 4*chi*gamt11*invgamt12*tDtDchi12 + 
    4*chi*gamt11*invgamt13*tDtDchi13 + 2*chi*gamt11*invgamt22*tDtDchi22 + 
    4*chi*gamt11*invgamt23*tDtDchi23 + 2*chi*gamt11*invgamt33*tDtDchi33)/
  (4.*Power(chi,2))
;

vreal Rtchi12
=
(-3*Power(dchi1,2)*gamt12*invgamt11 - 
    dchi1*(dchi2 + 6*dchi2*gamt12*invgamt12 + 6*dchi3*gamt12*invgamt13) - 
    3*Power(dchi2,2)*gamt12*invgamt22 - 6*dchi2*dchi3*gamt12*invgamt23 - 
    3*Power(dchi3,2)*gamt12*invgamt33 + 2*chi*gamt12*invgamt11*tDtDchi11 + 
    2*chi*tDtDchi12 + 4*chi*gamt12*invgamt12*tDtDchi12 + 
    4*chi*gamt12*invgamt13*tDtDchi13 + 2*chi*gamt12*invgamt22*tDtDchi22 + 
    4*chi*gamt12*invgamt23*tDtDchi23 + 2*chi*gamt12*invgamt33*tDtDchi33)/
  (4.*Power(chi,2))
;

vreal Rtchi13
=
(-3*Power(dchi1,2)*gamt13*invgamt11 - 
    dchi1*(dchi3 + 6*dchi2*gamt13*invgamt12 + 6*dchi3*gamt13*invgamt13) - 
    3*Power(dchi2,2)*gamt13*invgamt22 - 6*dchi2*dchi3*gamt13*invgamt23 - 
    3*Power(dchi3,2)*gamt13*invgamt33 + 2*chi*gamt13*invgamt11*tDtDchi11 + 
    4*chi*gamt13*invgamt12*tDtDchi12 + 2*chi*tDtDchi13 + 
    4*chi*gamt13*invgamt13*tDtDchi13 + 2*chi*gamt13*invgamt22*tDtDchi22 + 
    4*chi*gamt13*invgamt23*tDtDchi23 + 2*chi*gamt13*invgamt33*tDtDchi33)/
  (4.*Power(chi,2))
;

vreal Rtchi22
=
(-3*Power(dchi1,2)*gamt22*invgamt11 - 6*dchi1*dchi3*gamt22*invgamt13 - 
    Power(dchi2,2)*(1 + 3*gamt22*invgamt22) - 
    6*dchi2*gamt22*(dchi1*invgamt12 + dchi3*invgamt23) - 
    3*Power(dchi3,2)*gamt22*invgamt33 + 2*chi*gamt22*invgamt11*tDtDchi11 + 
    4*chi*gamt22*invgamt12*tDtDchi12 + 4*chi*gamt22*invgamt13*tDtDchi13 + 
    2*chi*tDtDchi22 + 2*chi*gamt22*invgamt22*tDtDchi22 + 
    4*chi*gamt22*invgamt23*tDtDchi23 + 2*chi*gamt22*invgamt33*tDtDchi33)/
  (4.*Power(chi,2))
;

vreal Rtchi23
=
(-3*Power(dchi1,2)*gamt23*invgamt11 - 6*dchi1*dchi3*gamt23*invgamt13 - 
    3*Power(dchi2,2)*gamt23*invgamt22 - 
    dchi2*(dchi3 + 6*dchi1*gamt23*invgamt12 + 6*dchi3*gamt23*invgamt23) - 
    3*Power(dchi3,2)*gamt23*invgamt33 + 2*chi*gamt23*invgamt11*tDtDchi11 + 
    4*chi*gamt23*invgamt12*tDtDchi12 + 4*chi*gamt23*invgamt13*tDtDchi13 + 
    2*chi*gamt23*invgamt22*tDtDchi22 + 2*chi*tDtDchi23 + 
    4*chi*gamt23*invgamt23*tDtDchi23 + 2*chi*gamt23*invgamt33*tDtDchi33)/
  (4.*Power(chi,2))
;

vreal Rtchi33
=
(-3*Power(dchi1,2)*gamt33*invgamt11 - 6*dchi1*dchi2*gamt33*invgamt12 - 
    3*Power(dchi2,2)*gamt33*invgamt22 - 
    6*dchi3*gamt33*(dchi1*invgamt13 + dchi2*invgamt23) - 
    Power(dchi3,2)*(1 + 3*gamt33*invgamt33) + 
    2*chi*gamt33*invgamt11*tDtDchi11 + 4*chi*gamt33*invgamt12*tDtDchi12 + 
    4*chi*gamt33*invgamt13*tDtDchi13 + 2*chi*gamt33*invgamt22*tDtDchi22 + 
    4*chi*gamt33*invgamt23*tDtDchi23 + 2*chi*tDtDchi33 + 
    2*chi*gamt33*invgamt33*tDtDchi33)/(4.*Power(chi,2))
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
Rt11 + Rtchi11
;

vreal R12
=
Rt12 + Rtchi12
;

vreal R13
=
Rt13 + Rtchi13
;

vreal R22
=
Rt22 + Rtchi22
;

vreal R23
=
Rt23 + Rtchi23
;

vreal R33
=
Rt33 + Rtchi33
;

vreal trR
=
invgam11*R11 + 2*invgam12*R12 + 2*invgam13*R13 + invgam22*R22 + 
  2*invgam23*R23 + invgam33*R33
;

vreal rho
=
(Power(beta1,2)*eT11 + beta1*beta2*(eT12 + eT21) + Power(beta2,2)*eT22 + 
    beta1*beta3*(eT13 + eT31) + beta2*beta3*(eT23 + eT32) + 
    Power(beta3,2)*eT33 - 2*beta1*eTt1 - 2*beta2*eTt2 - 2*beta3*eTt3 + eTtt)/
  Power(alpha,2)
;

vreal Sm1
=
(beta1*eT11 + beta2*eT21 + beta3*eT31 - eTt1)/alpha
;

vreal Sm2
=
(beta1*eT12 + beta2*eT22 + beta3*eT32 - eTt2)/alpha
;

vreal Sm3
=
(beta1*eT13 + beta2*eT23 + beta3*eT33 - eTt3)/alpha
;

vreal Ss11
=
eT11
;

vreal Ss12
=
eT12
;

vreal Ss13
=
eT13
;

vreal Ss22
=
eT22
;

vreal Ss23
=
eT23
;

vreal Ss33
=
eT33
;

vreal trSs
=
invgam11*Ss11 + 2*invgam12*Ss12 + 2*invgam13*Ss13 + invgam22*Ss22 + 
  2*invgam23*Ss23 + invgam33*Ss33
;


local_dtchi.store(mask, index2, 
(-2*chi*(dbeta11 + dbeta22 + dbeta33 - alpha*exKh - 2*alpha*Theta))/3.
);

local_dtgamt11.store(mask, index2, 
-2*alpha*exAt11 + 2*dbeta11*gamt11 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt11)/3. + 2*dbeta12*gamt12 + 
  2*dbeta13*gamt13
);

local_dtgamt12.store(mask, index2, 
-2*alpha*exAt12 + dbeta21*gamt11 + dbeta11*gamt12 + dbeta22*gamt12 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt12)/3. + dbeta23*gamt13 + 
  dbeta12*gamt22 + dbeta13*gamt23
);

local_dtgamt13.store(mask, index2, 
-2*alpha*exAt13 + dbeta31*gamt11 + dbeta32*gamt12 + dbeta11*gamt13 + 
  dbeta33*gamt13 - (2*(dbeta11 + dbeta22 + dbeta33)*gamt13)/3. + 
  dbeta12*gamt23 + dbeta13*gamt33
);

local_dtgamt22.store(mask, index2, 
-2*alpha*exAt22 + 2*dbeta21*gamt12 + 2*dbeta22*gamt22 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt22)/3. + 2*dbeta23*gamt23
);

local_dtgamt23.store(mask, index2, 
-2*alpha*exAt23 + dbeta31*gamt12 + dbeta21*gamt13 + dbeta32*gamt22 + 
  dbeta22*gamt23 + dbeta33*gamt23 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt23)/3. + dbeta23*gamt33
);

local_dtgamt33.store(mask, index2, 
-2*alpha*exAt33 + 2*dbeta31*gamt13 + 2*dbeta32*gamt23 + 2*dbeta33*gamt33 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt33)/3.
);

local_dtexKh.store(mask, index2, 
-(DDalpha11*invgam11) - 2*DDalpha12*invgam12 - 2*DDalpha13*invgam13 - 
  DDalpha22*invgam22 - 2*DDalpha23*invgam23 - DDalpha33*invgam33 + 
  alpha*(exAt11*exAtUU11 + 2*exAt12*exAtUU12 + 2*exAt13*exAtUU13 + 
     exAt22*exAtUU22 + 2*exAt23*exAtUU23 + exAt33*exAtUU33 + 
     Power(exKh,2)/3. + 4*cpi*rho + ckappa1*Theta - ckappa1*ckappa2*Theta + 
     (4*exKh*Theta)/3. + (4*Power(Theta,2))/3. + 4*cpi*trSs)
);

local_dtexAt11.store(mask, index2, 
(4*dbeta11*exAt11 - 2*dbeta22*exAt11 - 2*dbeta33*exAt11 + 
    6*dbeta12*exAt12 + 6*dbeta13*exAt13 + 3*alpha*exAt11*exKh - 
    6*alpha*Power(exAt11,2)*invgamt11 - 12*alpha*exAt11*exAt12*invgamt12 - 
    12*alpha*exAt11*exAt13*invgamt13 - 6*alpha*Power(exAt12,2)*invgamt22 - 
    12*alpha*exAt12*exAt13*invgamt23 - 6*alpha*Power(exAt13,2)*invgamt33 + 
    chi*(DDalpha11*(-3 + gam11*invgam11) + 2*DDalpha12*gam11*invgam12 + 
       2*DDalpha13*gam11*invgam13 + DDalpha22*gam11*invgam22 + 
       2*DDalpha23*gam11*invgam23 + DDalpha33*gam11*invgam33 + 
       3*alpha*R11 - alpha*gam11*invgam11*R11 - 
       2*alpha*gam11*invgam12*R12 - 2*alpha*gam11*invgam13*R13 - 
       alpha*gam11*invgam22*R22 - 2*alpha*gam11*invgam23*R23 - 
       alpha*gam11*invgam33*R33 - 24*alpha*cpi*Ss11 + 
       8*alpha*cpi*gam11*invgam11*Ss11 + 
       16*alpha*cpi*gam11*invgam12*Ss12 + 
       16*alpha*cpi*gam11*invgam13*Ss13 + 
       8*alpha*cpi*gam11*invgam22*Ss22 + 
       16*alpha*cpi*gam11*invgam23*Ss23 + 8*alpha*cpi*gam11*invgam33*Ss33) \
+ 6*alpha*exAt11*Theta)/3.
);

local_dtexAt12.store(mask, index2, 
(3*dbeta21*exAt11 + dbeta11*exAt12 + dbeta22*exAt12 - 2*dbeta33*exAt12 + 
    3*dbeta23*exAt13 + 3*dbeta12*exAt22 + 3*dbeta13*exAt23 + 
    3*alpha*exAt12*exKh - 6*alpha*exAt11*exAt12*invgamt11 - 
    6*alpha*Power(exAt12,2)*invgamt12 - 6*alpha*exAt11*exAt22*invgamt12 - 
    6*alpha*exAt12*exAt13*invgamt13 - 6*alpha*exAt11*exAt23*invgamt13 - 
    6*alpha*exAt12*exAt22*invgamt22 - 6*alpha*exAt13*exAt22*invgamt23 - 
    6*alpha*exAt12*exAt23*invgamt23 - 6*alpha*exAt13*exAt23*invgamt33 + 
    chi*(DDalpha11*gam12*invgam11 + DDalpha12*(-3 + 2*gam12*invgam12) + 
       2*DDalpha13*gam12*invgam13 + DDalpha22*gam12*invgam22 + 
       2*DDalpha23*gam12*invgam23 + DDalpha33*gam12*invgam33 - 
       alpha*gam12*invgam11*R11 + 3*alpha*R12 - 
       2*alpha*gam12*invgam12*R12 - 2*alpha*gam12*invgam13*R13 - 
       alpha*gam12*invgam22*R22 - 2*alpha*gam12*invgam23*R23 - 
       alpha*gam12*invgam33*R33 + 8*alpha*cpi*gam12*invgam11*Ss11 - 
       24*alpha*cpi*Ss12 + 16*alpha*cpi*gam12*invgam12*Ss12 + 
       16*alpha*cpi*gam12*invgam13*Ss13 + 
       8*alpha*cpi*gam12*invgam22*Ss22 + 
       16*alpha*cpi*gam12*invgam23*Ss23 + 8*alpha*cpi*gam12*invgam33*Ss33) \
+ 6*alpha*exAt12*Theta)/3.
);

local_dtexAt13.store(mask, index2, 
(3*dbeta31*exAt11 + 3*dbeta32*exAt12 + dbeta11*exAt13 - 2*dbeta22*exAt13 + 
    dbeta33*exAt13 + 3*dbeta12*exAt23 + 3*dbeta13*exAt33 + 
    3*alpha*exAt13*exKh - 6*alpha*exAt11*exAt13*invgamt11 - 
    6*alpha*exAt12*exAt13*invgamt12 - 6*alpha*exAt11*exAt23*invgamt12 - 
    6*alpha*Power(exAt13,2)*invgamt13 - 6*alpha*exAt11*exAt33*invgamt13 - 
    6*alpha*exAt12*exAt23*invgamt22 - 6*alpha*exAt13*exAt23*invgamt23 - 
    6*alpha*exAt12*exAt33*invgamt23 - 6*alpha*exAt13*exAt33*invgamt33 + 
    chi*(DDalpha11*gam13*invgam11 + 2*DDalpha12*gam13*invgam12 + 
       DDalpha13*(-3 + 2*gam13*invgam13) + DDalpha22*gam13*invgam22 + 
       2*DDalpha23*gam13*invgam23 + DDalpha33*gam13*invgam33 - 
       alpha*gam13*invgam11*R11 - 2*alpha*gam13*invgam12*R12 + 
       3*alpha*R13 - 2*alpha*gam13*invgam13*R13 - 
       alpha*gam13*invgam22*R22 - 2*alpha*gam13*invgam23*R23 - 
       alpha*gam13*invgam33*R33 + 8*alpha*cpi*gam13*invgam11*Ss11 + 
       16*alpha*cpi*gam13*invgam12*Ss12 - 24*alpha*cpi*Ss13 + 
       16*alpha*cpi*gam13*invgam13*Ss13 + 
       8*alpha*cpi*gam13*invgam22*Ss22 + 
       16*alpha*cpi*gam13*invgam23*Ss23 + 8*alpha*cpi*gam13*invgam33*Ss33) \
+ 6*alpha*exAt13*Theta)/3.
);

local_dtexAt22.store(mask, index2, 
(6*dbeta21*exAt12 - 2*dbeta11*exAt22 + 4*dbeta22*exAt22 - 
    2*dbeta33*exAt22 + 6*dbeta23*exAt23 + 3*alpha*exAt22*exKh - 
    6*alpha*Power(exAt12,2)*invgamt11 - 12*alpha*exAt12*exAt22*invgamt12 - 
    12*alpha*exAt12*exAt23*invgamt13 - 6*alpha*Power(exAt22,2)*invgamt22 - 
    12*alpha*exAt22*exAt23*invgamt23 - 6*alpha*Power(exAt23,2)*invgamt33 + 
    chi*(DDalpha11*gam22*invgam11 + 2*DDalpha12*gam22*invgam12 + 
       2*DDalpha13*gam22*invgam13 + DDalpha22*(-3 + gam22*invgam22) + 
       2*DDalpha23*gam22*invgam23 + DDalpha33*gam22*invgam33 - 
       alpha*gam22*invgam11*R11 - 2*alpha*gam22*invgam12*R12 - 
       2*alpha*gam22*invgam13*R13 + 3*alpha*R22 - 
       alpha*gam22*invgam22*R22 - 2*alpha*gam22*invgam23*R23 - 
       alpha*gam22*invgam33*R33 + 8*alpha*cpi*gam22*invgam11*Ss11 + 
       16*alpha*cpi*gam22*invgam12*Ss12 + 
       16*alpha*cpi*gam22*invgam13*Ss13 - 24*alpha*cpi*Ss22 + 
       8*alpha*cpi*gam22*invgam22*Ss22 + 
       16*alpha*cpi*gam22*invgam23*Ss23 + 8*alpha*cpi*gam22*invgam33*Ss33) \
+ 6*alpha*exAt22*Theta)/3.
);

local_dtexAt23.store(mask, index2, 
(3*dbeta31*exAt12 + 3*dbeta21*exAt13 + 3*dbeta32*exAt22 - 
    2*dbeta11*exAt23 + dbeta22*exAt23 + dbeta33*exAt23 + 3*dbeta23*exAt33 + 
    3*alpha*exAt23*exKh - 6*alpha*exAt12*exAt13*invgamt11 - 
    6*alpha*exAt13*exAt22*invgamt12 - 6*alpha*exAt12*exAt23*invgamt12 - 
    6*alpha*exAt13*exAt23*invgamt13 - 6*alpha*exAt12*exAt33*invgamt13 - 
    6*alpha*exAt22*exAt23*invgamt22 - 6*alpha*Power(exAt23,2)*invgamt23 - 
    6*alpha*exAt22*exAt33*invgamt23 - 6*alpha*exAt23*exAt33*invgamt33 + 
    chi*(DDalpha11*gam23*invgam11 + 2*DDalpha12*gam23*invgam12 + 
       2*DDalpha13*gam23*invgam13 + DDalpha22*gam23*invgam22 + 
       DDalpha23*(-3 + 2*gam23*invgam23) + DDalpha33*gam23*invgam33 - 
       alpha*gam23*invgam11*R11 - 2*alpha*gam23*invgam12*R12 - 
       2*alpha*gam23*invgam13*R13 - alpha*gam23*invgam22*R22 + 
       3*alpha*R23 - 2*alpha*gam23*invgam23*R23 - 
       alpha*gam23*invgam33*R33 + 8*alpha*cpi*gam23*invgam11*Ss11 + 
       16*alpha*cpi*gam23*invgam12*Ss12 + 
       16*alpha*cpi*gam23*invgam13*Ss13 + 
       8*alpha*cpi*gam23*invgam22*Ss22 - 24*alpha*cpi*Ss23 + 
       16*alpha*cpi*gam23*invgam23*Ss23 + 8*alpha*cpi*gam23*invgam33*Ss33) \
+ 6*alpha*exAt23*Theta)/3.
);

local_dtexAt33.store(mask, index2, 
(6*dbeta31*exAt13 + 6*dbeta32*exAt23 - 2*dbeta11*exAt33 - 
    2*dbeta22*exAt33 + 4*dbeta33*exAt33 + 3*alpha*exAt33*exKh - 
    6*alpha*Power(exAt13,2)*invgamt11 - 12*alpha*exAt13*exAt23*invgamt12 - 
    12*alpha*exAt13*exAt33*invgamt13 - 6*alpha*Power(exAt23,2)*invgamt22 - 
    12*alpha*exAt23*exAt33*invgamt23 - 6*alpha*Power(exAt33,2)*invgamt33 + 
    chi*(DDalpha11*gam33*invgam11 + 2*DDalpha12*gam33*invgam12 + 
       2*DDalpha13*gam33*invgam13 + DDalpha22*gam33*invgam22 + 
       2*DDalpha23*gam33*invgam23 + DDalpha33*(-3 + gam33*invgam33) - 
       alpha*gam33*invgam11*R11 - 2*alpha*gam33*invgam12*R12 - 
       2*alpha*gam33*invgam13*R13 - alpha*gam33*invgam22*R22 - 
       2*alpha*gam33*invgam23*R23 + 3*alpha*R33 - 
       alpha*gam33*invgam33*R33 + 8*alpha*cpi*gam33*invgam11*Ss11 + 
       16*alpha*cpi*gam33*invgam12*Ss12 + 
       16*alpha*cpi*gam33*invgam13*Ss13 + 
       8*alpha*cpi*gam33*invgam22*Ss22 + 
       16*alpha*cpi*gam33*invgam23*Ss23 - 24*alpha*cpi*Ss33 + 
       8*alpha*cpi*gam33*invgam33*Ss33) + 6*alpha*exAt33*Theta)/3.
);

local_dttrGt1.store(mask, index2, 
(-6*dalpha1*exAtUU11 - 6*dalpha2*exAtUU12 - 6*dalpha3*exAtUU13 - 
    (9*alpha*(dchi1*exAtUU11 + dchi2*exAtUU12 + dchi3*exAtUU13))/chi + 
    4*ddbeta111*invgamt11 + ddbeta122*invgamt11 + ddbeta133*invgamt11 + 
    7*ddbeta121*invgamt12 + ddbeta222*invgamt12 + ddbeta233*invgamt12 + 
    7*ddbeta131*invgamt13 + ddbeta232*invgamt13 + ddbeta333*invgamt13 + 
    3*ddbeta221*invgamt22 + 6*ddbeta231*invgamt23 + 3*ddbeta331*invgamt33 - 
    dbeta11*trGtd1 + 2*dbeta22*trGtd1 + 2*dbeta33*trGtd1 + 
    2*alpha*(3*exAtUU11*Gt111 + 6*exAtUU12*Gt112 + 6*exAtUU13*Gt113 + 
       3*exAtUU22*Gt122 + 6*exAtUU23*Gt123 + 3*exAtUU33*Gt133 - 
       2*dexKh1*invgamt11 - dTheta1*invgamt11 - 2*dexKh2*invgamt12 - 
       dTheta2*invgamt12 - 2*dexKh3*invgamt13 - dTheta3*invgamt13 - 
       24*cpi*invgamt11*Sm1 - 24*cpi*invgamt12*Sm2 - 
       24*cpi*invgamt13*Sm3 - 3*ckappa1*trGt1 + 3*ckappa1*trGtd1) - 
    3*dbeta21*trGtd2 - 3*dbeta31*trGtd3)/3.
);

local_dttrGt2.store(mask, index2, 
(-6*dalpha1*exAtUU12 - 6*dalpha2*exAtUU22 - 6*dalpha3*exAtUU23 - 
    (9*alpha*(dchi1*exAtUU12 + dchi2*exAtUU22 + dchi3*exAtUU23))/chi + 
    3*ddbeta112*invgamt11 + ddbeta111*invgamt12 + 7*ddbeta122*invgamt12 + 
    ddbeta133*invgamt12 + 6*ddbeta132*invgamt13 + ddbeta121*invgamt22 + 
    4*ddbeta222*invgamt22 + ddbeta233*invgamt22 + ddbeta131*invgamt23 + 
    7*ddbeta232*invgamt23 + ddbeta333*invgamt23 + 3*ddbeta332*invgamt33 - 
    3*dbeta12*trGtd1 + 2*dbeta11*trGtd2 - dbeta22*trGtd2 + 
    2*dbeta33*trGtd2 + 2*alpha*
     (3*exAtUU11*Gt211 + 6*exAtUU12*Gt212 + 6*exAtUU13*Gt213 + 
       3*exAtUU22*Gt222 + 6*exAtUU23*Gt223 + 3*exAtUU33*Gt233 - 
       2*dexKh1*invgamt12 - dTheta1*invgamt12 - 2*dexKh2*invgamt22 - 
       dTheta2*invgamt22 - 2*dexKh3*invgamt23 - dTheta3*invgamt23 - 
       24*cpi*invgamt12*Sm1 - 24*cpi*invgamt22*Sm2 - 
       24*cpi*invgamt23*Sm3 - 3*ckappa1*trGt2 + 3*ckappa1*trGtd2) - 
    3*dbeta32*trGtd3)/3.
);

local_dttrGt3.store(mask, index2, 
(-6*dalpha1*exAtUU13 - 6*dalpha2*exAtUU23 - 6*dalpha3*exAtUU33 - 
    (9*alpha*(dchi1*exAtUU13 + dchi2*exAtUU23 + dchi3*exAtUU33))/chi + 
    3*ddbeta113*invgamt11 + 6*ddbeta123*invgamt12 + ddbeta111*invgamt13 + 
    ddbeta122*invgamt13 + 7*ddbeta133*invgamt13 + 3*ddbeta223*invgamt22 + 
    ddbeta121*invgamt23 + ddbeta222*invgamt23 + 7*ddbeta233*invgamt23 + 
    ddbeta131*invgamt33 + ddbeta232*invgamt33 + 4*ddbeta333*invgamt33 - 
    3*dbeta13*trGtd1 - 3*dbeta23*trGtd2 + 2*dbeta11*trGtd3 + 
    2*dbeta22*trGtd3 - dbeta33*trGtd3 + 
    2*alpha*(3*exAtUU11*Gt311 + 6*exAtUU12*Gt312 + 6*exAtUU13*Gt313 + 
       3*exAtUU22*Gt322 + 6*exAtUU23*Gt323 + 3*exAtUU33*Gt333 - 
       2*dexKh1*invgamt13 - dTheta1*invgamt13 - 2*dexKh2*invgamt23 - 
       dTheta2*invgamt23 - 2*dexKh3*invgamt33 - dTheta3*invgamt33 - 
       24*cpi*invgamt13*Sm1 - 24*cpi*invgamt23*Sm2 - 24*cpi*invgamt33*Sm3 - 
       3*ckappa1*trGt3 + 3*ckappa1*trGtd3))/3.
);

local_dtTheta.store(mask, index2, 
-0.16666666666666666*(alpha*(3*exAt11*exAtUU11 + 6*exAt12*exAtUU12 + 
      6*exAt13*exAtUU13 + 3*exAt22*exAtUU22 + 6*exAt23*exAtUU23 + 
      3*exAt33*exAtUU33 - 2*Power(exKh,2) + 48*cpi*rho + 12*ckappa1*Theta + 
      6*ckappa1*ckappa2*Theta - 8*exKh*Theta - 8*Power(Theta,2) - 3*trR))
);

local_dtalpha.store(mask, index2, 
-(alpha*cmuL*exKh)
);

local_dtbeta1.store(mask, index2, 
-(beta1*ceta) + cmuS*trGt1
);

local_dtbeta2.store(mask, index2, 
-(beta2*ceta) + cmuS*trGt2
);

local_dtbeta3.store(mask, index2, 
-(beta3*ceta) + cmuS*trGt3
);


  });
});

/* Z4co_set_rhs.hxx */
