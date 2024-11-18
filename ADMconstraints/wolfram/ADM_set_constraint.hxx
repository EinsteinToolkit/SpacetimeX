/* ADM_set_constraint.hxx */
/* Produced with Mathematica */

#ifndef ADM_SET_CONSTRAINT_HXX
#define ADM_SET_CONSTRAINT_HXX

const GF3D2<CCTK_REAL> &local_HC = gf_HC;
const GF3D2<CCTK_REAL> &local_MC1 = gf_MC(0);
const GF3D2<CCTK_REAL> &local_MC2 = gf_MC(1);
const GF3D2<CCTK_REAL> &local_MC3 = gf_MC(2);

noinline([&]() __attribute__((__flatten__, __hot__)) {
  grid.loop_int_device<0, 0, 0, vsize>(
    grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
    const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
    const GF3D2index index2(layout2, p.I);
    const GF3D5index index5(layout5, p.I);

const auto &tmp_eTtt = gf_eTtt(mask, index2);
const auto &tmp_eTt = gf_eTt(mask, index2);
const auto &tmp_eT = gf_eT(mask, index2);
const auto &tmp_ADMgam = tl_ADMgam(mask, index5);
const auto &tmp_ADMK = tl_ADMK(mask, index5);
const auto &tmp_ADMalpha = tl_ADMalpha(mask, index5);
const auto &tmp_ADMbeta = tl_ADMbeta(mask, index5);
const auto &tmp_dADMgam = tl_dADMgam(mask, index5);
const auto &tmp_dADMK = tl_dADMK(mask, index5);
const auto &tmp_ddADMgam = tl_ddADMgam(mask, index5);

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
const vreal ADMgam11 = tmp_ADMgam(0,0);
const vreal ADMgam12 = tmp_ADMgam(0,1);
const vreal ADMgam13 = tmp_ADMgam(0,2);
const vreal ADMgam22 = tmp_ADMgam(1,1);
const vreal ADMgam23 = tmp_ADMgam(1,2);
const vreal ADMgam33 = tmp_ADMgam(2,2);
const vreal ADMK11 = tmp_ADMK(0,0);
const vreal ADMK12 = tmp_ADMK(0,1);
const vreal ADMK13 = tmp_ADMK(0,2);
const vreal ADMK22 = tmp_ADMK(1,1);
const vreal ADMK23 = tmp_ADMK(1,2);
const vreal ADMK33 = tmp_ADMK(2,2);
const vreal ADMalpha = tmp_ADMalpha;
const vreal ADMbeta1 = tmp_ADMbeta(0);
const vreal ADMbeta2 = tmp_ADMbeta(1);
const vreal ADMbeta3 = tmp_ADMbeta(2);
const vreal dADMgam111 = tmp_dADMgam(0,0)(0);
const vreal dADMgam112 = tmp_dADMgam(0,1)(0);
const vreal dADMgam113 = tmp_dADMgam(0,2)(0);
const vreal dADMgam122 = tmp_dADMgam(1,1)(0);
const vreal dADMgam123 = tmp_dADMgam(1,2)(0);
const vreal dADMgam133 = tmp_dADMgam(2,2)(0);
const vreal dADMgam211 = tmp_dADMgam(0,0)(1);
const vreal dADMgam212 = tmp_dADMgam(0,1)(1);
const vreal dADMgam213 = tmp_dADMgam(0,2)(1);
const vreal dADMgam222 = tmp_dADMgam(1,1)(1);
const vreal dADMgam223 = tmp_dADMgam(1,2)(1);
const vreal dADMgam233 = tmp_dADMgam(2,2)(1);
const vreal dADMgam311 = tmp_dADMgam(0,0)(2);
const vreal dADMgam312 = tmp_dADMgam(0,1)(2);
const vreal dADMgam313 = tmp_dADMgam(0,2)(2);
const vreal dADMgam322 = tmp_dADMgam(1,1)(2);
const vreal dADMgam323 = tmp_dADMgam(1,2)(2);
const vreal dADMgam333 = tmp_dADMgam(2,2)(2);
const vreal dADMK111 = tmp_dADMK(0,0)(0);
const vreal dADMK112 = tmp_dADMK(0,1)(0);
const vreal dADMK113 = tmp_dADMK(0,2)(0);
const vreal dADMK122 = tmp_dADMK(1,1)(0);
const vreal dADMK123 = tmp_dADMK(1,2)(0);
const vreal dADMK133 = tmp_dADMK(2,2)(0);
const vreal dADMK211 = tmp_dADMK(0,0)(1);
const vreal dADMK212 = tmp_dADMK(0,1)(1);
const vreal dADMK213 = tmp_dADMK(0,2)(1);
const vreal dADMK222 = tmp_dADMK(1,1)(1);
const vreal dADMK223 = tmp_dADMK(1,2)(1);
const vreal dADMK233 = tmp_dADMK(2,2)(1);
const vreal dADMK311 = tmp_dADMK(0,0)(2);
const vreal dADMK312 = tmp_dADMK(0,1)(2);
const vreal dADMK313 = tmp_dADMK(0,2)(2);
const vreal dADMK322 = tmp_dADMK(1,1)(2);
const vreal dADMK323 = tmp_dADMK(1,2)(2);
const vreal dADMK333 = tmp_dADMK(2,2)(2);
const vreal ddADMgam1111 = tmp_ddADMgam(0,0)(0,0);
const vreal ddADMgam1112 = tmp_ddADMgam(0,1)(0,0);
const vreal ddADMgam1113 = tmp_ddADMgam(0,2)(0,0);
const vreal ddADMgam1122 = tmp_ddADMgam(1,1)(0,0);
const vreal ddADMgam1123 = tmp_ddADMgam(1,2)(0,0);
const vreal ddADMgam1133 = tmp_ddADMgam(2,2)(0,0);
const vreal ddADMgam1211 = tmp_ddADMgam(0,0)(0,1);
const vreal ddADMgam1212 = tmp_ddADMgam(0,1)(0,1);
const vreal ddADMgam1213 = tmp_ddADMgam(0,2)(0,1);
const vreal ddADMgam1222 = tmp_ddADMgam(1,1)(0,1);
const vreal ddADMgam1223 = tmp_ddADMgam(1,2)(0,1);
const vreal ddADMgam1233 = tmp_ddADMgam(2,2)(0,1);
const vreal ddADMgam1311 = tmp_ddADMgam(0,0)(0,2);
const vreal ddADMgam1312 = tmp_ddADMgam(0,1)(0,2);
const vreal ddADMgam1313 = tmp_ddADMgam(0,2)(0,2);
const vreal ddADMgam1322 = tmp_ddADMgam(1,1)(0,2);
const vreal ddADMgam1323 = tmp_ddADMgam(1,2)(0,2);
const vreal ddADMgam1333 = tmp_ddADMgam(2,2)(0,2);
const vreal ddADMgam2211 = tmp_ddADMgam(0,0)(1,1);
const vreal ddADMgam2212 = tmp_ddADMgam(0,1)(1,1);
const vreal ddADMgam2213 = tmp_ddADMgam(0,2)(1,1);
const vreal ddADMgam2222 = tmp_ddADMgam(1,1)(1,1);
const vreal ddADMgam2223 = tmp_ddADMgam(1,2)(1,1);
const vreal ddADMgam2233 = tmp_ddADMgam(2,2)(1,1);
const vreal ddADMgam2311 = tmp_ddADMgam(0,0)(1,2);
const vreal ddADMgam2312 = tmp_ddADMgam(0,1)(1,2);
const vreal ddADMgam2313 = tmp_ddADMgam(0,2)(1,2);
const vreal ddADMgam2322 = tmp_ddADMgam(1,1)(1,2);
const vreal ddADMgam2323 = tmp_ddADMgam(1,2)(1,2);
const vreal ddADMgam2333 = tmp_ddADMgam(2,2)(1,2);
const vreal ddADMgam3311 = tmp_ddADMgam(0,0)(2,2);
const vreal ddADMgam3312 = tmp_ddADMgam(0,1)(2,2);
const vreal ddADMgam3313 = tmp_ddADMgam(0,2)(2,2);
const vreal ddADMgam3322 = tmp_ddADMgam(1,1)(2,2);
const vreal ddADMgam3323 = tmp_ddADMgam(1,2)(2,2);
const vreal ddADMgam3333 = tmp_ddADMgam(2,2)(2,2);

vreal detinvgam
=
1/(-(Power(ADMgam13,2)*ADMgam22) + 2*ADMgam12*ADMgam13*ADMgam23 -
    ADMgam11*Power(ADMgam23,2) - Power(ADMgam12,2)*ADMgam33 +
    ADMgam11*ADMgam22*ADMgam33)
;

vreal invgam11
=
(-Power(ADMgam23,2) + ADMgam22*ADMgam33)*detinvgam
;

vreal invgam12
=
(ADMgam13*ADMgam23 - ADMgam12*ADMgam33)*detinvgam
;

vreal invgam13
=
(-(ADMgam13*ADMgam22) + ADMgam12*ADMgam23)*detinvgam
;

vreal invgam22
=
(-Power(ADMgam13,2) + ADMgam11*ADMgam33)*detinvgam
;

vreal invgam23
=
(ADMgam12*ADMgam13 - ADMgam11*ADMgam23)*detinvgam
;

vreal invgam33
=
(-Power(ADMgam12,2) + ADMgam11*ADMgam22)*detinvgam
;

vreal GamDDD111
=
dADMgam111/2.
;

vreal GamDDD112
=
dADMgam211/2.
;

vreal GamDDD113
=
dADMgam311/2.
;

vreal GamDDD122
=
-0.5*dADMgam122 + dADMgam212
;

vreal GamDDD123
=
(-dADMgam123 + dADMgam213 + dADMgam312)/2.
;

vreal GamDDD133
=
-0.5*dADMgam133 + dADMgam313
;

vreal GamDDD211
=
dADMgam112 - dADMgam211/2.
;

vreal GamDDD212
=
dADMgam122/2.
;

vreal GamDDD213
=
(dADMgam123 - dADMgam213 + dADMgam312)/2.
;

vreal GamDDD222
=
dADMgam222/2.
;

vreal GamDDD223
=
dADMgam322/2.
;

vreal GamDDD233
=
-0.5*dADMgam233 + dADMgam323
;

vreal GamDDD311
=
dADMgam113 - dADMgam311/2.
;

vreal GamDDD312
=
(dADMgam123 + dADMgam213 - dADMgam312)/2.
;

vreal GamDDD313
=
dADMgam133/2.
;

vreal GamDDD322
=
dADMgam223 - dADMgam322/2.
;

vreal GamDDD323
=
dADMgam233/2.
;

vreal GamDDD333
=
dADMgam333/2.
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

vreal tr1dGam11
=
(ddADMgam1111*invgam11)/2. + ddADMgam1112*invgam12 -
  2*dADMgam112*GamDDD111*invgam11*invgam12 -
  dADMgam211*GamDDD111*invgam11*invgam12 -
  dADMgam122*GamDDD111*Power(invgam12,2) -
  dADMgam212*GamDDD111*Power(invgam12,2) -
  dADMgam112*GamDDD211*Power(invgam12,2) -
  dADMgam211*GamDDD211*Power(invgam12,2) + ddADMgam1113*invgam13 -
  2*dADMgam113*GamDDD111*invgam11*invgam13 -
  dADMgam311*GamDDD111*invgam11*invgam13 -
  2*dADMgam123*GamDDD111*invgam12*invgam13 -
  dADMgam213*GamDDD111*invgam12*invgam13 -
  dADMgam312*GamDDD111*invgam12*invgam13 -
  dADMgam113*GamDDD211*invgam12*invgam13 -
  dADMgam311*GamDDD211*invgam12*invgam13 -
  dADMgam112*GamDDD311*invgam12*invgam13 -
  dADMgam211*GamDDD311*invgam12*invgam13 -
  dADMgam133*GamDDD111*Power(invgam13,2) -
  dADMgam313*GamDDD111*Power(invgam13,2) -
  dADMgam113*GamDDD311*Power(invgam13,2) -
  dADMgam311*GamDDD311*Power(invgam13,2) -
  dADMgam111*invgam11*(GamDDD111*invgam11 + GamDDD211*invgam12 +
     GamDDD311*invgam13) + ddADMgam1212*invgam22 -
  (ddADMgam2211*invgam22)/2. - dADMgam212*GamDDD111*invgam11*invgam22 -
  dADMgam112*GamDDD211*invgam11*invgam22 -
  dADMgam222*GamDDD111*invgam12*invgam22 -
  dADMgam122*GamDDD211*invgam12*invgam22 -
  2*dADMgam212*GamDDD211*invgam12*invgam22 -
  dADMgam223*GamDDD111*invgam13*invgam22 -
  dADMgam123*GamDDD211*invgam13*invgam22 -
  dADMgam312*GamDDD211*invgam13*invgam22 -
  dADMgam212*GamDDD311*invgam13*invgam22 -
  dADMgam222*GamDDD211*Power(invgam22,2) + ddADMgam1213*invgam23 +
  ddADMgam1312*invgam23 - ddADMgam2311*invgam23 -
  dADMgam213*GamDDD111*invgam11*invgam23 -
  dADMgam312*GamDDD111*invgam11*invgam23 -
  dADMgam113*GamDDD211*invgam11*invgam23 -
  dADMgam112*GamDDD311*invgam11*invgam23 -
  dADMgam223*GamDDD111*invgam12*invgam23 -
  dADMgam322*GamDDD111*invgam12*invgam23 -
  dADMgam123*GamDDD211*invgam12*invgam23 -
  2*dADMgam213*GamDDD211*invgam12*invgam23 -
  dADMgam312*GamDDD211*invgam12*invgam23 -
  dADMgam122*GamDDD311*invgam12*invgam23 -
  dADMgam212*GamDDD311*invgam12*invgam23 -
  dADMgam233*GamDDD111*invgam13*invgam23 -
  dADMgam323*GamDDD111*invgam13*invgam23 -
  dADMgam133*GamDDD211*invgam13*invgam23 -
  dADMgam313*GamDDD211*invgam13*invgam23 -
  dADMgam123*GamDDD311*invgam13*invgam23 -
  dADMgam213*GamDDD311*invgam13*invgam23 -
  2*dADMgam312*GamDDD311*invgam13*invgam23 -
  2*dADMgam223*GamDDD211*invgam22*invgam23 -
  dADMgam322*GamDDD211*invgam22*invgam23 -
  dADMgam222*GamDDD311*invgam22*invgam23 -
  dADMgam233*GamDDD211*Power(invgam23,2) -
  dADMgam323*GamDDD211*Power(invgam23,2) -
  dADMgam223*GamDDD311*Power(invgam23,2) -
  dADMgam322*GamDDD311*Power(invgam23,2) + ddADMgam1313*invgam33 -
  (ddADMgam3311*invgam33)/2. - dADMgam313*GamDDD111*invgam11*invgam33 -
  dADMgam113*GamDDD311*invgam11*invgam33 -
  dADMgam323*GamDDD111*invgam12*invgam33 -
  dADMgam313*GamDDD211*invgam12*invgam33 -
  dADMgam123*GamDDD311*invgam12*invgam33 -
  dADMgam213*GamDDD311*invgam12*invgam33 -
  dADMgam333*GamDDD111*invgam13*invgam33 -
  dADMgam133*GamDDD311*invgam13*invgam33 -
  2*dADMgam313*GamDDD311*invgam13*invgam33 -
  dADMgam323*GamDDD211*invgam22*invgam33 -
  dADMgam223*GamDDD311*invgam22*invgam33 -
  dADMgam333*GamDDD211*invgam23*invgam33 -
  dADMgam233*GamDDD311*invgam23*invgam33 -
  2*dADMgam323*GamDDD311*invgam23*invgam33 -
  dADMgam333*GamDDD311*Power(invgam33,2)
;

vreal tr1dGam12
=
(ddADMgam1211*invgam11 + ddADMgam1122*invgam12 + ddADMgam2211*invgam12 -
    4*dADMgam112*GamDDD112*invgam11*invgam12 -
    2*dADMgam211*GamDDD112*invgam11*invgam12 -
    2*dADMgam122*GamDDD112*Power(invgam12,2) -
    2*dADMgam212*GamDDD112*Power(invgam12,2) -
    2*dADMgam112*GamDDD212*Power(invgam12,2) -
    2*dADMgam211*GamDDD212*Power(invgam12,2) + ddADMgam1123*invgam13 +
    ddADMgam1213*invgam13 - ddADMgam1312*invgam13 + ddADMgam2311*invgam13 -
    4*dADMgam113*GamDDD112*invgam11*invgam13 -
    2*dADMgam311*GamDDD112*invgam11*invgam13 -
    4*dADMgam123*GamDDD112*invgam12*invgam13 -
    2*dADMgam213*GamDDD112*invgam12*invgam13 -
    2*dADMgam312*GamDDD112*invgam12*invgam13 -
    2*dADMgam113*GamDDD212*invgam12*invgam13 -
    2*dADMgam311*GamDDD212*invgam12*invgam13 -
    2*dADMgam112*GamDDD312*invgam12*invgam13 -
    2*dADMgam211*GamDDD312*invgam12*invgam13 -
    2*dADMgam133*GamDDD112*Power(invgam13,2) -
    2*dADMgam313*GamDDD112*Power(invgam13,2) -
    2*dADMgam113*GamDDD312*Power(invgam13,2) -
    2*dADMgam311*GamDDD312*Power(invgam13,2) -
    2*dADMgam111*invgam11*(GamDDD112*invgam11 + GamDDD212*invgam12 +
       GamDDD312*invgam13) + ddADMgam1222*invgam22 -
    2*dADMgam212*GamDDD112*invgam11*invgam22 -
    2*dADMgam112*GamDDD212*invgam11*invgam22 -
    2*dADMgam222*GamDDD112*invgam12*invgam22 -
    2*dADMgam122*GamDDD212*invgam12*invgam22 -
    4*dADMgam212*GamDDD212*invgam12*invgam22 -
    2*dADMgam223*GamDDD112*invgam13*invgam22 -
    2*dADMgam123*GamDDD212*invgam13*invgam22 -
    2*dADMgam312*GamDDD212*invgam13*invgam22 -
    2*dADMgam212*GamDDD312*invgam13*invgam22 -
    2*dADMgam222*GamDDD212*Power(invgam22,2) + ddADMgam1223*invgam23 +
    ddADMgam1322*invgam23 + ddADMgam2213*invgam23 - ddADMgam2312*invgam23 -
    2*dADMgam213*GamDDD112*invgam11*invgam23 -
    2*dADMgam312*GamDDD112*invgam11*invgam23 -
    2*dADMgam113*GamDDD212*invgam11*invgam23 -
    2*dADMgam112*GamDDD312*invgam11*invgam23 -
    2*dADMgam223*GamDDD112*invgam12*invgam23 -
    2*dADMgam322*GamDDD112*invgam12*invgam23 -
    2*dADMgam123*GamDDD212*invgam12*invgam23 -
    4*dADMgam213*GamDDD212*invgam12*invgam23 -
    2*dADMgam312*GamDDD212*invgam12*invgam23 -
    2*dADMgam122*GamDDD312*invgam12*invgam23 -
    2*dADMgam212*GamDDD312*invgam12*invgam23 -
    2*dADMgam233*GamDDD112*invgam13*invgam23 -
    2*dADMgam323*GamDDD112*invgam13*invgam23 -
    2*dADMgam133*GamDDD212*invgam13*invgam23 -
    2*dADMgam313*GamDDD212*invgam13*invgam23 -
    2*dADMgam123*GamDDD312*invgam13*invgam23 -
    2*dADMgam213*GamDDD312*invgam13*invgam23 -
    4*dADMgam312*GamDDD312*invgam13*invgam23 -
    4*dADMgam223*GamDDD212*invgam22*invgam23 -
    2*dADMgam322*GamDDD212*invgam22*invgam23 -
    2*dADMgam222*GamDDD312*invgam22*invgam23 -
    2*dADMgam233*GamDDD212*Power(invgam23,2) -
    2*dADMgam323*GamDDD212*Power(invgam23,2) -
    2*dADMgam223*GamDDD312*Power(invgam23,2) -
    2*dADMgam322*GamDDD312*Power(invgam23,2) + ddADMgam1323*invgam33 +
    ddADMgam2313*invgam33 - ddADMgam3312*invgam33 -
    2*dADMgam313*GamDDD112*invgam11*invgam33 -
    2*dADMgam113*GamDDD312*invgam11*invgam33 -
    2*dADMgam323*GamDDD112*invgam12*invgam33 -
    2*dADMgam313*GamDDD212*invgam12*invgam33 -
    2*dADMgam123*GamDDD312*invgam12*invgam33 -
    2*dADMgam213*GamDDD312*invgam12*invgam33 -
    2*dADMgam333*GamDDD112*invgam13*invgam33 -
    2*dADMgam133*GamDDD312*invgam13*invgam33 -
    4*dADMgam313*GamDDD312*invgam13*invgam33 -
    2*dADMgam323*GamDDD212*invgam22*invgam33 -
    2*dADMgam223*GamDDD312*invgam22*invgam33 -
    2*dADMgam333*GamDDD212*invgam23*invgam33 -
    2*dADMgam233*GamDDD312*invgam23*invgam33 -
    4*dADMgam323*GamDDD312*invgam23*invgam33 -
    2*dADMgam333*GamDDD312*Power(invgam33,2))/2.
;

vreal tr1dGam13
=
(ddADMgam1311*invgam11 + ddADMgam1123*invgam12 - ddADMgam1213*invgam12 +
    ddADMgam1312*invgam12 + ddADMgam2311*invgam12 -
    4*dADMgam112*GamDDD113*invgam11*invgam12 -
    2*dADMgam211*GamDDD113*invgam11*invgam12 -
    2*dADMgam122*GamDDD113*Power(invgam12,2) -
    2*dADMgam212*GamDDD113*Power(invgam12,2) -
    2*dADMgam112*GamDDD213*Power(invgam12,2) -
    2*dADMgam211*GamDDD213*Power(invgam12,2) + ddADMgam1133*invgam13 +
    ddADMgam3311*invgam13 - 4*dADMgam113*GamDDD113*invgam11*invgam13 -
    2*dADMgam311*GamDDD113*invgam11*invgam13 -
    4*dADMgam123*GamDDD113*invgam12*invgam13 -
    2*dADMgam213*GamDDD113*invgam12*invgam13 -
    2*dADMgam312*GamDDD113*invgam12*invgam13 -
    2*dADMgam113*GamDDD213*invgam12*invgam13 -
    2*dADMgam311*GamDDD213*invgam12*invgam13 -
    2*dADMgam112*GamDDD313*invgam12*invgam13 -
    2*dADMgam211*GamDDD313*invgam12*invgam13 -
    2*dADMgam133*GamDDD113*Power(invgam13,2) -
    2*dADMgam313*GamDDD113*Power(invgam13,2) -
    2*dADMgam113*GamDDD313*Power(invgam13,2) -
    2*dADMgam311*GamDDD313*Power(invgam13,2) -
    2*dADMgam111*invgam11*(GamDDD113*invgam11 + GamDDD213*invgam12 +
       GamDDD313*invgam13) + ddADMgam1223*invgam22 -
    ddADMgam2213*invgam22 + ddADMgam2312*invgam22 -
    2*dADMgam212*GamDDD113*invgam11*invgam22 -
    2*dADMgam112*GamDDD213*invgam11*invgam22 -
    2*dADMgam222*GamDDD113*invgam12*invgam22 -
    2*dADMgam122*GamDDD213*invgam12*invgam22 -
    4*dADMgam212*GamDDD213*invgam12*invgam22 -
    2*dADMgam223*GamDDD113*invgam13*invgam22 -
    2*dADMgam123*GamDDD213*invgam13*invgam22 -
    2*dADMgam312*GamDDD213*invgam13*invgam22 -
    2*dADMgam212*GamDDD313*invgam13*invgam22 -
    2*dADMgam222*GamDDD213*Power(invgam22,2) + ddADMgam1233*invgam23 +
    ddADMgam1323*invgam23 - ddADMgam2313*invgam23 + ddADMgam3312*invgam23 -
    2*dADMgam213*GamDDD113*invgam11*invgam23 -
    2*dADMgam312*GamDDD113*invgam11*invgam23 -
    2*dADMgam113*GamDDD213*invgam11*invgam23 -
    2*dADMgam112*GamDDD313*invgam11*invgam23 -
    2*dADMgam223*GamDDD113*invgam12*invgam23 -
    2*dADMgam322*GamDDD113*invgam12*invgam23 -
    2*dADMgam123*GamDDD213*invgam12*invgam23 -
    4*dADMgam213*GamDDD213*invgam12*invgam23 -
    2*dADMgam312*GamDDD213*invgam12*invgam23 -
    2*dADMgam122*GamDDD313*invgam12*invgam23 -
    2*dADMgam212*GamDDD313*invgam12*invgam23 -
    2*dADMgam233*GamDDD113*invgam13*invgam23 -
    2*dADMgam323*GamDDD113*invgam13*invgam23 -
    2*dADMgam133*GamDDD213*invgam13*invgam23 -
    2*dADMgam313*GamDDD213*invgam13*invgam23 -
    2*dADMgam123*GamDDD313*invgam13*invgam23 -
    2*dADMgam213*GamDDD313*invgam13*invgam23 -
    4*dADMgam312*GamDDD313*invgam13*invgam23 -
    4*dADMgam223*GamDDD213*invgam22*invgam23 -
    2*dADMgam322*GamDDD213*invgam22*invgam23 -
    2*dADMgam222*GamDDD313*invgam22*invgam23 -
    2*dADMgam233*GamDDD213*Power(invgam23,2) -
    2*dADMgam323*GamDDD213*Power(invgam23,2) -
    2*dADMgam223*GamDDD313*Power(invgam23,2) -
    2*dADMgam322*GamDDD313*Power(invgam23,2) + ddADMgam1333*invgam33 -
    2*dADMgam313*GamDDD113*invgam11*invgam33 -
    2*dADMgam113*GamDDD313*invgam11*invgam33 -
    2*dADMgam323*GamDDD113*invgam12*invgam33 -
    2*dADMgam313*GamDDD213*invgam12*invgam33 -
    2*dADMgam123*GamDDD313*invgam12*invgam33 -
    2*dADMgam213*GamDDD313*invgam12*invgam33 -
    2*dADMgam333*GamDDD113*invgam13*invgam33 -
    2*dADMgam133*GamDDD313*invgam13*invgam33 -
    4*dADMgam313*GamDDD313*invgam13*invgam33 -
    2*dADMgam323*GamDDD213*invgam22*invgam33 -
    2*dADMgam223*GamDDD313*invgam22*invgam33 -
    2*dADMgam333*GamDDD213*invgam23*invgam33 -
    2*dADMgam233*GamDDD313*invgam23*invgam33 -
    4*dADMgam323*GamDDD313*invgam23*invgam33 -
    2*dADMgam333*GamDDD313*Power(invgam33,2))/2.
;

vreal tr1dGam22
=
-0.5*(ddADMgam1122*invgam11) + ddADMgam1212*invgam11 -
  dADMgam111*GamDDD122*Power(invgam11,2) + ddADMgam2212*invgam12 -
  2*dADMgam112*GamDDD122*invgam11*invgam12 -
  dADMgam211*GamDDD122*invgam11*invgam12 -
  dADMgam111*GamDDD222*invgam11*invgam12 -
  dADMgam122*GamDDD122*Power(invgam12,2) -
  dADMgam212*GamDDD122*Power(invgam12,2) -
  dADMgam112*GamDDD222*Power(invgam12,2) -
  dADMgam211*GamDDD222*Power(invgam12,2) + ddADMgam1223*invgam13 -
  ddADMgam1322*invgam13 + ddADMgam2312*invgam13 -
  2*dADMgam113*GamDDD122*invgam11*invgam13 -
  dADMgam311*GamDDD122*invgam11*invgam13 -
  dADMgam111*GamDDD322*invgam11*invgam13 -
  2*dADMgam123*GamDDD122*invgam12*invgam13 -
  dADMgam213*GamDDD122*invgam12*invgam13 -
  dADMgam312*GamDDD122*invgam12*invgam13 -
  dADMgam113*GamDDD222*invgam12*invgam13 -
  dADMgam311*GamDDD222*invgam12*invgam13 -
  dADMgam112*GamDDD322*invgam12*invgam13 -
  dADMgam211*GamDDD322*invgam12*invgam13 -
  dADMgam133*GamDDD122*Power(invgam13,2) -
  dADMgam313*GamDDD122*Power(invgam13,2) -
  dADMgam113*GamDDD322*Power(invgam13,2) -
  dADMgam311*GamDDD322*Power(invgam13,2) + (ddADMgam2222*invgam22)/2. -
  dADMgam212*GamDDD122*invgam11*invgam22 -
  dADMgam112*GamDDD222*invgam11*invgam22 -
  dADMgam222*GamDDD122*invgam12*invgam22 -
  dADMgam122*GamDDD222*invgam12*invgam22 -
  2*dADMgam212*GamDDD222*invgam12*invgam22 -
  dADMgam223*GamDDD122*invgam13*invgam22 -
  dADMgam123*GamDDD222*invgam13*invgam22 -
  dADMgam312*GamDDD222*invgam13*invgam22 -
  dADMgam212*GamDDD322*invgam13*invgam22 -
  dADMgam222*GamDDD222*Power(invgam22,2) + ddADMgam2223*invgam23 -
  dADMgam213*GamDDD122*invgam11*invgam23 -
  dADMgam312*GamDDD122*invgam11*invgam23 -
  dADMgam113*GamDDD222*invgam11*invgam23 -
  dADMgam112*GamDDD322*invgam11*invgam23 -
  dADMgam223*GamDDD122*invgam12*invgam23 -
  dADMgam322*GamDDD122*invgam12*invgam23 -
  dADMgam123*GamDDD222*invgam12*invgam23 -
  2*dADMgam213*GamDDD222*invgam12*invgam23 -
  dADMgam312*GamDDD222*invgam12*invgam23 -
  dADMgam122*GamDDD322*invgam12*invgam23 -
  dADMgam212*GamDDD322*invgam12*invgam23 -
  dADMgam233*GamDDD122*invgam13*invgam23 -
  dADMgam323*GamDDD122*invgam13*invgam23 -
  dADMgam133*GamDDD222*invgam13*invgam23 -
  dADMgam313*GamDDD222*invgam13*invgam23 -
  dADMgam123*GamDDD322*invgam13*invgam23 -
  dADMgam213*GamDDD322*invgam13*invgam23 -
  2*dADMgam312*GamDDD322*invgam13*invgam23 -
  2*dADMgam223*GamDDD222*invgam22*invgam23 -
  dADMgam322*GamDDD222*invgam22*invgam23 -
  dADMgam222*GamDDD322*invgam22*invgam23 -
  dADMgam233*GamDDD222*Power(invgam23,2) -
  dADMgam323*GamDDD222*Power(invgam23,2) -
  dADMgam223*GamDDD322*Power(invgam23,2) -
  dADMgam322*GamDDD322*Power(invgam23,2) + ddADMgam2323*invgam33 -
  (ddADMgam3322*invgam33)/2. - dADMgam313*GamDDD122*invgam11*invgam33 -
  dADMgam113*GamDDD322*invgam11*invgam33 -
  dADMgam323*GamDDD122*invgam12*invgam33 -
  dADMgam313*GamDDD222*invgam12*invgam33 -
  dADMgam123*GamDDD322*invgam12*invgam33 -
  dADMgam213*GamDDD322*invgam12*invgam33 -
  dADMgam333*GamDDD122*invgam13*invgam33 -
  dADMgam133*GamDDD322*invgam13*invgam33 -
  2*dADMgam313*GamDDD322*invgam13*invgam33 -
  dADMgam323*GamDDD222*invgam22*invgam33 -
  dADMgam223*GamDDD322*invgam22*invgam33 -
  dADMgam333*GamDDD222*invgam23*invgam33 -
  dADMgam233*GamDDD322*invgam23*invgam33 -
  2*dADMgam323*GamDDD322*invgam23*invgam33 -
  dADMgam333*GamDDD322*Power(invgam33,2)
;

vreal tr1dGam23
=
(-(ddADMgam1123*invgam11) + ddADMgam1213*invgam11 + ddADMgam1312*invgam11 -
    2*dADMgam111*GamDDD123*Power(invgam11,2) - ddADMgam1223*invgam12 +
    ddADMgam1322*invgam12 + ddADMgam2213*invgam12 + ddADMgam2312*invgam12 -
    4*dADMgam112*GamDDD123*invgam11*invgam12 -
    2*dADMgam211*GamDDD123*invgam11*invgam12 -
    2*dADMgam111*GamDDD223*invgam11*invgam12 -
    2*dADMgam122*GamDDD123*Power(invgam12,2) -
    2*dADMgam212*GamDDD123*Power(invgam12,2) -
    2*dADMgam112*GamDDD223*Power(invgam12,2) -
    2*dADMgam211*GamDDD223*Power(invgam12,2) + ddADMgam1233*invgam13 -
    ddADMgam1323*invgam13 + ddADMgam2313*invgam13 + ddADMgam3312*invgam13 -
    4*dADMgam113*GamDDD123*invgam11*invgam13 -
    2*dADMgam311*GamDDD123*invgam11*invgam13 -
    2*dADMgam111*GamDDD323*invgam11*invgam13 -
    4*dADMgam123*GamDDD123*invgam12*invgam13 -
    2*dADMgam213*GamDDD123*invgam12*invgam13 -
    2*dADMgam312*GamDDD123*invgam12*invgam13 -
    2*dADMgam113*GamDDD223*invgam12*invgam13 -
    2*dADMgam311*GamDDD223*invgam12*invgam13 -
    2*dADMgam112*GamDDD323*invgam12*invgam13 -
    2*dADMgam211*GamDDD323*invgam12*invgam13 -
    2*dADMgam133*GamDDD123*Power(invgam13,2) -
    2*dADMgam313*GamDDD123*Power(invgam13,2) -
    2*dADMgam113*GamDDD323*Power(invgam13,2) -
    2*dADMgam311*GamDDD323*Power(invgam13,2) + ddADMgam2322*invgam22 -
    2*dADMgam212*GamDDD123*invgam11*invgam22 -
    2*dADMgam112*GamDDD223*invgam11*invgam22 -
    2*dADMgam222*GamDDD123*invgam12*invgam22 -
    2*dADMgam122*GamDDD223*invgam12*invgam22 -
    4*dADMgam212*GamDDD223*invgam12*invgam22 -
    2*dADMgam223*GamDDD123*invgam13*invgam22 -
    2*dADMgam123*GamDDD223*invgam13*invgam22 -
    2*dADMgam312*GamDDD223*invgam13*invgam22 -
    2*dADMgam212*GamDDD323*invgam13*invgam22 -
    2*dADMgam222*GamDDD223*Power(invgam22,2) + ddADMgam2233*invgam23 +
    ddADMgam3322*invgam23 - 2*dADMgam213*GamDDD123*invgam11*invgam23 -
    2*dADMgam312*GamDDD123*invgam11*invgam23 -
    2*dADMgam113*GamDDD223*invgam11*invgam23 -
    2*dADMgam112*GamDDD323*invgam11*invgam23 -
    2*dADMgam223*GamDDD123*invgam12*invgam23 -
    2*dADMgam322*GamDDD123*invgam12*invgam23 -
    2*dADMgam123*GamDDD223*invgam12*invgam23 -
    4*dADMgam213*GamDDD223*invgam12*invgam23 -
    2*dADMgam312*GamDDD223*invgam12*invgam23 -
    2*dADMgam122*GamDDD323*invgam12*invgam23 -
    2*dADMgam212*GamDDD323*invgam12*invgam23 -
    2*dADMgam233*GamDDD123*invgam13*invgam23 -
    2*dADMgam323*GamDDD123*invgam13*invgam23 -
    2*dADMgam133*GamDDD223*invgam13*invgam23 -
    2*dADMgam313*GamDDD223*invgam13*invgam23 -
    2*dADMgam123*GamDDD323*invgam13*invgam23 -
    2*dADMgam213*GamDDD323*invgam13*invgam23 -
    4*dADMgam312*GamDDD323*invgam13*invgam23 -
    4*dADMgam223*GamDDD223*invgam22*invgam23 -
    2*dADMgam322*GamDDD223*invgam22*invgam23 -
    2*dADMgam222*GamDDD323*invgam22*invgam23 -
    2*dADMgam233*GamDDD223*Power(invgam23,2) -
    2*dADMgam323*GamDDD223*Power(invgam23,2) -
    2*dADMgam223*GamDDD323*Power(invgam23,2) -
    2*dADMgam322*GamDDD323*Power(invgam23,2) + ddADMgam2333*invgam33 -
    2*dADMgam313*GamDDD123*invgam11*invgam33 -
    2*dADMgam113*GamDDD323*invgam11*invgam33 -
    2*dADMgam323*GamDDD123*invgam12*invgam33 -
    2*dADMgam313*GamDDD223*invgam12*invgam33 -
    2*dADMgam123*GamDDD323*invgam12*invgam33 -
    2*dADMgam213*GamDDD323*invgam12*invgam33 -
    2*dADMgam333*GamDDD123*invgam13*invgam33 -
    2*dADMgam133*GamDDD323*invgam13*invgam33 -
    4*dADMgam313*GamDDD323*invgam13*invgam33 -
    2*dADMgam323*GamDDD223*invgam22*invgam33 -
    2*dADMgam223*GamDDD323*invgam22*invgam33 -
    2*dADMgam333*GamDDD223*invgam23*invgam33 -
    2*dADMgam233*GamDDD323*invgam23*invgam33 -
    4*dADMgam323*GamDDD323*invgam23*invgam33 -
    2*dADMgam333*GamDDD323*Power(invgam33,2))/2.
;

vreal tr1dGam33
=
-0.5*(ddADMgam1133*invgam11) + ddADMgam1313*invgam11 -
  dADMgam111*GamDDD133*Power(invgam11,2) - ddADMgam1233*invgam12 +
  ddADMgam1323*invgam12 + ddADMgam2313*invgam12 -
  2*dADMgam112*GamDDD133*invgam11*invgam12 -
  dADMgam211*GamDDD133*invgam11*invgam12 -
  dADMgam111*GamDDD233*invgam11*invgam12 -
  dADMgam122*GamDDD133*Power(invgam12,2) -
  dADMgam212*GamDDD133*Power(invgam12,2) -
  dADMgam112*GamDDD233*Power(invgam12,2) -
  dADMgam211*GamDDD233*Power(invgam12,2) + ddADMgam3313*invgam13 -
  2*dADMgam113*GamDDD133*invgam11*invgam13 -
  dADMgam311*GamDDD133*invgam11*invgam13 -
  dADMgam111*GamDDD333*invgam11*invgam13 -
  2*dADMgam123*GamDDD133*invgam12*invgam13 -
  dADMgam213*GamDDD133*invgam12*invgam13 -
  dADMgam312*GamDDD133*invgam12*invgam13 -
  dADMgam113*GamDDD233*invgam12*invgam13 -
  dADMgam311*GamDDD233*invgam12*invgam13 -
  dADMgam112*GamDDD333*invgam12*invgam13 -
  dADMgam211*GamDDD333*invgam12*invgam13 -
  dADMgam133*GamDDD133*Power(invgam13,2) -
  dADMgam313*GamDDD133*Power(invgam13,2) -
  dADMgam113*GamDDD333*Power(invgam13,2) -
  dADMgam311*GamDDD333*Power(invgam13,2) - (ddADMgam2233*invgam22)/2. +
  ddADMgam2323*invgam22 - dADMgam212*GamDDD133*invgam11*invgam22 -
  dADMgam112*GamDDD233*invgam11*invgam22 -
  dADMgam222*GamDDD133*invgam12*invgam22 -
  dADMgam122*GamDDD233*invgam12*invgam22 -
  2*dADMgam212*GamDDD233*invgam12*invgam22 -
  dADMgam223*GamDDD133*invgam13*invgam22 -
  dADMgam123*GamDDD233*invgam13*invgam22 -
  dADMgam312*GamDDD233*invgam13*invgam22 -
  dADMgam212*GamDDD333*invgam13*invgam22 -
  dADMgam222*GamDDD233*Power(invgam22,2) + ddADMgam3323*invgam23 -
  dADMgam213*GamDDD133*invgam11*invgam23 -
  dADMgam312*GamDDD133*invgam11*invgam23 -
  dADMgam113*GamDDD233*invgam11*invgam23 -
  dADMgam112*GamDDD333*invgam11*invgam23 -
  dADMgam223*GamDDD133*invgam12*invgam23 -
  dADMgam322*GamDDD133*invgam12*invgam23 -
  dADMgam123*GamDDD233*invgam12*invgam23 -
  2*dADMgam213*GamDDD233*invgam12*invgam23 -
  dADMgam312*GamDDD233*invgam12*invgam23 -
  dADMgam122*GamDDD333*invgam12*invgam23 -
  dADMgam212*GamDDD333*invgam12*invgam23 -
  dADMgam233*GamDDD133*invgam13*invgam23 -
  dADMgam323*GamDDD133*invgam13*invgam23 -
  dADMgam133*GamDDD233*invgam13*invgam23 -
  dADMgam313*GamDDD233*invgam13*invgam23 -
  dADMgam123*GamDDD333*invgam13*invgam23 -
  dADMgam213*GamDDD333*invgam13*invgam23 -
  2*dADMgam312*GamDDD333*invgam13*invgam23 -
  2*dADMgam223*GamDDD233*invgam22*invgam23 -
  dADMgam322*GamDDD233*invgam22*invgam23 -
  dADMgam222*GamDDD333*invgam22*invgam23 -
  dADMgam233*GamDDD233*Power(invgam23,2) -
  dADMgam323*GamDDD233*Power(invgam23,2) -
  dADMgam223*GamDDD333*Power(invgam23,2) -
  dADMgam322*GamDDD333*Power(invgam23,2) + (ddADMgam3333*invgam33)/2. -
  dADMgam313*GamDDD133*invgam11*invgam33 -
  dADMgam113*GamDDD333*invgam11*invgam33 -
  dADMgam323*GamDDD133*invgam12*invgam33 -
  dADMgam313*GamDDD233*invgam12*invgam33 -
  dADMgam123*GamDDD333*invgam12*invgam33 -
  dADMgam213*GamDDD333*invgam12*invgam33 -
  dADMgam333*GamDDD133*invgam13*invgam33 -
  dADMgam133*GamDDD333*invgam13*invgam33 -
  2*dADMgam313*GamDDD333*invgam13*invgam33 -
  dADMgam323*GamDDD233*invgam22*invgam33 -
  dADMgam223*GamDDD333*invgam22*invgam33 -
  dADMgam333*GamDDD233*invgam23*invgam33 -
  dADMgam233*GamDDD333*invgam23*invgam33 -
  2*dADMgam323*GamDDD333*invgam23*invgam33 -
  dADMgam333*GamDDD333*Power(invgam33,2)
;

vreal tr2dGam11
=
(ddADMgam1111*invgam11 - Power(dADMgam111,2)*Power(invgam11,2) +
    2*ddADMgam1112*invgam12 - 2*Power(dADMgam112,2)*Power(invgam12,2) +
    2*ddADMgam1113*invgam13 - 4*dADMgam112*dADMgam113*invgam12*invgam13 -
    2*Power(dADMgam113,2)*Power(invgam13,2) -
    2*dADMgam111*(2*dADMgam112*invgam11*invgam12 +
       dADMgam122*Power(invgam12,2) +
       invgam13*(2*dADMgam113*invgam11 + 2*dADMgam123*invgam12 +
          dADMgam133*invgam13)) + ddADMgam1122*invgam22 -
    2*Power(dADMgam112,2)*invgam11*invgam22 -
    4*dADMgam112*dADMgam122*invgam12*invgam22 -
    4*dADMgam112*dADMgam123*invgam13*invgam22 -
    Power(dADMgam122,2)*Power(invgam22,2) + 2*ddADMgam1123*invgam23 -
    4*dADMgam112*dADMgam113*invgam11*invgam23 -
    4*dADMgam113*dADMgam122*invgam12*invgam23 -
    4*dADMgam112*dADMgam123*invgam12*invgam23 -
    4*dADMgam113*dADMgam123*invgam13*invgam23 -
    4*dADMgam112*dADMgam133*invgam13*invgam23 -
    4*dADMgam122*dADMgam123*invgam22*invgam23 -
    2*Power(dADMgam123,2)*Power(invgam23,2) -
    2*dADMgam122*dADMgam133*Power(invgam23,2) + ddADMgam1133*invgam33 -
    2*Power(dADMgam113,2)*invgam11*invgam33 -
    4*dADMgam113*dADMgam123*invgam12*invgam33 -
    4*dADMgam113*dADMgam133*invgam13*invgam33 -
    2*Power(dADMgam123,2)*invgam22*invgam33 -
    4*dADMgam123*dADMgam133*invgam23*invgam33 -
    Power(dADMgam133,2)*Power(invgam33,2))/2.
;

vreal tr2dGam12
=
(ddADMgam1211*invgam11 + 2*ddADMgam1212*invgam12 -
    2*dADMgam112*dADMgam211*invgam11*invgam12 -
    dADMgam122*dADMgam211*Power(invgam12,2) -
    2*dADMgam112*dADMgam212*Power(invgam12,2) + 2*ddADMgam1213*invgam13 -
    2*dADMgam113*dADMgam211*invgam11*invgam13 -
    2*dADMgam123*dADMgam211*invgam12*invgam13 -
    2*dADMgam113*dADMgam212*invgam12*invgam13 -
    2*dADMgam112*dADMgam213*invgam12*invgam13 -
    dADMgam133*dADMgam211*Power(invgam13,2) -
    2*dADMgam113*dADMgam213*Power(invgam13,2) -
    dADMgam111*(dADMgam211*Power(invgam11,2) +
       2*dADMgam212*invgam11*invgam12 + dADMgam222*Power(invgam12,2) +
       2*dADMgam213*invgam11*invgam13 + 2*dADMgam223*invgam12*invgam13 +
       dADMgam233*Power(invgam13,2)) + ddADMgam1222*invgam22 -
    2*dADMgam112*dADMgam212*invgam11*invgam22 -
    2*dADMgam122*dADMgam212*invgam12*invgam22 -
    2*dADMgam112*dADMgam222*invgam12*invgam22 -
    2*dADMgam123*dADMgam212*invgam13*invgam22 -
    2*dADMgam112*dADMgam223*invgam13*invgam22 -
    dADMgam122*dADMgam222*Power(invgam22,2) + 2*ddADMgam1223*invgam23 -
    2*dADMgam113*dADMgam212*invgam11*invgam23 -
    2*dADMgam112*dADMgam213*invgam11*invgam23 -
    2*dADMgam123*dADMgam212*invgam12*invgam23 -
    2*dADMgam122*dADMgam213*invgam12*invgam23 -
    2*dADMgam113*dADMgam222*invgam12*invgam23 -
    2*dADMgam112*dADMgam223*invgam12*invgam23 -
    2*dADMgam133*dADMgam212*invgam13*invgam23 -
    2*dADMgam123*dADMgam213*invgam13*invgam23 -
    2*dADMgam113*dADMgam223*invgam13*invgam23 -
    2*dADMgam112*dADMgam233*invgam13*invgam23 -
    2*dADMgam123*dADMgam222*invgam22*invgam23 -
    2*dADMgam122*dADMgam223*invgam22*invgam23 -
    dADMgam133*dADMgam222*Power(invgam23,2) -
    2*dADMgam123*dADMgam223*Power(invgam23,2) -
    dADMgam122*dADMgam233*Power(invgam23,2) + ddADMgam1233*invgam33 -
    2*dADMgam113*dADMgam213*invgam11*invgam33 -
    2*dADMgam123*dADMgam213*invgam12*invgam33 -
    2*dADMgam113*dADMgam223*invgam12*invgam33 -
    2*dADMgam133*dADMgam213*invgam13*invgam33 -
    2*dADMgam113*dADMgam233*invgam13*invgam33 -
    2*dADMgam123*dADMgam223*invgam22*invgam33 -
    2*dADMgam133*dADMgam223*invgam23*invgam33 -
    2*dADMgam123*dADMgam233*invgam23*invgam33 -
    dADMgam133*dADMgam233*Power(invgam33,2))/2.
;

vreal tr2dGam13
=
(ddADMgam1311*invgam11 + 2*ddADMgam1312*invgam12 -
    2*dADMgam112*dADMgam311*invgam11*invgam12 -
    dADMgam122*dADMgam311*Power(invgam12,2) -
    2*dADMgam112*dADMgam312*Power(invgam12,2) + 2*ddADMgam1313*invgam13 -
    2*dADMgam113*dADMgam311*invgam11*invgam13 -
    2*dADMgam123*dADMgam311*invgam12*invgam13 -
    2*dADMgam113*dADMgam312*invgam12*invgam13 -
    2*dADMgam112*dADMgam313*invgam12*invgam13 -
    dADMgam133*dADMgam311*Power(invgam13,2) -
    2*dADMgam113*dADMgam313*Power(invgam13,2) -
    dADMgam111*(dADMgam311*Power(invgam11,2) +
       2*dADMgam312*invgam11*invgam12 + dADMgam322*Power(invgam12,2) +
       2*dADMgam313*invgam11*invgam13 + 2*dADMgam323*invgam12*invgam13 +
       dADMgam333*Power(invgam13,2)) + ddADMgam1322*invgam22 -
    2*dADMgam112*dADMgam312*invgam11*invgam22 -
    2*dADMgam122*dADMgam312*invgam12*invgam22 -
    2*dADMgam112*dADMgam322*invgam12*invgam22 -
    2*dADMgam123*dADMgam312*invgam13*invgam22 -
    2*dADMgam112*dADMgam323*invgam13*invgam22 -
    dADMgam122*dADMgam322*Power(invgam22,2) + 2*ddADMgam1323*invgam23 -
    2*dADMgam113*dADMgam312*invgam11*invgam23 -
    2*dADMgam112*dADMgam313*invgam11*invgam23 -
    2*dADMgam123*dADMgam312*invgam12*invgam23 -
    2*dADMgam122*dADMgam313*invgam12*invgam23 -
    2*dADMgam113*dADMgam322*invgam12*invgam23 -
    2*dADMgam112*dADMgam323*invgam12*invgam23 -
    2*dADMgam133*dADMgam312*invgam13*invgam23 -
    2*dADMgam123*dADMgam313*invgam13*invgam23 -
    2*dADMgam113*dADMgam323*invgam13*invgam23 -
    2*dADMgam112*dADMgam333*invgam13*invgam23 -
    2*dADMgam123*dADMgam322*invgam22*invgam23 -
    2*dADMgam122*dADMgam323*invgam22*invgam23 -
    dADMgam133*dADMgam322*Power(invgam23,2) -
    2*dADMgam123*dADMgam323*Power(invgam23,2) -
    dADMgam122*dADMgam333*Power(invgam23,2) + ddADMgam1333*invgam33 -
    2*dADMgam113*dADMgam313*invgam11*invgam33 -
    2*dADMgam123*dADMgam313*invgam12*invgam33 -
    2*dADMgam113*dADMgam323*invgam12*invgam33 -
    2*dADMgam133*dADMgam313*invgam13*invgam33 -
    2*dADMgam113*dADMgam333*invgam13*invgam33 -
    2*dADMgam123*dADMgam323*invgam22*invgam33 -
    2*dADMgam133*dADMgam323*invgam23*invgam33 -
    2*dADMgam123*dADMgam333*invgam23*invgam33 -
    dADMgam133*dADMgam333*Power(invgam33,2))/2.
;

vreal tr2dGam22
=
(ddADMgam2211*invgam11 - Power(dADMgam211,2)*Power(invgam11,2) +
    2*ddADMgam2212*invgam12 - 2*Power(dADMgam212,2)*Power(invgam12,2) +
    2*ddADMgam2213*invgam13 - 4*dADMgam212*dADMgam213*invgam12*invgam13 -
    2*Power(dADMgam213,2)*Power(invgam13,2) -
    2*dADMgam211*(2*dADMgam212*invgam11*invgam12 +
       dADMgam222*Power(invgam12,2) +
       invgam13*(2*dADMgam213*invgam11 + 2*dADMgam223*invgam12 +
          dADMgam233*invgam13)) + ddADMgam2222*invgam22 -
    2*Power(dADMgam212,2)*invgam11*invgam22 -
    4*dADMgam212*dADMgam222*invgam12*invgam22 -
    4*dADMgam212*dADMgam223*invgam13*invgam22 -
    Power(dADMgam222,2)*Power(invgam22,2) + 2*ddADMgam2223*invgam23 -
    4*dADMgam212*dADMgam213*invgam11*invgam23 -
    4*dADMgam213*dADMgam222*invgam12*invgam23 -
    4*dADMgam212*dADMgam223*invgam12*invgam23 -
    4*dADMgam213*dADMgam223*invgam13*invgam23 -
    4*dADMgam212*dADMgam233*invgam13*invgam23 -
    4*dADMgam222*dADMgam223*invgam22*invgam23 -
    2*Power(dADMgam223,2)*Power(invgam23,2) -
    2*dADMgam222*dADMgam233*Power(invgam23,2) + ddADMgam2233*invgam33 -
    2*Power(dADMgam213,2)*invgam11*invgam33 -
    4*dADMgam213*dADMgam223*invgam12*invgam33 -
    4*dADMgam213*dADMgam233*invgam13*invgam33 -
    2*Power(dADMgam223,2)*invgam22*invgam33 -
    4*dADMgam223*dADMgam233*invgam23*invgam33 -
    Power(dADMgam233,2)*Power(invgam33,2))/2.
;

vreal tr2dGam23
=
(ddADMgam2311*invgam11 + 2*ddADMgam2312*invgam12 -
    2*dADMgam212*dADMgam311*invgam11*invgam12 -
    dADMgam222*dADMgam311*Power(invgam12,2) -
    2*dADMgam212*dADMgam312*Power(invgam12,2) + 2*ddADMgam2313*invgam13 -
    2*dADMgam213*dADMgam311*invgam11*invgam13 -
    2*dADMgam223*dADMgam311*invgam12*invgam13 -
    2*dADMgam213*dADMgam312*invgam12*invgam13 -
    2*dADMgam212*dADMgam313*invgam12*invgam13 -
    dADMgam233*dADMgam311*Power(invgam13,2) -
    2*dADMgam213*dADMgam313*Power(invgam13,2) -
    dADMgam211*(dADMgam311*Power(invgam11,2) +
       2*dADMgam312*invgam11*invgam12 + dADMgam322*Power(invgam12,2) +
       2*dADMgam313*invgam11*invgam13 + 2*dADMgam323*invgam12*invgam13 +
       dADMgam333*Power(invgam13,2)) + ddADMgam2322*invgam22 -
    2*dADMgam212*dADMgam312*invgam11*invgam22 -
    2*dADMgam222*dADMgam312*invgam12*invgam22 -
    2*dADMgam212*dADMgam322*invgam12*invgam22 -
    2*dADMgam223*dADMgam312*invgam13*invgam22 -
    2*dADMgam212*dADMgam323*invgam13*invgam22 -
    dADMgam222*dADMgam322*Power(invgam22,2) + 2*ddADMgam2323*invgam23 -
    2*dADMgam213*dADMgam312*invgam11*invgam23 -
    2*dADMgam212*dADMgam313*invgam11*invgam23 -
    2*dADMgam223*dADMgam312*invgam12*invgam23 -
    2*dADMgam222*dADMgam313*invgam12*invgam23 -
    2*dADMgam213*dADMgam322*invgam12*invgam23 -
    2*dADMgam212*dADMgam323*invgam12*invgam23 -
    2*dADMgam233*dADMgam312*invgam13*invgam23 -
    2*dADMgam223*dADMgam313*invgam13*invgam23 -
    2*dADMgam213*dADMgam323*invgam13*invgam23 -
    2*dADMgam212*dADMgam333*invgam13*invgam23 -
    2*dADMgam223*dADMgam322*invgam22*invgam23 -
    2*dADMgam222*dADMgam323*invgam22*invgam23 -
    dADMgam233*dADMgam322*Power(invgam23,2) -
    2*dADMgam223*dADMgam323*Power(invgam23,2) -
    dADMgam222*dADMgam333*Power(invgam23,2) + ddADMgam2333*invgam33 -
    2*dADMgam213*dADMgam313*invgam11*invgam33 -
    2*dADMgam223*dADMgam313*invgam12*invgam33 -
    2*dADMgam213*dADMgam323*invgam12*invgam33 -
    2*dADMgam233*dADMgam313*invgam13*invgam33 -
    2*dADMgam213*dADMgam333*invgam13*invgam33 -
    2*dADMgam223*dADMgam323*invgam22*invgam33 -
    2*dADMgam233*dADMgam323*invgam23*invgam33 -
    2*dADMgam223*dADMgam333*invgam23*invgam33 -
    dADMgam233*dADMgam333*Power(invgam33,2))/2.
;

vreal tr2dGam33
=
(ddADMgam3311*invgam11 - Power(dADMgam311,2)*Power(invgam11,2) +
    2*ddADMgam3312*invgam12 - 2*Power(dADMgam312,2)*Power(invgam12,2) +
    2*ddADMgam3313*invgam13 - 4*dADMgam312*dADMgam313*invgam12*invgam13 -
    2*Power(dADMgam313,2)*Power(invgam13,2) -
    2*dADMgam311*(2*dADMgam312*invgam11*invgam12 +
       dADMgam322*Power(invgam12,2) +
       invgam13*(2*dADMgam313*invgam11 + 2*dADMgam323*invgam12 +
          dADMgam333*invgam13)) + ddADMgam3322*invgam22 -
    2*Power(dADMgam312,2)*invgam11*invgam22 -
    4*dADMgam312*dADMgam322*invgam12*invgam22 -
    4*dADMgam312*dADMgam323*invgam13*invgam22 -
    Power(dADMgam322,2)*Power(invgam22,2) + 2*ddADMgam3323*invgam23 -
    4*dADMgam312*dADMgam313*invgam11*invgam23 -
    4*dADMgam313*dADMgam322*invgam12*invgam23 -
    4*dADMgam312*dADMgam323*invgam12*invgam23 -
    4*dADMgam313*dADMgam323*invgam13*invgam23 -
    4*dADMgam312*dADMgam333*invgam13*invgam23 -
    4*dADMgam322*dADMgam323*invgam22*invgam23 -
    2*Power(dADMgam323,2)*Power(invgam23,2) -
    2*dADMgam322*dADMgam333*Power(invgam23,2) + ddADMgam3333*invgam33 -
    2*Power(dADMgam313,2)*invgam11*invgam33 -
    4*dADMgam313*dADMgam323*invgam12*invgam33 -
    4*dADMgam313*dADMgam333*invgam13*invgam33 -
    2*Power(dADMgam323,2)*invgam22*invgam33 -
    4*dADMgam323*dADMgam333*invgam23*invgam33 -
    Power(dADMgam333,2)*Power(invgam33,2))/2.
;

vreal R11
=
-(Gam112*Gam211) - Power(Gam212,2) + Gam211*Gam222 - Gam113*Gam311 +
  Gam223*Gam311 - 2*Gam213*Gam312 - Power(Gam313,2) +
  Gam111*(Gam212 + Gam313) + Gam211*Gam323 + Gam311*Gam333 + tr1dGam11 -
  tr2dGam11
;

vreal R12
=
-(Gam122*Gam211) - Gam123*Gam311 + Gam112*(Gam212 + Gam313) -
  Gam213*Gam322 + Gam212*Gam323 - Gam313*Gam323 + Gam312*Gam333 +
  tr1dGam12 - tr2dGam12
;

vreal R13
=
-(Gam123*Gam211) + Gam213*Gam222 - Gam212*Gam223 - Gam133*Gam311 -
  Gam233*Gam312 + Gam223*Gam313 + Gam113*(Gam212 + Gam313) + tr1dGam13 -
  tr2dGam13
;

vreal R22
=
-Power(Gam112,2) + Gam111*Gam122 - Gam122*Gam212 + Gam112*Gam222 -
  2*Gam123*Gam312 + Gam122*Gam313 + Gam113*Gam322 - Gam223*Gam322 +
  Gam222*Gam323 - Power(Gam323,2) + Gam322*Gam333 + tr1dGam22 - tr2dGam22
;

vreal R23
=
Gam111*Gam123 - Gam122*Gam213 + Gam112*(-Gam113 + Gam223) - Gam133*Gam312 -
  Gam233*Gam322 + Gam113*Gam323 + Gam223*Gam323 + tr1dGam23 - tr2dGam23
;

vreal R33
=
-Power(Gam113,2) + Gam111*Gam133 + Gam133*Gam212 - 2*Gam123*Gam213 -
  Power(Gam223,2) + Gam112*Gam233 + Gam222*Gam233 - Gam133*Gam313 -
  Gam233*Gam323 + Gam113*Gam333 + Gam223*Gam333 + tr1dGam33 - tr2dGam33
;

vreal trK
=
ADMK11*invgam11 + 2*ADMK12*invgam12 + 2*ADMK13*invgam13 + ADMK22*invgam22 +
  2*ADMK23*invgam23 + ADMK33*invgam33
;

vreal DADMK111
=
dADMK111 - 2*(ADMK11*Gam111 + ADMK12*Gam211 + ADMK13*Gam311)
;

vreal DADMK112
=
dADMK112 - ADMK11*Gam112 - ADMK22*Gam211 - ADMK12*(Gam111 + Gam212) -
  ADMK23*Gam311 - ADMK13*Gam312
;

vreal DADMK113
=
dADMK113 - ADMK11*Gam113 - ADMK23*Gam211 - ADMK12*Gam213 - ADMK33*Gam311 -
  ADMK13*(Gam111 + Gam313)
;

vreal DADMK122
=
dADMK122 - 2*(ADMK12*Gam112 + ADMK22*Gam212 + ADMK23*Gam312)
;

vreal DADMK123
=
dADMK123 - ADMK13*Gam112 - ADMK12*Gam113 - ADMK23*Gam212 - ADMK22*Gam213 -
  ADMK33*Gam312 - ADMK23*Gam313
;

vreal DADMK133
=
dADMK133 - 2*(ADMK13*Gam113 + ADMK23*Gam213 + ADMK33*Gam313)
;

vreal DADMK211
=
dADMK211 - 2*(ADMK11*Gam112 + ADMK12*Gam212 + ADMK13*Gam312)
;

vreal DADMK212
=
dADMK212 - ADMK11*Gam122 - ADMK22*Gam212 - ADMK12*(Gam112 + Gam222) -
  ADMK23*Gam312 - ADMK13*Gam322
;

vreal DADMK213
=
dADMK213 - ADMK11*Gam123 - ADMK23*Gam212 - ADMK12*Gam223 - ADMK33*Gam312 -
  ADMK13*(Gam112 + Gam323)
;

vreal DADMK222
=
dADMK222 - 2*(ADMK12*Gam122 + ADMK22*Gam222 + ADMK23*Gam322)
;

vreal DADMK223
=
dADMK223 - ADMK13*Gam122 - ADMK12*Gam123 - ADMK23*Gam222 - ADMK22*Gam223 -
  ADMK33*Gam322 - ADMK23*Gam323
;

vreal DADMK233
=
dADMK233 - 2*(ADMK13*Gam123 + ADMK23*Gam223 + ADMK33*Gam323)
;

vreal DADMK311
=
dADMK311 - 2*(ADMK11*Gam113 + ADMK12*Gam213 + ADMK13*Gam313)
;

vreal DADMK312
=
dADMK312 - ADMK11*Gam123 - ADMK22*Gam213 - ADMK12*(Gam113 + Gam223) -
  ADMK23*Gam313 - ADMK13*Gam323
;

vreal DADMK313
=
dADMK313 - ADMK11*Gam133 - ADMK23*Gam213 - ADMK12*Gam233 - ADMK33*Gam313 -
  ADMK13*(Gam113 + Gam333)
;

vreal DADMK322
=
dADMK322 - 2*(ADMK12*Gam123 + ADMK22*Gam223 + ADMK23*Gam323)
;

vreal DADMK323
=
dADMK323 - ADMK13*Gam123 - ADMK12*Gam133 - ADMK23*Gam223 - ADMK22*Gam233 -
  ADMK33*Gam323 - ADMK23*Gam333
;

vreal DADMK333
=
dADMK333 - 2*(ADMK13*Gam133 + ADMK23*Gam233 + ADMK33*Gam333)
;

vreal rho
=
(Power(ADMbeta1,2)*eT11 + Power(ADMbeta2,2)*eT22 +
    2*ADMbeta2*ADMbeta3*eT23 + Power(ADMbeta3,2)*eT33 +
    2*ADMbeta1*(ADMbeta2*eT12 + ADMbeta3*eT13 - eTt1) - 2*ADMbeta2*eTt2 -
    2*ADMbeta3*eTt3 + eTtt)/Power(ADMalpha,2)
;

vreal Sm1
=
(ADMbeta1*eT11 + ADMbeta2*eT12 + ADMbeta3*eT13 - eTt1)/ADMalpha
;

vreal Sm2
=
(ADMbeta1*eT12 + ADMbeta2*eT22 + ADMbeta3*eT23 - eTt2)/ADMalpha
;

vreal Sm3
=
(ADMbeta1*eT13 + ADMbeta2*eT23 + ADMbeta3*eT33 - eTt3)/ADMalpha
;


local_HC.store(mask, index2,
-(Power(ADMK11,2)*Power(invgam11,2)) - 2*Power(ADMK13,2)*Power(invgam13,2) -
  2*ADMK11*(2*ADMK12*invgam11*invgam12 + ADMK22*Power(invgam12,2) +
     invgam13*(2*ADMK13*invgam11 + 2*ADMK23*invgam12 + ADMK33*invgam13)) -
  Power(ADMK22,2)*Power(invgam22,2) -
  2*Power(ADMK12,2)*(Power(invgam12,2) + invgam11*invgam22) -
  4*ADMK13*ADMK22*invgam12*invgam23 - 4*ADMK13*ADMK23*invgam13*invgam23 -
  4*ADMK22*ADMK23*invgam22*invgam23 - 2*Power(ADMK23,2)*Power(invgam23,2) -
  2*ADMK22*ADMK33*Power(invgam23,2) -
  4*ADMK12*(ADMK13*invgam12*invgam13 + ADMK22*invgam12*invgam22 +
     ADMK23*invgam13*invgam22 + ADMK13*invgam11*invgam23 +
     ADMK23*invgam12*invgam23 + ADMK33*invgam13*invgam23) -
  2*Power(ADMK13,2)*invgam11*invgam33 - 4*ADMK13*ADMK23*invgam12*invgam33 -
  4*ADMK13*ADMK33*invgam13*invgam33 - 2*Power(ADMK23,2)*invgam22*invgam33 -
  4*ADMK23*ADMK33*invgam23*invgam33 - Power(ADMK33,2)*Power(invgam33,2) +
  invgam11*R11 + 2*invgam12*R12 + 2*invgam13*R13 + invgam22*R22 +
  2*invgam23*R23 + invgam33*R33 - 16*cpi*rho + Power(trK,2)
);

local_MC1.store(mask, index2,
-(DADMK212*Power(invgam12,2)) + 2*DADMK123*invgam12*invgam13 -
  DADMK213*invgam12*invgam13 - DADMK312*invgam12*invgam13 +
  DADMK133*Power(invgam13,2) - DADMK313*Power(invgam13,2) +
  DADMK212*invgam11*invgam22 + DADMK223*invgam13*invgam22 -
  DADMK322*invgam13*invgam22 +
  DADMK122*(Power(invgam12,2) - invgam11*invgam22) -
  2*DADMK123*invgam11*invgam23 + DADMK213*invgam11*invgam23 +
  DADMK312*invgam11*invgam23 - DADMK223*invgam12*invgam23 +
  DADMK322*invgam12*invgam23 + DADMK233*invgam13*invgam23 -
  DADMK323*invgam13*invgam23 - DADMK133*invgam11*invgam33 +
  DADMK313*invgam11*invgam33 - DADMK233*invgam12*invgam33 +
  DADMK323*invgam12*invgam33 - 8*cpi*invgam11*Sm1 - 8*cpi*invgam12*Sm2 -
  8*cpi*invgam13*Sm3
);

local_MC2.store(mask, index2,
-(DADMK112*Power(invgam12,2)) + DADMK211*Power(invgam12,2) -
  DADMK113*invgam12*invgam13 + DADMK311*invgam12*invgam13 +
  DADMK112*invgam11*invgam22 - DADMK211*invgam11*invgam22 +
  DADMK123*invgam13*invgam22 - 2*DADMK213*invgam13*invgam22 +
  DADMK312*invgam13*invgam22 + DADMK113*invgam11*invgam23 -
  DADMK311*invgam11*invgam23 - DADMK123*invgam12*invgam23 +
  2*DADMK213*invgam12*invgam23 - DADMK312*invgam12*invgam23 +
  DADMK133*invgam13*invgam23 - DADMK313*invgam13*invgam23 +
  DADMK233*Power(invgam23,2) - DADMK323*Power(invgam23,2) -
  DADMK133*invgam12*invgam33 + DADMK313*invgam12*invgam33 -
  DADMK233*invgam22*invgam33 + DADMK323*invgam22*invgam33 -
  8*cpi*invgam12*Sm1 - 8*cpi*invgam22*Sm2 - 8*cpi*invgam23*Sm3
);

local_MC3.store(mask, index2,
-(DADMK112*invgam12*invgam13) + DADMK211*invgam12*invgam13 -
  DADMK113*Power(invgam13,2) + DADMK311*Power(invgam13,2) -
  DADMK122*invgam13*invgam22 + DADMK212*invgam13*invgam22 +
  DADMK112*invgam11*invgam23 - DADMK211*invgam11*invgam23 +
  DADMK122*invgam12*invgam23 - DADMK212*invgam12*invgam23 -
  DADMK123*invgam13*invgam23 - DADMK213*invgam13*invgam23 +
  2*DADMK312*invgam13*invgam23 - DADMK223*Power(invgam23,2) +
  DADMK322*Power(invgam23,2) + DADMK113*invgam11*invgam33 -
  DADMK311*invgam11*invgam33 + DADMK123*invgam12*invgam33 +
  DADMK213*invgam12*invgam33 - 2*DADMK312*invgam12*invgam33 +
  DADMK223*invgam22*invgam33 - DADMK322*invgam22*invgam33 -
  8*cpi*invgam13*Sm1 - 8*cpi*invgam23*Sm2 - 8*cpi*invgam33*Sm3
);


  });
});

#endif // #ifndef ADM_SET_CONSTRAINT_HXX

/* ADM_set_constraint.hxx */
