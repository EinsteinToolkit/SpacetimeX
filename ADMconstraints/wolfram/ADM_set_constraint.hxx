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
const auto &tmp_gam = tl_gam(mask, index5);
const auto &tmp_exK = tl_exK(mask, index5);
const auto &tmp_alpha = tl_alpha(mask, index5);
const auto &tmp_beta = tl_beta(mask, index5);
const auto &tmp_dgam = tl_dgam(mask, index5);
const auto &tmp_dexK = tl_dexK(mask, index5);
const auto &tmp_ddgam = tl_ddgam(mask, index5);

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
const vreal gam11 = tmp_gam(0,0);
const vreal gam12 = tmp_gam(0,1);
const vreal gam13 = tmp_gam(0,2);
const vreal gam22 = tmp_gam(1,1);
const vreal gam23 = tmp_gam(1,2);
const vreal gam33 = tmp_gam(2,2);
const vreal exK11 = tmp_exK(0,0);
const vreal exK12 = tmp_exK(0,1);
const vreal exK13 = tmp_exK(0,2);
const vreal exK22 = tmp_exK(1,1);
const vreal exK23 = tmp_exK(1,2);
const vreal exK33 = tmp_exK(2,2);
const vreal alpha = tmp_alpha;
const vreal beta1 = tmp_beta(0);
const vreal beta2 = tmp_beta(1);
const vreal beta3 = tmp_beta(2);
const vreal dgam111 = tmp_dgam(0,0)(0);
const vreal dgam112 = tmp_dgam(0,1)(0);
const vreal dgam113 = tmp_dgam(0,2)(0);
const vreal dgam122 = tmp_dgam(1,1)(0);
const vreal dgam123 = tmp_dgam(1,2)(0);
const vreal dgam133 = tmp_dgam(2,2)(0);
const vreal dgam211 = tmp_dgam(0,0)(1);
const vreal dgam212 = tmp_dgam(0,1)(1);
const vreal dgam213 = tmp_dgam(0,2)(1);
const vreal dgam222 = tmp_dgam(1,1)(1);
const vreal dgam223 = tmp_dgam(1,2)(1);
const vreal dgam233 = tmp_dgam(2,2)(1);
const vreal dgam311 = tmp_dgam(0,0)(2);
const vreal dgam312 = tmp_dgam(0,1)(2);
const vreal dgam313 = tmp_dgam(0,2)(2);
const vreal dgam322 = tmp_dgam(1,1)(2);
const vreal dgam323 = tmp_dgam(1,2)(2);
const vreal dgam333 = tmp_dgam(2,2)(2);
const vreal dexK111 = tmp_dexK(0,0)(0);
const vreal dexK112 = tmp_dexK(0,1)(0);
const vreal dexK113 = tmp_dexK(0,2)(0);
const vreal dexK122 = tmp_dexK(1,1)(0);
const vreal dexK123 = tmp_dexK(1,2)(0);
const vreal dexK133 = tmp_dexK(2,2)(0);
const vreal dexK211 = tmp_dexK(0,0)(1);
const vreal dexK212 = tmp_dexK(0,1)(1);
const vreal dexK213 = tmp_dexK(0,2)(1);
const vreal dexK222 = tmp_dexK(1,1)(1);
const vreal dexK223 = tmp_dexK(1,2)(1);
const vreal dexK233 = tmp_dexK(2,2)(1);
const vreal dexK311 = tmp_dexK(0,0)(2);
const vreal dexK312 = tmp_dexK(0,1)(2);
const vreal dexK313 = tmp_dexK(0,2)(2);
const vreal dexK322 = tmp_dexK(1,1)(2);
const vreal dexK323 = tmp_dexK(1,2)(2);
const vreal dexK333 = tmp_dexK(2,2)(2);
const vreal ddgam1111 = tmp_ddgam(0,0)(0,0);
const vreal ddgam1112 = tmp_ddgam(0,1)(0,0);
const vreal ddgam1113 = tmp_ddgam(0,2)(0,0);
const vreal ddgam1122 = tmp_ddgam(1,1)(0,0);
const vreal ddgam1123 = tmp_ddgam(1,2)(0,0);
const vreal ddgam1133 = tmp_ddgam(2,2)(0,0);
const vreal ddgam1211 = tmp_ddgam(0,0)(0,1);
const vreal ddgam1212 = tmp_ddgam(0,1)(0,1);
const vreal ddgam1213 = tmp_ddgam(0,2)(0,1);
const vreal ddgam1222 = tmp_ddgam(1,1)(0,1);
const vreal ddgam1223 = tmp_ddgam(1,2)(0,1);
const vreal ddgam1233 = tmp_ddgam(2,2)(0,1);
const vreal ddgam1311 = tmp_ddgam(0,0)(0,2);
const vreal ddgam1312 = tmp_ddgam(0,1)(0,2);
const vreal ddgam1313 = tmp_ddgam(0,2)(0,2);
const vreal ddgam1322 = tmp_ddgam(1,1)(0,2);
const vreal ddgam1323 = tmp_ddgam(1,2)(0,2);
const vreal ddgam1333 = tmp_ddgam(2,2)(0,2);
const vreal ddgam2211 = tmp_ddgam(0,0)(1,1);
const vreal ddgam2212 = tmp_ddgam(0,1)(1,1);
const vreal ddgam2213 = tmp_ddgam(0,2)(1,1);
const vreal ddgam2222 = tmp_ddgam(1,1)(1,1);
const vreal ddgam2223 = tmp_ddgam(1,2)(1,1);
const vreal ddgam2233 = tmp_ddgam(2,2)(1,1);
const vreal ddgam2311 = tmp_ddgam(0,0)(1,2);
const vreal ddgam2312 = tmp_ddgam(0,1)(1,2);
const vreal ddgam2313 = tmp_ddgam(0,2)(1,2);
const vreal ddgam2322 = tmp_ddgam(1,1)(1,2);
const vreal ddgam2323 = tmp_ddgam(1,2)(1,2);
const vreal ddgam2333 = tmp_ddgam(2,2)(1,2);
const vreal ddgam3311 = tmp_ddgam(0,0)(2,2);
const vreal ddgam3312 = tmp_ddgam(0,1)(2,2);
const vreal ddgam3313 = tmp_ddgam(0,2)(2,2);
const vreal ddgam3322 = tmp_ddgam(1,1)(2,2);
const vreal ddgam3323 = tmp_ddgam(1,2)(2,2);
const vreal ddgam3333 = tmp_ddgam(2,2)(2,2);

vreal detinvgam
=
1/(-(Power(gam13,2)*gam22) + 2*gam12*gam13*gam23 - gam11*Power(gam23,2) -
    Power(gam12,2)*gam33 + gam11*gam22*gam33)
;

vreal invgam11
=
detinvgam*(-Power(gam23,2) + gam22*gam33)
;

vreal invgam12
=
detinvgam*(gam13*gam23 - gam12*gam33)
;

vreal invgam13
=
detinvgam*(-(gam13*gam22) + gam12*gam23)
;

vreal invgam22
=
detinvgam*(-Power(gam13,2) + gam11*gam33)
;

vreal invgam23
=
detinvgam*(gam12*gam13 - gam11*gam23)
;

vreal invgam33
=
detinvgam*(-Power(gam12,2) + gam11*gam22)
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

vreal tr1dGam11
=
(ddgam1111*invgam11)/2. + ddgam1112*invgam12 -
  2*dgam112*GamDDD111*invgam11*invgam12 -
  dgam211*GamDDD111*invgam11*invgam12 -
  dgam122*GamDDD111*Power(invgam12,2) -
  dgam212*GamDDD111*Power(invgam12,2) -
  dgam112*GamDDD211*Power(invgam12,2) -
  dgam211*GamDDD211*Power(invgam12,2) + ddgam1113*invgam13 -
  2*dgam113*GamDDD111*invgam11*invgam13 -
  dgam311*GamDDD111*invgam11*invgam13 -
  2*dgam123*GamDDD111*invgam12*invgam13 -
  dgam213*GamDDD111*invgam12*invgam13 -
  dgam312*GamDDD111*invgam12*invgam13 -
  dgam113*GamDDD211*invgam12*invgam13 -
  dgam311*GamDDD211*invgam12*invgam13 -
  dgam112*GamDDD311*invgam12*invgam13 -
  dgam211*GamDDD311*invgam12*invgam13 -
  dgam133*GamDDD111*Power(invgam13,2) -
  dgam313*GamDDD111*Power(invgam13,2) -
  dgam113*GamDDD311*Power(invgam13,2) -
  dgam311*GamDDD311*Power(invgam13,2) -
  dgam111*invgam11*(GamDDD111*invgam11 + GamDDD211*invgam12 +
     GamDDD311*invgam13) + ddgam1212*invgam22 - (ddgam2211*invgam22)/2. -
  dgam212*GamDDD111*invgam11*invgam22 -
  dgam112*GamDDD211*invgam11*invgam22 -
  dgam222*GamDDD111*invgam12*invgam22 -
  dgam122*GamDDD211*invgam12*invgam22 -
  2*dgam212*GamDDD211*invgam12*invgam22 -
  dgam223*GamDDD111*invgam13*invgam22 -
  dgam123*GamDDD211*invgam13*invgam22 -
  dgam312*GamDDD211*invgam13*invgam22 -
  dgam212*GamDDD311*invgam13*invgam22 -
  dgam222*GamDDD211*Power(invgam22,2) + ddgam1213*invgam23 +
  ddgam1312*invgam23 - ddgam2311*invgam23 -
  dgam213*GamDDD111*invgam11*invgam23 -
  dgam312*GamDDD111*invgam11*invgam23 -
  dgam113*GamDDD211*invgam11*invgam23 -
  dgam112*GamDDD311*invgam11*invgam23 -
  dgam223*GamDDD111*invgam12*invgam23 -
  dgam322*GamDDD111*invgam12*invgam23 -
  dgam123*GamDDD211*invgam12*invgam23 -
  2*dgam213*GamDDD211*invgam12*invgam23 -
  dgam312*GamDDD211*invgam12*invgam23 -
  dgam122*GamDDD311*invgam12*invgam23 -
  dgam212*GamDDD311*invgam12*invgam23 -
  dgam233*GamDDD111*invgam13*invgam23 -
  dgam323*GamDDD111*invgam13*invgam23 -
  dgam133*GamDDD211*invgam13*invgam23 -
  dgam313*GamDDD211*invgam13*invgam23 -
  dgam123*GamDDD311*invgam13*invgam23 -
  dgam213*GamDDD311*invgam13*invgam23 -
  2*dgam312*GamDDD311*invgam13*invgam23 -
  2*dgam223*GamDDD211*invgam22*invgam23 -
  dgam322*GamDDD211*invgam22*invgam23 -
  dgam222*GamDDD311*invgam22*invgam23 -
  dgam233*GamDDD211*Power(invgam23,2) -
  dgam323*GamDDD211*Power(invgam23,2) -
  dgam223*GamDDD311*Power(invgam23,2) -
  dgam322*GamDDD311*Power(invgam23,2) + ddgam1313*invgam33 -
  (ddgam3311*invgam33)/2. - dgam313*GamDDD111*invgam11*invgam33 -
  dgam113*GamDDD311*invgam11*invgam33 -
  dgam323*GamDDD111*invgam12*invgam33 -
  dgam313*GamDDD211*invgam12*invgam33 -
  dgam123*GamDDD311*invgam12*invgam33 -
  dgam213*GamDDD311*invgam12*invgam33 -
  dgam333*GamDDD111*invgam13*invgam33 -
  dgam133*GamDDD311*invgam13*invgam33 -
  2*dgam313*GamDDD311*invgam13*invgam33 -
  dgam323*GamDDD211*invgam22*invgam33 -
  dgam223*GamDDD311*invgam22*invgam33 -
  dgam333*GamDDD211*invgam23*invgam33 -
  dgam233*GamDDD311*invgam23*invgam33 -
  2*dgam323*GamDDD311*invgam23*invgam33 - dgam333*GamDDD311*Power(invgam33,2)
;

vreal tr1dGam12
=
(ddgam1211*invgam11 + ddgam1122*invgam12 + ddgam2211*invgam12 -
    4*dgam112*GamDDD112*invgam11*invgam12 -
    2*dgam211*GamDDD112*invgam11*invgam12 -
    2*dgam122*GamDDD112*Power(invgam12,2) -
    2*dgam212*GamDDD112*Power(invgam12,2) -
    2*dgam112*GamDDD212*Power(invgam12,2) -
    2*dgam211*GamDDD212*Power(invgam12,2) + ddgam1123*invgam13 +
    ddgam1213*invgam13 - ddgam1312*invgam13 + ddgam2311*invgam13 -
    4*dgam113*GamDDD112*invgam11*invgam13 -
    2*dgam311*GamDDD112*invgam11*invgam13 -
    4*dgam123*GamDDD112*invgam12*invgam13 -
    2*dgam213*GamDDD112*invgam12*invgam13 -
    2*dgam312*GamDDD112*invgam12*invgam13 -
    2*dgam113*GamDDD212*invgam12*invgam13 -
    2*dgam311*GamDDD212*invgam12*invgam13 -
    2*dgam112*GamDDD312*invgam12*invgam13 -
    2*dgam211*GamDDD312*invgam12*invgam13 -
    2*dgam133*GamDDD112*Power(invgam13,2) -
    2*dgam313*GamDDD112*Power(invgam13,2) -
    2*dgam113*GamDDD312*Power(invgam13,2) -
    2*dgam311*GamDDD312*Power(invgam13,2) -
    2*dgam111*invgam11*(GamDDD112*invgam11 + GamDDD212*invgam12 +
       GamDDD312*invgam13) + ddgam1222*invgam22 -
    2*dgam212*GamDDD112*invgam11*invgam22 -
    2*dgam112*GamDDD212*invgam11*invgam22 -
    2*dgam222*GamDDD112*invgam12*invgam22 -
    2*dgam122*GamDDD212*invgam12*invgam22 -
    4*dgam212*GamDDD212*invgam12*invgam22 -
    2*dgam223*GamDDD112*invgam13*invgam22 -
    2*dgam123*GamDDD212*invgam13*invgam22 -
    2*dgam312*GamDDD212*invgam13*invgam22 -
    2*dgam212*GamDDD312*invgam13*invgam22 -
    2*dgam222*GamDDD212*Power(invgam22,2) + ddgam1223*invgam23 +
    ddgam1322*invgam23 + ddgam2213*invgam23 - ddgam2312*invgam23 -
    2*dgam213*GamDDD112*invgam11*invgam23 -
    2*dgam312*GamDDD112*invgam11*invgam23 -
    2*dgam113*GamDDD212*invgam11*invgam23 -
    2*dgam112*GamDDD312*invgam11*invgam23 -
    2*dgam223*GamDDD112*invgam12*invgam23 -
    2*dgam322*GamDDD112*invgam12*invgam23 -
    2*dgam123*GamDDD212*invgam12*invgam23 -
    4*dgam213*GamDDD212*invgam12*invgam23 -
    2*dgam312*GamDDD212*invgam12*invgam23 -
    2*dgam122*GamDDD312*invgam12*invgam23 -
    2*dgam212*GamDDD312*invgam12*invgam23 -
    2*dgam233*GamDDD112*invgam13*invgam23 -
    2*dgam323*GamDDD112*invgam13*invgam23 -
    2*dgam133*GamDDD212*invgam13*invgam23 -
    2*dgam313*GamDDD212*invgam13*invgam23 -
    2*dgam123*GamDDD312*invgam13*invgam23 -
    2*dgam213*GamDDD312*invgam13*invgam23 -
    4*dgam312*GamDDD312*invgam13*invgam23 -
    4*dgam223*GamDDD212*invgam22*invgam23 -
    2*dgam322*GamDDD212*invgam22*invgam23 -
    2*dgam222*GamDDD312*invgam22*invgam23 -
    2*dgam233*GamDDD212*Power(invgam23,2) -
    2*dgam323*GamDDD212*Power(invgam23,2) -
    2*dgam223*GamDDD312*Power(invgam23,2) -
    2*dgam322*GamDDD312*Power(invgam23,2) + ddgam1323*invgam33 +
    ddgam2313*invgam33 - ddgam3312*invgam33 -
    2*dgam313*GamDDD112*invgam11*invgam33 -
    2*dgam113*GamDDD312*invgam11*invgam33 -
    2*dgam323*GamDDD112*invgam12*invgam33 -
    2*dgam313*GamDDD212*invgam12*invgam33 -
    2*dgam123*GamDDD312*invgam12*invgam33 -
    2*dgam213*GamDDD312*invgam12*invgam33 -
    2*dgam333*GamDDD112*invgam13*invgam33 -
    2*dgam133*GamDDD312*invgam13*invgam33 -
    4*dgam313*GamDDD312*invgam13*invgam33 -
    2*dgam323*GamDDD212*invgam22*invgam33 -
    2*dgam223*GamDDD312*invgam22*invgam33 -
    2*dgam333*GamDDD212*invgam23*invgam33 -
    2*dgam233*GamDDD312*invgam23*invgam33 -
    4*dgam323*GamDDD312*invgam23*invgam33 -
    2*dgam333*GamDDD312*Power(invgam33,2))/2.
;

vreal tr1dGam13
=
(ddgam1311*invgam11 + ddgam1123*invgam12 - ddgam1213*invgam12 +
    ddgam1312*invgam12 + ddgam2311*invgam12 -
    4*dgam112*GamDDD113*invgam11*invgam12 -
    2*dgam211*GamDDD113*invgam11*invgam12 -
    2*dgam122*GamDDD113*Power(invgam12,2) -
    2*dgam212*GamDDD113*Power(invgam12,2) -
    2*dgam112*GamDDD213*Power(invgam12,2) -
    2*dgam211*GamDDD213*Power(invgam12,2) + ddgam1133*invgam13 +
    ddgam3311*invgam13 - 4*dgam113*GamDDD113*invgam11*invgam13 -
    2*dgam311*GamDDD113*invgam11*invgam13 -
    4*dgam123*GamDDD113*invgam12*invgam13 -
    2*dgam213*GamDDD113*invgam12*invgam13 -
    2*dgam312*GamDDD113*invgam12*invgam13 -
    2*dgam113*GamDDD213*invgam12*invgam13 -
    2*dgam311*GamDDD213*invgam12*invgam13 -
    2*dgam112*GamDDD313*invgam12*invgam13 -
    2*dgam211*GamDDD313*invgam12*invgam13 -
    2*dgam133*GamDDD113*Power(invgam13,2) -
    2*dgam313*GamDDD113*Power(invgam13,2) -
    2*dgam113*GamDDD313*Power(invgam13,2) -
    2*dgam311*GamDDD313*Power(invgam13,2) -
    2*dgam111*invgam11*(GamDDD113*invgam11 + GamDDD213*invgam12 +
       GamDDD313*invgam13) + ddgam1223*invgam22 - ddgam2213*invgam22 +
    ddgam2312*invgam22 - 2*dgam212*GamDDD113*invgam11*invgam22 -
    2*dgam112*GamDDD213*invgam11*invgam22 -
    2*dgam222*GamDDD113*invgam12*invgam22 -
    2*dgam122*GamDDD213*invgam12*invgam22 -
    4*dgam212*GamDDD213*invgam12*invgam22 -
    2*dgam223*GamDDD113*invgam13*invgam22 -
    2*dgam123*GamDDD213*invgam13*invgam22 -
    2*dgam312*GamDDD213*invgam13*invgam22 -
    2*dgam212*GamDDD313*invgam13*invgam22 -
    2*dgam222*GamDDD213*Power(invgam22,2) + ddgam1233*invgam23 +
    ddgam1323*invgam23 - ddgam2313*invgam23 + ddgam3312*invgam23 -
    2*dgam213*GamDDD113*invgam11*invgam23 -
    2*dgam312*GamDDD113*invgam11*invgam23 -
    2*dgam113*GamDDD213*invgam11*invgam23 -
    2*dgam112*GamDDD313*invgam11*invgam23 -
    2*dgam223*GamDDD113*invgam12*invgam23 -
    2*dgam322*GamDDD113*invgam12*invgam23 -
    2*dgam123*GamDDD213*invgam12*invgam23 -
    4*dgam213*GamDDD213*invgam12*invgam23 -
    2*dgam312*GamDDD213*invgam12*invgam23 -
    2*dgam122*GamDDD313*invgam12*invgam23 -
    2*dgam212*GamDDD313*invgam12*invgam23 -
    2*dgam233*GamDDD113*invgam13*invgam23 -
    2*dgam323*GamDDD113*invgam13*invgam23 -
    2*dgam133*GamDDD213*invgam13*invgam23 -
    2*dgam313*GamDDD213*invgam13*invgam23 -
    2*dgam123*GamDDD313*invgam13*invgam23 -
    2*dgam213*GamDDD313*invgam13*invgam23 -
    4*dgam312*GamDDD313*invgam13*invgam23 -
    4*dgam223*GamDDD213*invgam22*invgam23 -
    2*dgam322*GamDDD213*invgam22*invgam23 -
    2*dgam222*GamDDD313*invgam22*invgam23 -
    2*dgam233*GamDDD213*Power(invgam23,2) -
    2*dgam323*GamDDD213*Power(invgam23,2) -
    2*dgam223*GamDDD313*Power(invgam23,2) -
    2*dgam322*GamDDD313*Power(invgam23,2) + ddgam1333*invgam33 -
    2*dgam313*GamDDD113*invgam11*invgam33 -
    2*dgam113*GamDDD313*invgam11*invgam33 -
    2*dgam323*GamDDD113*invgam12*invgam33 -
    2*dgam313*GamDDD213*invgam12*invgam33 -
    2*dgam123*GamDDD313*invgam12*invgam33 -
    2*dgam213*GamDDD313*invgam12*invgam33 -
    2*dgam333*GamDDD113*invgam13*invgam33 -
    2*dgam133*GamDDD313*invgam13*invgam33 -
    4*dgam313*GamDDD313*invgam13*invgam33 -
    2*dgam323*GamDDD213*invgam22*invgam33 -
    2*dgam223*GamDDD313*invgam22*invgam33 -
    2*dgam333*GamDDD213*invgam23*invgam33 -
    2*dgam233*GamDDD313*invgam23*invgam33 -
    4*dgam323*GamDDD313*invgam23*invgam33 -
    2*dgam333*GamDDD313*Power(invgam33,2))/2.
;

vreal tr1dGam22
=
-0.5*(ddgam1122*invgam11) + ddgam1212*invgam11 -
  dgam111*GamDDD122*Power(invgam11,2) + ddgam2212*invgam12 -
  2*dgam112*GamDDD122*invgam11*invgam12 -
  dgam211*GamDDD122*invgam11*invgam12 -
  dgam111*GamDDD222*invgam11*invgam12 -
  dgam122*GamDDD122*Power(invgam12,2) -
  dgam212*GamDDD122*Power(invgam12,2) -
  dgam112*GamDDD222*Power(invgam12,2) -
  dgam211*GamDDD222*Power(invgam12,2) + ddgam1223*invgam13 -
  ddgam1322*invgam13 + ddgam2312*invgam13 -
  2*dgam113*GamDDD122*invgam11*invgam13 -
  dgam311*GamDDD122*invgam11*invgam13 -
  dgam111*GamDDD322*invgam11*invgam13 -
  2*dgam123*GamDDD122*invgam12*invgam13 -
  dgam213*GamDDD122*invgam12*invgam13 -
  dgam312*GamDDD122*invgam12*invgam13 -
  dgam113*GamDDD222*invgam12*invgam13 -
  dgam311*GamDDD222*invgam12*invgam13 -
  dgam112*GamDDD322*invgam12*invgam13 -
  dgam211*GamDDD322*invgam12*invgam13 -
  dgam133*GamDDD122*Power(invgam13,2) -
  dgam313*GamDDD122*Power(invgam13,2) -
  dgam113*GamDDD322*Power(invgam13,2) -
  dgam311*GamDDD322*Power(invgam13,2) + (ddgam2222*invgam22)/2. -
  dgam212*GamDDD122*invgam11*invgam22 -
  dgam112*GamDDD222*invgam11*invgam22 -
  dgam222*GamDDD122*invgam12*invgam22 -
  dgam122*GamDDD222*invgam12*invgam22 -
  2*dgam212*GamDDD222*invgam12*invgam22 -
  dgam223*GamDDD122*invgam13*invgam22 -
  dgam123*GamDDD222*invgam13*invgam22 -
  dgam312*GamDDD222*invgam13*invgam22 -
  dgam212*GamDDD322*invgam13*invgam22 -
  dgam222*GamDDD222*Power(invgam22,2) + ddgam2223*invgam23 -
  dgam213*GamDDD122*invgam11*invgam23 -
  dgam312*GamDDD122*invgam11*invgam23 -
  dgam113*GamDDD222*invgam11*invgam23 -
  dgam112*GamDDD322*invgam11*invgam23 -
  dgam223*GamDDD122*invgam12*invgam23 -
  dgam322*GamDDD122*invgam12*invgam23 -
  dgam123*GamDDD222*invgam12*invgam23 -
  2*dgam213*GamDDD222*invgam12*invgam23 -
  dgam312*GamDDD222*invgam12*invgam23 -
  dgam122*GamDDD322*invgam12*invgam23 -
  dgam212*GamDDD322*invgam12*invgam23 -
  dgam233*GamDDD122*invgam13*invgam23 -
  dgam323*GamDDD122*invgam13*invgam23 -
  dgam133*GamDDD222*invgam13*invgam23 -
  dgam313*GamDDD222*invgam13*invgam23 -
  dgam123*GamDDD322*invgam13*invgam23 -
  dgam213*GamDDD322*invgam13*invgam23 -
  2*dgam312*GamDDD322*invgam13*invgam23 -
  2*dgam223*GamDDD222*invgam22*invgam23 -
  dgam322*GamDDD222*invgam22*invgam23 -
  dgam222*GamDDD322*invgam22*invgam23 -
  dgam233*GamDDD222*Power(invgam23,2) -
  dgam323*GamDDD222*Power(invgam23,2) -
  dgam223*GamDDD322*Power(invgam23,2) -
  dgam322*GamDDD322*Power(invgam23,2) + ddgam2323*invgam33 -
  (ddgam3322*invgam33)/2. - dgam313*GamDDD122*invgam11*invgam33 -
  dgam113*GamDDD322*invgam11*invgam33 -
  dgam323*GamDDD122*invgam12*invgam33 -
  dgam313*GamDDD222*invgam12*invgam33 -
  dgam123*GamDDD322*invgam12*invgam33 -
  dgam213*GamDDD322*invgam12*invgam33 -
  dgam333*GamDDD122*invgam13*invgam33 -
  dgam133*GamDDD322*invgam13*invgam33 -
  2*dgam313*GamDDD322*invgam13*invgam33 -
  dgam323*GamDDD222*invgam22*invgam33 -
  dgam223*GamDDD322*invgam22*invgam33 -
  dgam333*GamDDD222*invgam23*invgam33 -
  dgam233*GamDDD322*invgam23*invgam33 -
  2*dgam323*GamDDD322*invgam23*invgam33 - dgam333*GamDDD322*Power(invgam33,2)
;

vreal tr1dGam23
=
(-(ddgam1123*invgam11) + ddgam1213*invgam11 + ddgam1312*invgam11 -
    2*dgam111*GamDDD123*Power(invgam11,2) - ddgam1223*invgam12 +
    ddgam1322*invgam12 + ddgam2213*invgam12 + ddgam2312*invgam12 -
    4*dgam112*GamDDD123*invgam11*invgam12 -
    2*dgam211*GamDDD123*invgam11*invgam12 -
    2*dgam111*GamDDD223*invgam11*invgam12 -
    2*dgam122*GamDDD123*Power(invgam12,2) -
    2*dgam212*GamDDD123*Power(invgam12,2) -
    2*dgam112*GamDDD223*Power(invgam12,2) -
    2*dgam211*GamDDD223*Power(invgam12,2) + ddgam1233*invgam13 -
    ddgam1323*invgam13 + ddgam2313*invgam13 + ddgam3312*invgam13 -
    4*dgam113*GamDDD123*invgam11*invgam13 -
    2*dgam311*GamDDD123*invgam11*invgam13 -
    2*dgam111*GamDDD323*invgam11*invgam13 -
    4*dgam123*GamDDD123*invgam12*invgam13 -
    2*dgam213*GamDDD123*invgam12*invgam13 -
    2*dgam312*GamDDD123*invgam12*invgam13 -
    2*dgam113*GamDDD223*invgam12*invgam13 -
    2*dgam311*GamDDD223*invgam12*invgam13 -
    2*dgam112*GamDDD323*invgam12*invgam13 -
    2*dgam211*GamDDD323*invgam12*invgam13 -
    2*dgam133*GamDDD123*Power(invgam13,2) -
    2*dgam313*GamDDD123*Power(invgam13,2) -
    2*dgam113*GamDDD323*Power(invgam13,2) -
    2*dgam311*GamDDD323*Power(invgam13,2) + ddgam2322*invgam22 -
    2*dgam212*GamDDD123*invgam11*invgam22 -
    2*dgam112*GamDDD223*invgam11*invgam22 -
    2*dgam222*GamDDD123*invgam12*invgam22 -
    2*dgam122*GamDDD223*invgam12*invgam22 -
    4*dgam212*GamDDD223*invgam12*invgam22 -
    2*dgam223*GamDDD123*invgam13*invgam22 -
    2*dgam123*GamDDD223*invgam13*invgam22 -
    2*dgam312*GamDDD223*invgam13*invgam22 -
    2*dgam212*GamDDD323*invgam13*invgam22 -
    2*dgam222*GamDDD223*Power(invgam22,2) + ddgam2233*invgam23 +
    ddgam3322*invgam23 - 2*dgam213*GamDDD123*invgam11*invgam23 -
    2*dgam312*GamDDD123*invgam11*invgam23 -
    2*dgam113*GamDDD223*invgam11*invgam23 -
    2*dgam112*GamDDD323*invgam11*invgam23 -
    2*dgam223*GamDDD123*invgam12*invgam23 -
    2*dgam322*GamDDD123*invgam12*invgam23 -
    2*dgam123*GamDDD223*invgam12*invgam23 -
    4*dgam213*GamDDD223*invgam12*invgam23 -
    2*dgam312*GamDDD223*invgam12*invgam23 -
    2*dgam122*GamDDD323*invgam12*invgam23 -
    2*dgam212*GamDDD323*invgam12*invgam23 -
    2*dgam233*GamDDD123*invgam13*invgam23 -
    2*dgam323*GamDDD123*invgam13*invgam23 -
    2*dgam133*GamDDD223*invgam13*invgam23 -
    2*dgam313*GamDDD223*invgam13*invgam23 -
    2*dgam123*GamDDD323*invgam13*invgam23 -
    2*dgam213*GamDDD323*invgam13*invgam23 -
    4*dgam312*GamDDD323*invgam13*invgam23 -
    4*dgam223*GamDDD223*invgam22*invgam23 -
    2*dgam322*GamDDD223*invgam22*invgam23 -
    2*dgam222*GamDDD323*invgam22*invgam23 -
    2*dgam233*GamDDD223*Power(invgam23,2) -
    2*dgam323*GamDDD223*Power(invgam23,2) -
    2*dgam223*GamDDD323*Power(invgam23,2) -
    2*dgam322*GamDDD323*Power(invgam23,2) + ddgam2333*invgam33 -
    2*dgam313*GamDDD123*invgam11*invgam33 -
    2*dgam113*GamDDD323*invgam11*invgam33 -
    2*dgam323*GamDDD123*invgam12*invgam33 -
    2*dgam313*GamDDD223*invgam12*invgam33 -
    2*dgam123*GamDDD323*invgam12*invgam33 -
    2*dgam213*GamDDD323*invgam12*invgam33 -
    2*dgam333*GamDDD123*invgam13*invgam33 -
    2*dgam133*GamDDD323*invgam13*invgam33 -
    4*dgam313*GamDDD323*invgam13*invgam33 -
    2*dgam323*GamDDD223*invgam22*invgam33 -
    2*dgam223*GamDDD323*invgam22*invgam33 -
    2*dgam333*GamDDD223*invgam23*invgam33 -
    2*dgam233*GamDDD323*invgam23*invgam33 -
    4*dgam323*GamDDD323*invgam23*invgam33 -
    2*dgam333*GamDDD323*Power(invgam33,2))/2.
;

vreal tr1dGam33
=
-0.5*(ddgam1133*invgam11) + ddgam1313*invgam11 -
  dgam111*GamDDD133*Power(invgam11,2) - ddgam1233*invgam12 +
  ddgam1323*invgam12 + ddgam2313*invgam12 -
  2*dgam112*GamDDD133*invgam11*invgam12 -
  dgam211*GamDDD133*invgam11*invgam12 -
  dgam111*GamDDD233*invgam11*invgam12 -
  dgam122*GamDDD133*Power(invgam12,2) -
  dgam212*GamDDD133*Power(invgam12,2) -
  dgam112*GamDDD233*Power(invgam12,2) -
  dgam211*GamDDD233*Power(invgam12,2) + ddgam3313*invgam13 -
  2*dgam113*GamDDD133*invgam11*invgam13 -
  dgam311*GamDDD133*invgam11*invgam13 -
  dgam111*GamDDD333*invgam11*invgam13 -
  2*dgam123*GamDDD133*invgam12*invgam13 -
  dgam213*GamDDD133*invgam12*invgam13 -
  dgam312*GamDDD133*invgam12*invgam13 -
  dgam113*GamDDD233*invgam12*invgam13 -
  dgam311*GamDDD233*invgam12*invgam13 -
  dgam112*GamDDD333*invgam12*invgam13 -
  dgam211*GamDDD333*invgam12*invgam13 -
  dgam133*GamDDD133*Power(invgam13,2) -
  dgam313*GamDDD133*Power(invgam13,2) -
  dgam113*GamDDD333*Power(invgam13,2) -
  dgam311*GamDDD333*Power(invgam13,2) - (ddgam2233*invgam22)/2. +
  ddgam2323*invgam22 - dgam212*GamDDD133*invgam11*invgam22 -
  dgam112*GamDDD233*invgam11*invgam22 -
  dgam222*GamDDD133*invgam12*invgam22 -
  dgam122*GamDDD233*invgam12*invgam22 -
  2*dgam212*GamDDD233*invgam12*invgam22 -
  dgam223*GamDDD133*invgam13*invgam22 -
  dgam123*GamDDD233*invgam13*invgam22 -
  dgam312*GamDDD233*invgam13*invgam22 -
  dgam212*GamDDD333*invgam13*invgam22 -
  dgam222*GamDDD233*Power(invgam22,2) + ddgam3323*invgam23 -
  dgam213*GamDDD133*invgam11*invgam23 -
  dgam312*GamDDD133*invgam11*invgam23 -
  dgam113*GamDDD233*invgam11*invgam23 -
  dgam112*GamDDD333*invgam11*invgam23 -
  dgam223*GamDDD133*invgam12*invgam23 -
  dgam322*GamDDD133*invgam12*invgam23 -
  dgam123*GamDDD233*invgam12*invgam23 -
  2*dgam213*GamDDD233*invgam12*invgam23 -
  dgam312*GamDDD233*invgam12*invgam23 -
  dgam122*GamDDD333*invgam12*invgam23 -
  dgam212*GamDDD333*invgam12*invgam23 -
  dgam233*GamDDD133*invgam13*invgam23 -
  dgam323*GamDDD133*invgam13*invgam23 -
  dgam133*GamDDD233*invgam13*invgam23 -
  dgam313*GamDDD233*invgam13*invgam23 -
  dgam123*GamDDD333*invgam13*invgam23 -
  dgam213*GamDDD333*invgam13*invgam23 -
  2*dgam312*GamDDD333*invgam13*invgam23 -
  2*dgam223*GamDDD233*invgam22*invgam23 -
  dgam322*GamDDD233*invgam22*invgam23 -
  dgam222*GamDDD333*invgam22*invgam23 -
  dgam233*GamDDD233*Power(invgam23,2) -
  dgam323*GamDDD233*Power(invgam23,2) -
  dgam223*GamDDD333*Power(invgam23,2) -
  dgam322*GamDDD333*Power(invgam23,2) + (ddgam3333*invgam33)/2. -
  dgam313*GamDDD133*invgam11*invgam33 -
  dgam113*GamDDD333*invgam11*invgam33 -
  dgam323*GamDDD133*invgam12*invgam33 -
  dgam313*GamDDD233*invgam12*invgam33 -
  dgam123*GamDDD333*invgam12*invgam33 -
  dgam213*GamDDD333*invgam12*invgam33 -
  dgam333*GamDDD133*invgam13*invgam33 -
  dgam133*GamDDD333*invgam13*invgam33 -
  2*dgam313*GamDDD333*invgam13*invgam33 -
  dgam323*GamDDD233*invgam22*invgam33 -
  dgam223*GamDDD333*invgam22*invgam33 -
  dgam333*GamDDD233*invgam23*invgam33 -
  dgam233*GamDDD333*invgam23*invgam33 -
  2*dgam323*GamDDD333*invgam23*invgam33 - dgam333*GamDDD333*Power(invgam33,2)
;

vreal tr2dGam11
=
(ddgam1111*invgam11 - Power(dgam111,2)*Power(invgam11,2) +
    2*ddgam1112*invgam12 - 2*Power(dgam112,2)*Power(invgam12,2) +
    2*ddgam1113*invgam13 - 4*dgam112*dgam113*invgam12*invgam13 -
    2*Power(dgam113,2)*Power(invgam13,2) -
    2*dgam111*(2*dgam112*invgam11*invgam12 + dgam122*Power(invgam12,2) +
       invgam13*(2*dgam113*invgam11 + 2*dgam123*invgam12 +
          dgam133*invgam13)) + ddgam1122*invgam22 -
    2*Power(dgam112,2)*invgam11*invgam22 -
    4*dgam112*dgam122*invgam12*invgam22 -
    4*dgam112*dgam123*invgam13*invgam22 -
    Power(dgam122,2)*Power(invgam22,2) + 2*ddgam1123*invgam23 -
    4*dgam112*dgam113*invgam11*invgam23 -
    4*dgam113*dgam122*invgam12*invgam23 -
    4*dgam112*dgam123*invgam12*invgam23 -
    4*dgam113*dgam123*invgam13*invgam23 -
    4*dgam112*dgam133*invgam13*invgam23 -
    4*dgam122*dgam123*invgam22*invgam23 -
    2*Power(dgam123,2)*Power(invgam23,2) -
    2*dgam122*dgam133*Power(invgam23,2) + ddgam1133*invgam33 -
    2*Power(dgam113,2)*invgam11*invgam33 -
    4*dgam113*dgam123*invgam12*invgam33 -
    4*dgam113*dgam133*invgam13*invgam33 -
    2*Power(dgam123,2)*invgam22*invgam33 -
    4*dgam123*dgam133*invgam23*invgam33 - Power(dgam133,2)*Power(invgam33,2)\
)/2.
;

vreal tr2dGam12
=
(ddgam1211*invgam11 + 2*ddgam1212*invgam12 -
    2*dgam112*dgam211*invgam11*invgam12 -
    dgam122*dgam211*Power(invgam12,2) -
    2*dgam112*dgam212*Power(invgam12,2) + 2*ddgam1213*invgam13 -
    2*dgam113*dgam211*invgam11*invgam13 -
    2*dgam123*dgam211*invgam12*invgam13 -
    2*dgam113*dgam212*invgam12*invgam13 -
    2*dgam112*dgam213*invgam12*invgam13 -
    dgam133*dgam211*Power(invgam13,2) -
    2*dgam113*dgam213*Power(invgam13,2) -
    dgam111*(dgam211*Power(invgam11,2) + 2*dgam212*invgam11*invgam12 +
       dgam222*Power(invgam12,2) + 2*dgam213*invgam11*invgam13 +
       2*dgam223*invgam12*invgam13 + dgam233*Power(invgam13,2)) +
    ddgam1222*invgam22 - 2*dgam112*dgam212*invgam11*invgam22 -
    2*dgam122*dgam212*invgam12*invgam22 -
    2*dgam112*dgam222*invgam12*invgam22 -
    2*dgam123*dgam212*invgam13*invgam22 -
    2*dgam112*dgam223*invgam13*invgam22 -
    dgam122*dgam222*Power(invgam22,2) + 2*ddgam1223*invgam23 -
    2*dgam113*dgam212*invgam11*invgam23 -
    2*dgam112*dgam213*invgam11*invgam23 -
    2*dgam123*dgam212*invgam12*invgam23 -
    2*dgam122*dgam213*invgam12*invgam23 -
    2*dgam113*dgam222*invgam12*invgam23 -
    2*dgam112*dgam223*invgam12*invgam23 -
    2*dgam133*dgam212*invgam13*invgam23 -
    2*dgam123*dgam213*invgam13*invgam23 -
    2*dgam113*dgam223*invgam13*invgam23 -
    2*dgam112*dgam233*invgam13*invgam23 -
    2*dgam123*dgam222*invgam22*invgam23 -
    2*dgam122*dgam223*invgam22*invgam23 -
    dgam133*dgam222*Power(invgam23,2) -
    2*dgam123*dgam223*Power(invgam23,2) -
    dgam122*dgam233*Power(invgam23,2) + ddgam1233*invgam33 -
    2*dgam113*dgam213*invgam11*invgam33 -
    2*dgam123*dgam213*invgam12*invgam33 -
    2*dgam113*dgam223*invgam12*invgam33 -
    2*dgam133*dgam213*invgam13*invgam33 -
    2*dgam113*dgam233*invgam13*invgam33 -
    2*dgam123*dgam223*invgam22*invgam33 -
    2*dgam133*dgam223*invgam23*invgam33 -
    2*dgam123*dgam233*invgam23*invgam33 - dgam133*dgam233*Power(invgam33,2))/
  2.
;

vreal tr2dGam13
=
(ddgam1311*invgam11 + 2*ddgam1312*invgam12 -
    2*dgam112*dgam311*invgam11*invgam12 -
    dgam122*dgam311*Power(invgam12,2) -
    2*dgam112*dgam312*Power(invgam12,2) + 2*ddgam1313*invgam13 -
    2*dgam113*dgam311*invgam11*invgam13 -
    2*dgam123*dgam311*invgam12*invgam13 -
    2*dgam113*dgam312*invgam12*invgam13 -
    2*dgam112*dgam313*invgam12*invgam13 -
    dgam133*dgam311*Power(invgam13,2) -
    2*dgam113*dgam313*Power(invgam13,2) -
    dgam111*(dgam311*Power(invgam11,2) + 2*dgam312*invgam11*invgam12 +
       dgam322*Power(invgam12,2) + 2*dgam313*invgam11*invgam13 +
       2*dgam323*invgam12*invgam13 + dgam333*Power(invgam13,2)) +
    ddgam1322*invgam22 - 2*dgam112*dgam312*invgam11*invgam22 -
    2*dgam122*dgam312*invgam12*invgam22 -
    2*dgam112*dgam322*invgam12*invgam22 -
    2*dgam123*dgam312*invgam13*invgam22 -
    2*dgam112*dgam323*invgam13*invgam22 -
    dgam122*dgam322*Power(invgam22,2) + 2*ddgam1323*invgam23 -
    2*dgam113*dgam312*invgam11*invgam23 -
    2*dgam112*dgam313*invgam11*invgam23 -
    2*dgam123*dgam312*invgam12*invgam23 -
    2*dgam122*dgam313*invgam12*invgam23 -
    2*dgam113*dgam322*invgam12*invgam23 -
    2*dgam112*dgam323*invgam12*invgam23 -
    2*dgam133*dgam312*invgam13*invgam23 -
    2*dgam123*dgam313*invgam13*invgam23 -
    2*dgam113*dgam323*invgam13*invgam23 -
    2*dgam112*dgam333*invgam13*invgam23 -
    2*dgam123*dgam322*invgam22*invgam23 -
    2*dgam122*dgam323*invgam22*invgam23 -
    dgam133*dgam322*Power(invgam23,2) -
    2*dgam123*dgam323*Power(invgam23,2) -
    dgam122*dgam333*Power(invgam23,2) + ddgam1333*invgam33 -
    2*dgam113*dgam313*invgam11*invgam33 -
    2*dgam123*dgam313*invgam12*invgam33 -
    2*dgam113*dgam323*invgam12*invgam33 -
    2*dgam133*dgam313*invgam13*invgam33 -
    2*dgam113*dgam333*invgam13*invgam33 -
    2*dgam123*dgam323*invgam22*invgam33 -
    2*dgam133*dgam323*invgam23*invgam33 -
    2*dgam123*dgam333*invgam23*invgam33 - dgam133*dgam333*Power(invgam33,2))/
  2.
;

vreal tr2dGam22
=
(ddgam2211*invgam11 - Power(dgam211,2)*Power(invgam11,2) +
    2*ddgam2212*invgam12 - 2*Power(dgam212,2)*Power(invgam12,2) +
    2*ddgam2213*invgam13 - 4*dgam212*dgam213*invgam12*invgam13 -
    2*Power(dgam213,2)*Power(invgam13,2) -
    2*dgam211*(2*dgam212*invgam11*invgam12 + dgam222*Power(invgam12,2) +
       invgam13*(2*dgam213*invgam11 + 2*dgam223*invgam12 +
          dgam233*invgam13)) + ddgam2222*invgam22 -
    2*Power(dgam212,2)*invgam11*invgam22 -
    4*dgam212*dgam222*invgam12*invgam22 -
    4*dgam212*dgam223*invgam13*invgam22 -
    Power(dgam222,2)*Power(invgam22,2) + 2*ddgam2223*invgam23 -
    4*dgam212*dgam213*invgam11*invgam23 -
    4*dgam213*dgam222*invgam12*invgam23 -
    4*dgam212*dgam223*invgam12*invgam23 -
    4*dgam213*dgam223*invgam13*invgam23 -
    4*dgam212*dgam233*invgam13*invgam23 -
    4*dgam222*dgam223*invgam22*invgam23 -
    2*Power(dgam223,2)*Power(invgam23,2) -
    2*dgam222*dgam233*Power(invgam23,2) + ddgam2233*invgam33 -
    2*Power(dgam213,2)*invgam11*invgam33 -
    4*dgam213*dgam223*invgam12*invgam33 -
    4*dgam213*dgam233*invgam13*invgam33 -
    2*Power(dgam223,2)*invgam22*invgam33 -
    4*dgam223*dgam233*invgam23*invgam33 - Power(dgam233,2)*Power(invgam33,2)\
)/2.
;

vreal tr2dGam23
=
(ddgam2311*invgam11 + 2*ddgam2312*invgam12 -
    2*dgam212*dgam311*invgam11*invgam12 -
    dgam222*dgam311*Power(invgam12,2) -
    2*dgam212*dgam312*Power(invgam12,2) + 2*ddgam2313*invgam13 -
    2*dgam213*dgam311*invgam11*invgam13 -
    2*dgam223*dgam311*invgam12*invgam13 -
    2*dgam213*dgam312*invgam12*invgam13 -
    2*dgam212*dgam313*invgam12*invgam13 -
    dgam233*dgam311*Power(invgam13,2) -
    2*dgam213*dgam313*Power(invgam13,2) -
    dgam211*(dgam311*Power(invgam11,2) + 2*dgam312*invgam11*invgam12 +
       dgam322*Power(invgam12,2) + 2*dgam313*invgam11*invgam13 +
       2*dgam323*invgam12*invgam13 + dgam333*Power(invgam13,2)) +
    ddgam2322*invgam22 - 2*dgam212*dgam312*invgam11*invgam22 -
    2*dgam222*dgam312*invgam12*invgam22 -
    2*dgam212*dgam322*invgam12*invgam22 -
    2*dgam223*dgam312*invgam13*invgam22 -
    2*dgam212*dgam323*invgam13*invgam22 -
    dgam222*dgam322*Power(invgam22,2) + 2*ddgam2323*invgam23 -
    2*dgam213*dgam312*invgam11*invgam23 -
    2*dgam212*dgam313*invgam11*invgam23 -
    2*dgam223*dgam312*invgam12*invgam23 -
    2*dgam222*dgam313*invgam12*invgam23 -
    2*dgam213*dgam322*invgam12*invgam23 -
    2*dgam212*dgam323*invgam12*invgam23 -
    2*dgam233*dgam312*invgam13*invgam23 -
    2*dgam223*dgam313*invgam13*invgam23 -
    2*dgam213*dgam323*invgam13*invgam23 -
    2*dgam212*dgam333*invgam13*invgam23 -
    2*dgam223*dgam322*invgam22*invgam23 -
    2*dgam222*dgam323*invgam22*invgam23 -
    dgam233*dgam322*Power(invgam23,2) -
    2*dgam223*dgam323*Power(invgam23,2) -
    dgam222*dgam333*Power(invgam23,2) + ddgam2333*invgam33 -
    2*dgam213*dgam313*invgam11*invgam33 -
    2*dgam223*dgam313*invgam12*invgam33 -
    2*dgam213*dgam323*invgam12*invgam33 -
    2*dgam233*dgam313*invgam13*invgam33 -
    2*dgam213*dgam333*invgam13*invgam33 -
    2*dgam223*dgam323*invgam22*invgam33 -
    2*dgam233*dgam323*invgam23*invgam33 -
    2*dgam223*dgam333*invgam23*invgam33 - dgam233*dgam333*Power(invgam33,2))/
  2.
;

vreal tr2dGam33
=
(ddgam3311*invgam11 - Power(dgam311,2)*Power(invgam11,2) +
    2*ddgam3312*invgam12 - 2*Power(dgam312,2)*Power(invgam12,2) +
    2*ddgam3313*invgam13 - 4*dgam312*dgam313*invgam12*invgam13 -
    2*Power(dgam313,2)*Power(invgam13,2) -
    2*dgam311*(2*dgam312*invgam11*invgam12 + dgam322*Power(invgam12,2) +
       invgam13*(2*dgam313*invgam11 + 2*dgam323*invgam12 +
          dgam333*invgam13)) + ddgam3322*invgam22 -
    2*Power(dgam312,2)*invgam11*invgam22 -
    4*dgam312*dgam322*invgam12*invgam22 -
    4*dgam312*dgam323*invgam13*invgam22 -
    Power(dgam322,2)*Power(invgam22,2) + 2*ddgam3323*invgam23 -
    4*dgam312*dgam313*invgam11*invgam23 -
    4*dgam313*dgam322*invgam12*invgam23 -
    4*dgam312*dgam323*invgam12*invgam23 -
    4*dgam313*dgam323*invgam13*invgam23 -
    4*dgam312*dgam333*invgam13*invgam23 -
    4*dgam322*dgam323*invgam22*invgam23 -
    2*Power(dgam323,2)*Power(invgam23,2) -
    2*dgam322*dgam333*Power(invgam23,2) + ddgam3333*invgam33 -
    2*Power(dgam313,2)*invgam11*invgam33 -
    4*dgam313*dgam323*invgam12*invgam33 -
    4*dgam313*dgam333*invgam13*invgam33 -
    2*Power(dgam323,2)*invgam22*invgam33 -
    4*dgam323*dgam333*invgam23*invgam33 - Power(dgam333,2)*Power(invgam33,2)\
)/2.
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

vreal trexK
=
exK11*invgam11 + 2*exK12*invgam12 + 2*exK13*invgam13 + exK22*invgam22 +
  2*exK23*invgam23 + exK33*invgam33
;


vreal DexK112
=
dexK112 - exK11*Gam112 - exK22*Gam211 - exK12*(Gam111 + Gam212) -
  exK23*Gam311 - exK13*Gam312
;

vreal DexK113
=
dexK113 - exK11*Gam113 - exK23*Gam211 - exK12*Gam213 - exK33*Gam311 -
  exK13*(Gam111 + Gam313)
;

vreal DexK122
=
dexK122 - 2*(exK12*Gam112 + exK22*Gam212 + exK23*Gam312)
;

vreal DexK123
=
dexK123 - exK13*Gam112 - exK12*Gam113 - exK23*Gam212 - exK22*Gam213 -
  exK33*Gam312 - exK23*Gam313
;

vreal DexK133
=
dexK133 - 2*(exK13*Gam113 + exK23*Gam213 + exK33*Gam313)
;

vreal DexK211
=
dexK211 - 2*(exK11*Gam112 + exK12*Gam212 + exK13*Gam312)
;

vreal DexK212
=
dexK212 - exK11*Gam122 - exK22*Gam212 - exK12*(Gam112 + Gam222) -
  exK23*Gam312 - exK13*Gam322
;

vreal DexK213
=
dexK213 - exK11*Gam123 - exK23*Gam212 - exK12*Gam223 - exK33*Gam312 -
  exK13*(Gam112 + Gam323)
;

vreal DexK223
=
dexK223 - exK13*Gam122 - exK12*Gam123 - exK23*Gam222 - exK22*Gam223 -
  exK33*Gam322 - exK23*Gam323
;

vreal DexK233
=
dexK233 - 2*(exK13*Gam123 + exK23*Gam223 + exK33*Gam323)
;

vreal DexK311
=
dexK311 - 2*(exK11*Gam113 + exK12*Gam213 + exK13*Gam313)
;

vreal DexK312
=
dexK312 - exK11*Gam123 - exK22*Gam213 - exK12*(Gam113 + Gam223) -
  exK23*Gam313 - exK13*Gam323
;

vreal DexK313
=
dexK313 - exK11*Gam133 - exK23*Gam213 - exK12*Gam233 - exK33*Gam313 -
  exK13*(Gam113 + Gam333)
;

vreal DexK322
=
dexK322 - 2*(exK12*Gam123 + exK22*Gam223 + exK23*Gam323)
;

vreal DexK323
=
dexK323 - exK13*Gam123 - exK12*Gam133 - exK23*Gam223 - exK22*Gam233 -
  exK33*Gam323 - exK23*Gam333
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


local_HC.store(mask, index2,
-(Power(exK11,2)*Power(invgam11,2)) - 2*Power(exK13,2)*Power(invgam13,2) -
  2*exK11*(2*exK12*invgam11*invgam12 + exK22*Power(invgam12,2) +
     invgam13*(2*exK13*invgam11 + 2*exK23*invgam12 + exK33*invgam13)) -
  Power(exK22,2)*Power(invgam22,2) -
  2*Power(exK12,2)*(Power(invgam12,2) + invgam11*invgam22) -
  4*exK13*exK22*invgam12*invgam23 - 4*exK13*exK23*invgam13*invgam23 -
  4*exK22*exK23*invgam22*invgam23 - 2*Power(exK23,2)*Power(invgam23,2) -
  2*exK22*exK33*Power(invgam23,2) -
  4*exK12*(exK13*invgam12*invgam13 + exK22*invgam12*invgam22 +
     exK23*invgam13*invgam22 + exK13*invgam11*invgam23 +
     exK23*invgam12*invgam23 + exK33*invgam13*invgam23) -
  2*Power(exK13,2)*invgam11*invgam33 - 4*exK13*exK23*invgam12*invgam33 -
  4*exK13*exK33*invgam13*invgam33 - 2*Power(exK23,2)*invgam22*invgam33 -
  4*exK23*exK33*invgam23*invgam33 - Power(exK33,2)*Power(invgam33,2) +
  invgam11*R11 + 2*invgam12*R12 + 2*invgam13*R13 + invgam22*R22 +
  2*invgam23*R23 + invgam33*R33 - 16*cpi*rho + Power(trexK,2)
);

local_MC1.store(mask, index2,
-(DexK212*Power(invgam12,2)) + 2*DexK123*invgam12*invgam13 -
  DexK213*invgam12*invgam13 - DexK312*invgam12*invgam13 +
  DexK133*Power(invgam13,2) - DexK313*Power(invgam13,2) +
  DexK212*invgam11*invgam22 + DexK223*invgam13*invgam22 -
  DexK322*invgam13*invgam22 + DexK122*
   (Power(invgam12,2) - invgam11*invgam22) - 2*DexK123*invgam11*invgam23 +
  DexK213*invgam11*invgam23 + DexK312*invgam11*invgam23 -
  DexK223*invgam12*invgam23 + DexK322*invgam12*invgam23 +
  DexK233*invgam13*invgam23 - DexK323*invgam13*invgam23 -
  DexK133*invgam11*invgam33 + DexK313*invgam11*invgam33 -
  DexK233*invgam12*invgam33 + DexK323*invgam12*invgam33 -
  8*cpi*invgam11*Sm1 - 8*cpi*invgam12*Sm2 - 8*cpi*invgam13*Sm3
);

local_MC2.store(mask, index2,
-(DexK112*Power(invgam12,2)) + DexK211*Power(invgam12,2) -
  DexK113*invgam12*invgam13 + DexK311*invgam12*invgam13 +
  DexK112*invgam11*invgam22 - DexK211*invgam11*invgam22 +
  DexK123*invgam13*invgam22 - 2*DexK213*invgam13*invgam22 +
  DexK312*invgam13*invgam22 + DexK113*invgam11*invgam23 -
  DexK311*invgam11*invgam23 - DexK123*invgam12*invgam23 +
  2*DexK213*invgam12*invgam23 - DexK312*invgam12*invgam23 +
  DexK133*invgam13*invgam23 - DexK313*invgam13*invgam23 +
  DexK233*Power(invgam23,2) - DexK323*Power(invgam23,2) -
  DexK133*invgam12*invgam33 + DexK313*invgam12*invgam33 -
  DexK233*invgam22*invgam33 + DexK323*invgam22*invgam33 -
  8*cpi*invgam12*Sm1 - 8*cpi*invgam22*Sm2 - 8*cpi*invgam23*Sm3
);

local_MC3.store(mask, index2,
-(DexK112*invgam12*invgam13) + DexK211*invgam12*invgam13 -
  DexK113*Power(invgam13,2) + DexK311*Power(invgam13,2) -
  DexK122*invgam13*invgam22 + DexK212*invgam13*invgam22 +
  DexK112*invgam11*invgam23 - DexK211*invgam11*invgam23 +
  DexK122*invgam12*invgam23 - DexK212*invgam12*invgam23 -
  DexK123*invgam13*invgam23 - DexK213*invgam13*invgam23 +
  2*DexK312*invgam13*invgam23 - DexK223*Power(invgam23,2) +
  DexK322*Power(invgam23,2) + DexK113*invgam11*invgam33 -
  DexK311*invgam11*invgam33 + DexK123*invgam12*invgam33 +
  DexK213*invgam12*invgam33 - 2*DexK312*invgam12*invgam33 +
  DexK223*invgam22*invgam33 - DexK322*invgam22*invgam33 -
  8*cpi*invgam13*Sm1 - 8*cpi*invgam23*Sm2 - 8*cpi*invgam33*Sm3
);


  });
});

#endif // #ifndef ADM_SET_CONSTRAINT_HXX

/* ADM_set_constraint.hxx */
