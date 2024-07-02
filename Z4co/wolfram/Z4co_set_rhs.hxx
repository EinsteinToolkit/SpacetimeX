/* Z4co_set_rhs.hxx */
/* Produced with Mathematica */

vreal &dtchi = gf_dtchi(mask, index5);
vreal &dtgamt11 = gf_dtgamt(mask, index5)(0,0);
vreal &dtgamt12 = gf_dtgamt(mask, index5)(0,1);
vreal &dtgamt13 = gf_dtgamt(mask, index5)(0,2);
vreal &dtgamt22 = gf_dtgamt(mask, index5)(1,1);
vreal &dtgamt23 = gf_dtgamt(mask, index5)(1,2);
vreal &dtgamt33 = gf_dtgamt(mask, index5)(2,2);
vreal &dtexKh = gf_dtexKh(mask, index5);
vreal &dtexAt11 = gf_dtexAt(mask, index5)(0,0);
vreal &dtexAt12 = gf_dtexAt(mask, index5)(0,1);
vreal &dtexAt13 = gf_dtexAt(mask, index5)(0,2);
vreal &dtexAt22 = gf_dtexAt(mask, index5)(1,1);
vreal &dtexAt23 = gf_dtexAt(mask, index5)(1,2);
vreal &dtexAt33 = gf_dtexAt(mask, index5)(2,2);
vreal &dttrGt1 = gf_dttrGt(mask, index5)(0);
vreal &dttrGt2 = gf_dttrGt(mask, index5)(1);
vreal &dttrGt3 = gf_dttrGt(mask, index5)(2);
vreal &dtTheta = gf_dtTheta(mask, index5);
vreal &dtalpha = gf_dtalpha(mask, index5);
vreal &dtbeta1 = gf_dtbeta(mask, index5)(0);
vreal &dtbeta2 = gf_dtbeta(mask, index5)(1);
vreal &dtbeta3 = gf_dtbeta(mask, index5)(2);
vreal &ZtC1 = gf_ZtC(mask, index5)(0);
vreal &ZtC2 = gf_ZtC(mask, index5)(1);
vreal &ZtC3 = gf_ZtC(mask, index5)(2);
vreal &HC = gf_HC(mask, index5);
vreal &MtC1 = gf_MtC(mask, index5)(0);
vreal &MtC2 = gf_MtC(mask, index5)(1);
vreal &MtC3 = gf_MtC(mask, index5)(2);
vreal &eTtt = gf_eTtt(mask, index5);
vreal &eTt1 = gf_eTt(mask, index5)(0);
vreal &eTt2 = gf_eTt(mask, index5)(1);
vreal &eTt3 = gf_eTt(mask, index5)(2);
vreal &eT11 = gf_eT(mask, index5)(0,0);
vreal &eT12 = gf_eT(mask, index5)(0,1);
vreal &eT13 = gf_eT(mask, index5)(0,2);
vreal &eT21 = gf_eT(mask, index5)(1,0);
vreal &eT22 = gf_eT(mask, index5)(1,1);
vreal &eT23 = gf_eT(mask, index5)(1,2);
vreal &eT31 = gf_eT(mask, index5)(2,0);
vreal &eT32 = gf_eT(mask, index5)(2,1);
vreal &eT33 = gf_eT(mask, index5)(2,2);
vreal &chi = gf_chi(mask, index2);
vreal &gamt11 = gf_gamt(mask, index2)(0,0);
vreal &gamt12 = gf_gamt(mask, index2)(0,1);
vreal &gamt13 = gf_gamt(mask, index2)(0,2);
vreal &gamt22 = gf_gamt(mask, index2)(1,1);
vreal &gamt23 = gf_gamt(mask, index2)(1,2);
vreal &gamt33 = gf_gamt(mask, index2)(2,2);
vreal &exKh = gf_exKh(mask, index2);
vreal &exAt11 = gf_exAt(mask, index2)(0,0);
vreal &exAt12 = gf_exAt(mask, index2)(0,1);
vreal &exAt13 = gf_exAt(mask, index2)(0,2);
vreal &exAt22 = gf_exAt(mask, index2)(1,1);
vreal &exAt23 = gf_exAt(mask, index2)(1,2);
vreal &exAt33 = gf_exAt(mask, index2)(2,2);
vreal &trGt1 = gf_trGt(mask, index2)(0);
vreal &trGt2 = gf_trGt(mask, index2)(1);
vreal &trGt3 = gf_trGt(mask, index2)(2);
vreal &Theta = gf_Theta(mask, index2);
vreal &alpha = gf_alpha(mask, index2);
vreal &beta1 = gf_beta(mask, index2)(0);
vreal &beta2 = gf_beta(mask, index2)(1);
vreal &beta3 = gf_beta(mask, index2)(2);
vreal &dchi1 = gf_dchi(mask, index2)(0);
vreal &dchi2 = gf_dchi(mask, index2)(1);
vreal &dchi3 = gf_dchi(mask, index2)(2);
vreal &dgamt111 = gf_dgamt(mask, index2)(0,0)(0);
vreal &dgamt112 = gf_dgamt(mask, index2)(0,1)(0);
vreal &dgamt113 = gf_dgamt(mask, index2)(0,2)(0);
vreal &dgamt122 = gf_dgamt(mask, index2)(1,1)(0);
vreal &dgamt123 = gf_dgamt(mask, index2)(1,2)(0);
vreal &dgamt133 = gf_dgamt(mask, index2)(2,2)(0);
vreal &dgamt211 = gf_dgamt(mask, index2)(0,0)(1);
vreal &dgamt212 = gf_dgamt(mask, index2)(0,1)(1);
vreal &dgamt213 = gf_dgamt(mask, index2)(0,2)(1);
vreal &dgamt222 = gf_dgamt(mask, index2)(1,1)(1);
vreal &dgamt223 = gf_dgamt(mask, index2)(1,2)(1);
vreal &dgamt233 = gf_dgamt(mask, index2)(2,2)(1);
vreal &dgamt311 = gf_dgamt(mask, index2)(0,0)(2);
vreal &dgamt312 = gf_dgamt(mask, index2)(0,1)(2);
vreal &dgamt313 = gf_dgamt(mask, index2)(0,2)(2);
vreal &dgamt322 = gf_dgamt(mask, index2)(1,1)(2);
vreal &dgamt323 = gf_dgamt(mask, index2)(1,2)(2);
vreal &dgamt333 = gf_dgamt(mask, index2)(2,2)(2);
vreal &dexKh1 = gf_dexKh(mask, index2)(0);
vreal &dexKh2 = gf_dexKh(mask, index2)(1);
vreal &dexKh3 = gf_dexKh(mask, index2)(2);
vreal &dexAt111 = gf_dexAt(mask, index2)(0,0)(0);
vreal &dexAt112 = gf_dexAt(mask, index2)(0,1)(0);
vreal &dexAt113 = gf_dexAt(mask, index2)(0,2)(0);
vreal &dexAt122 = gf_dexAt(mask, index2)(1,1)(0);
vreal &dexAt123 = gf_dexAt(mask, index2)(1,2)(0);
vreal &dexAt133 = gf_dexAt(mask, index2)(2,2)(0);
vreal &dexAt211 = gf_dexAt(mask, index2)(0,0)(1);
vreal &dexAt212 = gf_dexAt(mask, index2)(0,1)(1);
vreal &dexAt213 = gf_dexAt(mask, index2)(0,2)(1);
vreal &dexAt222 = gf_dexAt(mask, index2)(1,1)(1);
vreal &dexAt223 = gf_dexAt(mask, index2)(1,2)(1);
vreal &dexAt233 = gf_dexAt(mask, index2)(2,2)(1);
vreal &dexAt311 = gf_dexAt(mask, index2)(0,0)(2);
vreal &dexAt312 = gf_dexAt(mask, index2)(0,1)(2);
vreal &dexAt313 = gf_dexAt(mask, index2)(0,2)(2);
vreal &dexAt322 = gf_dexAt(mask, index2)(1,1)(2);
vreal &dexAt323 = gf_dexAt(mask, index2)(1,2)(2);
vreal &dexAt333 = gf_dexAt(mask, index2)(2,2)(2);
vreal &dtrGt11 = gf_dtrGt(mask, index2)(0)(0);
vreal &dtrGt12 = gf_dtrGt(mask, index2)(1)(0);
vreal &dtrGt13 = gf_dtrGt(mask, index2)(2)(0);
vreal &dtrGt21 = gf_dtrGt(mask, index2)(0)(1);
vreal &dtrGt22 = gf_dtrGt(mask, index2)(1)(1);
vreal &dtrGt23 = gf_dtrGt(mask, index2)(2)(1);
vreal &dtrGt31 = gf_dtrGt(mask, index2)(0)(2);
vreal &dtrGt32 = gf_dtrGt(mask, index2)(1)(2);
vreal &dtrGt33 = gf_dtrGt(mask, index2)(2)(2);
vreal &dTheta1 = gf_dTheta(mask, index2)(0);
vreal &dTheta2 = gf_dTheta(mask, index2)(1);
vreal &dTheta3 = gf_dTheta(mask, index2)(2);
vreal &dalpha1 = gf_dalpha(mask, index2)(0);
vreal &dalpha2 = gf_dalpha(mask, index2)(1);
vreal &dalpha3 = gf_dalpha(mask, index2)(2);
vreal &dbeta11 = gf_dbeta(mask, index2)(0)(0);
vreal &dbeta12 = gf_dbeta(mask, index2)(1)(0);
vreal &dbeta13 = gf_dbeta(mask, index2)(2)(0);
vreal &dbeta21 = gf_dbeta(mask, index2)(0)(1);
vreal &dbeta22 = gf_dbeta(mask, index2)(1)(1);
vreal &dbeta23 = gf_dbeta(mask, index2)(2)(1);
vreal &dbeta31 = gf_dbeta(mask, index2)(0)(2);
vreal &dbeta32 = gf_dbeta(mask, index2)(1)(2);
vreal &dbeta33 = gf_dbeta(mask, index2)(2)(2);
vreal &ddchi11 = gf_ddchi(mask, index2)(0,0);
vreal &ddchi12 = gf_ddchi(mask, index2)(0,1);
vreal &ddchi13 = gf_ddchi(mask, index2)(0,2);
vreal &ddchi22 = gf_ddchi(mask, index2)(1,1);
vreal &ddchi23 = gf_ddchi(mask, index2)(1,2);
vreal &ddchi33 = gf_ddchi(mask, index2)(2,2);
vreal &ddgamt1111 = gf_ddgamt(mask, index2)(0,0)(0,0);
vreal &ddgamt1112 = gf_ddgamt(mask, index2)(0,1)(0,0);
vreal &ddgamt1113 = gf_ddgamt(mask, index2)(0,2)(0,0);
vreal &ddgamt1122 = gf_ddgamt(mask, index2)(1,1)(0,0);
vreal &ddgamt1123 = gf_ddgamt(mask, index2)(1,2)(0,0);
vreal &ddgamt1133 = gf_ddgamt(mask, index2)(2,2)(0,0);
vreal &ddgamt1211 = gf_ddgamt(mask, index2)(0,0)(0,1);
vreal &ddgamt1212 = gf_ddgamt(mask, index2)(0,1)(0,1);
vreal &ddgamt1213 = gf_ddgamt(mask, index2)(0,2)(0,1);
vreal &ddgamt1222 = gf_ddgamt(mask, index2)(1,1)(0,1);
vreal &ddgamt1223 = gf_ddgamt(mask, index2)(1,2)(0,1);
vreal &ddgamt1233 = gf_ddgamt(mask, index2)(2,2)(0,1);
vreal &ddgamt1311 = gf_ddgamt(mask, index2)(0,0)(0,2);
vreal &ddgamt1312 = gf_ddgamt(mask, index2)(0,1)(0,2);
vreal &ddgamt1313 = gf_ddgamt(mask, index2)(0,2)(0,2);
vreal &ddgamt1322 = gf_ddgamt(mask, index2)(1,1)(0,2);
vreal &ddgamt1323 = gf_ddgamt(mask, index2)(1,2)(0,2);
vreal &ddgamt1333 = gf_ddgamt(mask, index2)(2,2)(0,2);
vreal &ddgamt2211 = gf_ddgamt(mask, index2)(0,0)(1,1);
vreal &ddgamt2212 = gf_ddgamt(mask, index2)(0,1)(1,1);
vreal &ddgamt2213 = gf_ddgamt(mask, index2)(0,2)(1,1);
vreal &ddgamt2222 = gf_ddgamt(mask, index2)(1,1)(1,1);
vreal &ddgamt2223 = gf_ddgamt(mask, index2)(1,2)(1,1);
vreal &ddgamt2233 = gf_ddgamt(mask, index2)(2,2)(1,1);
vreal &ddgamt2311 = gf_ddgamt(mask, index2)(0,0)(1,2);
vreal &ddgamt2312 = gf_ddgamt(mask, index2)(0,1)(1,2);
vreal &ddgamt2313 = gf_ddgamt(mask, index2)(0,2)(1,2);
vreal &ddgamt2322 = gf_ddgamt(mask, index2)(1,1)(1,2);
vreal &ddgamt2323 = gf_ddgamt(mask, index2)(1,2)(1,2);
vreal &ddgamt2333 = gf_ddgamt(mask, index2)(2,2)(1,2);
vreal &ddgamt3311 = gf_ddgamt(mask, index2)(0,0)(2,2);
vreal &ddgamt3312 = gf_ddgamt(mask, index2)(0,1)(2,2);
vreal &ddgamt3313 = gf_ddgamt(mask, index2)(0,2)(2,2);
vreal &ddgamt3322 = gf_ddgamt(mask, index2)(1,1)(2,2);
vreal &ddgamt3323 = gf_ddgamt(mask, index2)(1,2)(2,2);
vreal &ddgamt3333 = gf_ddgamt(mask, index2)(2,2)(2,2);
vreal &ddalpha11 = gf_ddalpha(mask, index2)(0,0);
vreal &ddalpha12 = gf_ddalpha(mask, index2)(0,1);
vreal &ddalpha13 = gf_ddalpha(mask, index2)(0,2);
vreal &ddalpha22 = gf_ddalpha(mask, index2)(1,1);
vreal &ddalpha23 = gf_ddalpha(mask, index2)(1,2);
vreal &ddalpha33 = gf_ddalpha(mask, index2)(2,2);
vreal &ddbeta111 = gf_ddbeta(mask, index2)(0)(0,0);
vreal &ddbeta112 = gf_ddbeta(mask, index2)(1)(0,0);
vreal &ddbeta113 = gf_ddbeta(mask, index2)(2)(0,0);
vreal &ddbeta121 = gf_ddbeta(mask, index2)(0)(0,1);
vreal &ddbeta122 = gf_ddbeta(mask, index2)(1)(0,1);
vreal &ddbeta123 = gf_ddbeta(mask, index2)(2)(0,1);
vreal &ddbeta131 = gf_ddbeta(mask, index2)(0)(0,2);
vreal &ddbeta132 = gf_ddbeta(mask, index2)(1)(0,2);
vreal &ddbeta133 = gf_ddbeta(mask, index2)(2)(0,2);
vreal &ddbeta221 = gf_ddbeta(mask, index2)(0)(1,1);
vreal &ddbeta222 = gf_ddbeta(mask, index2)(1)(1,1);
vreal &ddbeta223 = gf_ddbeta(mask, index2)(2)(1,1);
vreal &ddbeta231 = gf_ddbeta(mask, index2)(0)(1,2);
vreal &ddbeta232 = gf_ddbeta(mask, index2)(1)(1,2);
vreal &ddbeta233 = gf_ddbeta(mask, index2)(2)(1,2);
vreal &ddbeta331 = gf_ddbeta(mask, index2)(0)(2,2);
vreal &ddbeta332 = gf_ddbeta(mask, index2)(1)(2,2);
vreal &ddbeta333 = gf_ddbeta(mask, index2)(2)(2,2);

vreal detinvgamt
=
1/(-(Power(gamt13,2)*gamt22) + 2*gamt12*gamt13*gamt23 - 
    gamt11*Power(gamt23,2) - Power(gamt12,2)*gamt33 + gamt11*gamt22*gamt33)
;

vreal invgamt11
=
detinvgamt*(-Power(gamt23,2) + gamt22*gamt33)
;

vreal invgamt12
=
detinvgamt*(gamt13*gamt23 - gamt12*gamt33)
;

vreal invgamt13
=
detinvgamt*(-(gamt13*gamt22) + gamt12*gamt23)
;

vreal invgamt22
=
detinvgamt*(-Power(gamt13,2) + gamt11*gamt33)
;

vreal invgamt23
=
detinvgamt*(gamt12*gamt13 - gamt11*gamt23)
;

vreal invgamt33
=
detinvgamt*(-Power(gamt12,2) + gamt11*gamt22)
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
-0.5*(2*chi*(-ddalpha11 + dalpha2*Gt211 + dalpha3*Gt311) + 
     dalpha1*(2*chi*Gt111 + dchi1*(-2 + gamt11*invgamt11) + 
        dchi2*gamt11*invgamt12 + dchi3*gamt11*invgamt13) + 
     gamt11*(dalpha2*dchi1*invgamt12 + dalpha3*dchi1*invgamt13 + 
        dalpha2*dchi2*invgamt22 + dalpha3*dchi2*invgamt23 + 
        dalpha2*dchi3*invgamt23 + dalpha3*dchi3*invgamt33))/chi
;

vreal DDalpha12
=
-0.5*(-2*chi*ddalpha12 + 2*chi*dalpha3*Gt312 + 
     dalpha3*dchi1*gamt12*invgamt13 + 
     dalpha1*(-dchi2 + 2*chi*Gt112 + dchi1*gamt12*invgamt11 + 
        dchi2*gamt12*invgamt12 + dchi3*gamt12*invgamt13) + 
     dalpha3*dchi2*gamt12*invgamt23 + 
     dalpha2*(2*chi*Gt212 + dchi1*(-1 + gamt12*invgamt12) + 
        dchi2*gamt12*invgamt22 + dchi3*gamt12*invgamt23) + 
     dalpha3*dchi3*gamt12*invgamt33)/chi
;

vreal DDalpha13
=
-0.5*(-2*chi*ddalpha13 + 2*chi*dalpha2*Gt213 + 
     dalpha2*dchi1*gamt13*invgamt12 + 
     dalpha1*(-dchi3 + 2*chi*Gt113 + dchi1*gamt13*invgamt11 + 
        dchi2*gamt13*invgamt12 + dchi3*gamt13*invgamt13) + 
     dalpha2*dchi2*gamt13*invgamt22 + dalpha2*dchi3*gamt13*invgamt23 + 
     dalpha3*(2*chi*Gt313 + dchi1*(-1 + gamt13*invgamt13) + 
        dchi2*gamt13*invgamt23 + dchi3*gamt13*invgamt33))/chi
;

vreal DDalpha22
=
-0.5*(2*chi*(-ddalpha22 + dalpha1*Gt122 + dalpha3*Gt322) + 
     dalpha2*(2*chi*Gt222 + dchi1*gamt22*invgamt12 + 
        dchi2*(-2 + gamt22*invgamt22) + dchi3*gamt22*invgamt23) + 
     gamt22*(dalpha1*dchi1*invgamt11 + dalpha1*dchi2*invgamt12 + 
        dalpha3*dchi1*invgamt13 + dalpha1*dchi3*invgamt13 + 
        dalpha3*dchi2*invgamt23 + dalpha3*dchi3*invgamt33))/chi
;

vreal DDalpha23
=
-0.5*(-2*chi*ddalpha23 + 2*chi*dalpha1*Gt123 + 
     dalpha1*dchi1*gamt23*invgamt11 + dalpha1*dchi2*gamt23*invgamt12 + 
     dalpha1*dchi3*gamt23*invgamt13 + 
     dalpha2*(-dchi3 + 2*chi*Gt223 + dchi1*gamt23*invgamt12 + 
        dchi2*gamt23*invgamt22 + dchi3*gamt23*invgamt23) + 
     dalpha3*(2*chi*Gt323 + dchi1*gamt23*invgamt13 + 
        dchi2*(-1 + gamt23*invgamt23) + dchi3*gamt23*invgamt33))/chi
;

vreal DDalpha33
=
-0.5*(2*chi*(-ddalpha33 + dalpha1*Gt133 + dalpha2*Gt233) + 
     gamt33*(dalpha1*dchi1*invgamt11 + dalpha2*dchi1*invgamt12 + 
        dalpha1*dchi2*invgamt12 + dalpha1*dchi3*invgamt13 + 
        dalpha2*dchi2*invgamt22 + dalpha2*dchi3*invgamt23) + 
     dalpha3*(2*chi*Gt333 + dchi1*gamt33*invgamt13 + 
        dchi2*gamt33*invgamt23 + dchi3*(-2 + gamt33*invgamt33)))/chi
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

vreal dexAtUU111
=
-2*dgamt111*exAtUU11*invgamt11 - 2*dgamt113*exAtUU13*invgamt11 + 
  dexAt111*Power(invgamt11,2) - 2*dgamt122*exAtUU12*invgamt12 - 
  2*dgamt123*exAtUU13*invgamt12 + 2*dexAt112*invgamt11*invgamt12 + 
  dexAt122*Power(invgamt12,2) - 
  2*dgamt112*(exAtUU12*invgamt11 + exAtUU11*invgamt12) - 
  2*dgamt113*exAtUU11*invgamt13 - 2*dgamt123*exAtUU12*invgamt13 - 
  2*dgamt133*exAtUU13*invgamt13 + 2*dexAt113*invgamt11*invgamt13 + 
  2*dexAt123*invgamt12*invgamt13 + dexAt133*Power(invgamt13,2)
;

vreal dexAtUU112
=
-2*dgamt111*exAtUU12*invgamt11 - 2*dgamt113*exAtUU23*invgamt11 - 
  2*dgamt122*exAtUU22*invgamt12 - 2*dgamt123*exAtUU23*invgamt12 + 
  dexAt111*invgamt11*invgamt12 + dexAt112*Power(invgamt12,2) - 
  2*dgamt112*(exAtUU22*invgamt11 + exAtUU12*invgamt12) - 
  2*dgamt113*exAtUU12*invgamt13 - 2*dgamt123*exAtUU22*invgamt13 - 
  2*dgamt133*exAtUU23*invgamt13 + dexAt113*invgamt12*invgamt13 + 
  dexAt112*invgamt11*invgamt22 + dexAt122*invgamt12*invgamt22 + 
  dexAt123*invgamt13*invgamt22 + dexAt113*invgamt11*invgamt23 + 
  dexAt123*invgamt12*invgamt23 + dexAt133*invgamt13*invgamt23
;

vreal dexAtUU113
=
-2*dgamt111*exAtUU13*invgamt11 - 2*dgamt113*exAtUU33*invgamt11 - 
  2*dgamt122*exAtUU23*invgamt12 - 2*dgamt123*exAtUU33*invgamt12 - 
  2*dgamt112*(exAtUU23*invgamt11 + exAtUU13*invgamt12) - 
  2*dgamt113*exAtUU13*invgamt13 - 2*dgamt123*exAtUU23*invgamt13 - 
  2*dgamt133*exAtUU33*invgamt13 + dexAt111*invgamt11*invgamt13 + 
  dexAt112*invgamt12*invgamt13 + dexAt113*Power(invgamt13,2) + 
  dexAt112*invgamt11*invgamt23 + dexAt122*invgamt12*invgamt23 + 
  dexAt123*invgamt13*invgamt23 + dexAt113*invgamt11*invgamt33 + 
  dexAt123*invgamt12*invgamt33 + dexAt133*invgamt13*invgamt33
;

vreal dexAtUU122
=
-2*dgamt111*exAtUU12*invgamt12 - 2*dgamt113*exAtUU23*invgamt12 + 
  dexAt111*Power(invgamt12,2) - 2*dgamt122*exAtUU22*invgamt22 - 
  2*dgamt123*exAtUU23*invgamt22 + 2*dexAt112*invgamt12*invgamt22 + 
  dexAt122*Power(invgamt22,2) - 
  2*dgamt112*(exAtUU22*invgamt12 + exAtUU12*invgamt22) - 
  2*dgamt113*exAtUU12*invgamt23 - 2*dgamt123*exAtUU22*invgamt23 - 
  2*dgamt133*exAtUU23*invgamt23 + 2*dexAt113*invgamt12*invgamt23 + 
  2*dexAt123*invgamt22*invgamt23 + dexAt133*Power(invgamt23,2)
;

vreal dexAtUU123
=
-2*dgamt111*exAtUU13*invgamt12 - 2*dgamt113*exAtUU33*invgamt12 + 
  dexAt111*invgamt12*invgamt13 - 2*dgamt122*exAtUU23*invgamt22 - 
  2*dgamt123*exAtUU33*invgamt22 + dexAt112*invgamt13*invgamt22 - 
  2*dgamt112*(exAtUU23*invgamt12 + exAtUU13*invgamt22) - 
  2*dgamt113*exAtUU13*invgamt23 - 2*dgamt123*exAtUU23*invgamt23 - 
  2*dgamt133*exAtUU33*invgamt23 + dexAt112*invgamt12*invgamt23 + 
  dexAt113*invgamt13*invgamt23 + dexAt122*invgamt22*invgamt23 + 
  dexAt123*Power(invgamt23,2) + dexAt113*invgamt12*invgamt33 + 
  dexAt123*invgamt22*invgamt33 + dexAt133*invgamt23*invgamt33
;

vreal dexAtUU133
=
-2*dgamt111*exAtUU13*invgamt13 - 2*dgamt113*exAtUU33*invgamt13 + 
  dexAt111*Power(invgamt13,2) - 2*dgamt122*exAtUU23*invgamt23 - 
  2*dgamt123*exAtUU33*invgamt23 + 2*dexAt112*invgamt13*invgamt23 + 
  dexAt122*Power(invgamt23,2) - 
  2*dgamt112*(exAtUU23*invgamt13 + exAtUU13*invgamt23) - 
  2*dgamt113*exAtUU13*invgamt33 - 2*dgamt123*exAtUU23*invgamt33 - 
  2*dgamt133*exAtUU33*invgamt33 + 2*dexAt113*invgamt13*invgamt33 + 
  2*dexAt123*invgamt23*invgamt33 + dexAt133*Power(invgamt33,2)
;

vreal dexAtUU211
=
-2*dgamt211*exAtUU11*invgamt11 - 2*dgamt213*exAtUU13*invgamt11 + 
  dexAt211*Power(invgamt11,2) - 2*dgamt222*exAtUU12*invgamt12 - 
  2*dgamt223*exAtUU13*invgamt12 + 2*dexAt212*invgamt11*invgamt12 + 
  dexAt222*Power(invgamt12,2) - 
  2*dgamt212*(exAtUU12*invgamt11 + exAtUU11*invgamt12) - 
  2*dgamt213*exAtUU11*invgamt13 - 2*dgamt223*exAtUU12*invgamt13 - 
  2*dgamt233*exAtUU13*invgamt13 + 2*dexAt213*invgamt11*invgamt13 + 
  2*dexAt223*invgamt12*invgamt13 + dexAt233*Power(invgamt13,2)
;

vreal dexAtUU212
=
-2*dgamt211*exAtUU12*invgamt11 - 2*dgamt213*exAtUU23*invgamt11 - 
  2*dgamt222*exAtUU22*invgamt12 - 2*dgamt223*exAtUU23*invgamt12 + 
  dexAt211*invgamt11*invgamt12 + dexAt212*Power(invgamt12,2) - 
  2*dgamt212*(exAtUU22*invgamt11 + exAtUU12*invgamt12) - 
  2*dgamt213*exAtUU12*invgamt13 - 2*dgamt223*exAtUU22*invgamt13 - 
  2*dgamt233*exAtUU23*invgamt13 + dexAt213*invgamt12*invgamt13 + 
  dexAt212*invgamt11*invgamt22 + dexAt222*invgamt12*invgamt22 + 
  dexAt223*invgamt13*invgamt22 + dexAt213*invgamt11*invgamt23 + 
  dexAt223*invgamt12*invgamt23 + dexAt233*invgamt13*invgamt23
;

vreal dexAtUU213
=
-2*dgamt211*exAtUU13*invgamt11 - 2*dgamt213*exAtUU33*invgamt11 - 
  2*dgamt222*exAtUU23*invgamt12 - 2*dgamt223*exAtUU33*invgamt12 - 
  2*dgamt212*(exAtUU23*invgamt11 + exAtUU13*invgamt12) - 
  2*dgamt213*exAtUU13*invgamt13 - 2*dgamt223*exAtUU23*invgamt13 - 
  2*dgamt233*exAtUU33*invgamt13 + dexAt211*invgamt11*invgamt13 + 
  dexAt212*invgamt12*invgamt13 + dexAt213*Power(invgamt13,2) + 
  dexAt212*invgamt11*invgamt23 + dexAt222*invgamt12*invgamt23 + 
  dexAt223*invgamt13*invgamt23 + dexAt213*invgamt11*invgamt33 + 
  dexAt223*invgamt12*invgamt33 + dexAt233*invgamt13*invgamt33
;

vreal dexAtUU222
=
-2*dgamt211*exAtUU12*invgamt12 - 2*dgamt213*exAtUU23*invgamt12 + 
  dexAt211*Power(invgamt12,2) - 2*dgamt222*exAtUU22*invgamt22 - 
  2*dgamt223*exAtUU23*invgamt22 + 2*dexAt212*invgamt12*invgamt22 + 
  dexAt222*Power(invgamt22,2) - 
  2*dgamt212*(exAtUU22*invgamt12 + exAtUU12*invgamt22) - 
  2*dgamt213*exAtUU12*invgamt23 - 2*dgamt223*exAtUU22*invgamt23 - 
  2*dgamt233*exAtUU23*invgamt23 + 2*dexAt213*invgamt12*invgamt23 + 
  2*dexAt223*invgamt22*invgamt23 + dexAt233*Power(invgamt23,2)
;

vreal dexAtUU223
=
-2*dgamt211*exAtUU13*invgamt12 - 2*dgamt213*exAtUU33*invgamt12 + 
  dexAt211*invgamt12*invgamt13 - 2*dgamt222*exAtUU23*invgamt22 - 
  2*dgamt223*exAtUU33*invgamt22 + dexAt212*invgamt13*invgamt22 - 
  2*dgamt212*(exAtUU23*invgamt12 + exAtUU13*invgamt22) - 
  2*dgamt213*exAtUU13*invgamt23 - 2*dgamt223*exAtUU23*invgamt23 - 
  2*dgamt233*exAtUU33*invgamt23 + dexAt212*invgamt12*invgamt23 + 
  dexAt213*invgamt13*invgamt23 + dexAt222*invgamt22*invgamt23 + 
  dexAt223*Power(invgamt23,2) + dexAt213*invgamt12*invgamt33 + 
  dexAt223*invgamt22*invgamt33 + dexAt233*invgamt23*invgamt33
;

vreal dexAtUU233
=
-2*dgamt211*exAtUU13*invgamt13 - 2*dgamt213*exAtUU33*invgamt13 + 
  dexAt211*Power(invgamt13,2) - 2*dgamt222*exAtUU23*invgamt23 - 
  2*dgamt223*exAtUU33*invgamt23 + 2*dexAt212*invgamt13*invgamt23 + 
  dexAt222*Power(invgamt23,2) - 
  2*dgamt212*(exAtUU23*invgamt13 + exAtUU13*invgamt23) - 
  2*dgamt213*exAtUU13*invgamt33 - 2*dgamt223*exAtUU23*invgamt33 - 
  2*dgamt233*exAtUU33*invgamt33 + 2*dexAt213*invgamt13*invgamt33 + 
  2*dexAt223*invgamt23*invgamt33 + dexAt233*Power(invgamt33,2)
;

vreal dexAtUU311
=
-2*dgamt311*exAtUU11*invgamt11 - 2*dgamt313*exAtUU13*invgamt11 + 
  dexAt311*Power(invgamt11,2) - 2*dgamt322*exAtUU12*invgamt12 - 
  2*dgamt323*exAtUU13*invgamt12 + 2*dexAt312*invgamt11*invgamt12 + 
  dexAt322*Power(invgamt12,2) - 
  2*dgamt312*(exAtUU12*invgamt11 + exAtUU11*invgamt12) - 
  2*dgamt313*exAtUU11*invgamt13 - 2*dgamt323*exAtUU12*invgamt13 - 
  2*dgamt333*exAtUU13*invgamt13 + 2*dexAt313*invgamt11*invgamt13 + 
  2*dexAt323*invgamt12*invgamt13 + dexAt333*Power(invgamt13,2)
;

vreal dexAtUU312
=
-2*dgamt311*exAtUU12*invgamt11 - 2*dgamt313*exAtUU23*invgamt11 - 
  2*dgamt322*exAtUU22*invgamt12 - 2*dgamt323*exAtUU23*invgamt12 + 
  dexAt311*invgamt11*invgamt12 + dexAt312*Power(invgamt12,2) - 
  2*dgamt312*(exAtUU22*invgamt11 + exAtUU12*invgamt12) - 
  2*dgamt313*exAtUU12*invgamt13 - 2*dgamt323*exAtUU22*invgamt13 - 
  2*dgamt333*exAtUU23*invgamt13 + dexAt313*invgamt12*invgamt13 + 
  dexAt312*invgamt11*invgamt22 + dexAt322*invgamt12*invgamt22 + 
  dexAt323*invgamt13*invgamt22 + dexAt313*invgamt11*invgamt23 + 
  dexAt323*invgamt12*invgamt23 + dexAt333*invgamt13*invgamt23
;

vreal dexAtUU313
=
-2*dgamt311*exAtUU13*invgamt11 - 2*dgamt313*exAtUU33*invgamt11 - 
  2*dgamt322*exAtUU23*invgamt12 - 2*dgamt323*exAtUU33*invgamt12 - 
  2*dgamt312*(exAtUU23*invgamt11 + exAtUU13*invgamt12) - 
  2*dgamt313*exAtUU13*invgamt13 - 2*dgamt323*exAtUU23*invgamt13 - 
  2*dgamt333*exAtUU33*invgamt13 + dexAt311*invgamt11*invgamt13 + 
  dexAt312*invgamt12*invgamt13 + dexAt313*Power(invgamt13,2) + 
  dexAt312*invgamt11*invgamt23 + dexAt322*invgamt12*invgamt23 + 
  dexAt323*invgamt13*invgamt23 + dexAt313*invgamt11*invgamt33 + 
  dexAt323*invgamt12*invgamt33 + dexAt333*invgamt13*invgamt33
;

vreal dexAtUU322
=
-2*dgamt311*exAtUU12*invgamt12 - 2*dgamt313*exAtUU23*invgamt12 + 
  dexAt311*Power(invgamt12,2) - 2*dgamt322*exAtUU22*invgamt22 - 
  2*dgamt323*exAtUU23*invgamt22 + 2*dexAt312*invgamt12*invgamt22 + 
  dexAt322*Power(invgamt22,2) - 
  2*dgamt312*(exAtUU22*invgamt12 + exAtUU12*invgamt22) - 
  2*dgamt313*exAtUU12*invgamt23 - 2*dgamt323*exAtUU22*invgamt23 - 
  2*dgamt333*exAtUU23*invgamt23 + 2*dexAt313*invgamt12*invgamt23 + 
  2*dexAt323*invgamt22*invgamt23 + dexAt333*Power(invgamt23,2)
;

vreal dexAtUU323
=
-2*dgamt311*exAtUU13*invgamt12 - 2*dgamt313*exAtUU33*invgamt12 + 
  dexAt311*invgamt12*invgamt13 - 2*dgamt322*exAtUU23*invgamt22 - 
  2*dgamt323*exAtUU33*invgamt22 + dexAt312*invgamt13*invgamt22 - 
  2*dgamt312*(exAtUU23*invgamt12 + exAtUU13*invgamt22) - 
  2*dgamt313*exAtUU13*invgamt23 - 2*dgamt323*exAtUU23*invgamt23 - 
  2*dgamt333*exAtUU33*invgamt23 + dexAt312*invgamt12*invgamt23 + 
  dexAt313*invgamt13*invgamt23 + dexAt322*invgamt22*invgamt23 + 
  dexAt323*Power(invgamt23,2) + dexAt313*invgamt12*invgamt33 + 
  dexAt323*invgamt22*invgamt33 + dexAt333*invgamt23*invgamt33
;

vreal dexAtUU333
=
-2*dgamt311*exAtUU13*invgamt13 - 2*dgamt313*exAtUU33*invgamt13 + 
  dexAt311*Power(invgamt13,2) - 2*dgamt322*exAtUU23*invgamt23 - 
  2*dgamt323*exAtUU33*invgamt23 + 2*dexAt312*invgamt13*invgamt23 + 
  dexAt322*Power(invgamt23,2) - 
  2*dgamt312*(exAtUU23*invgamt13 + exAtUU13*invgamt23) - 
  2*dgamt313*exAtUU13*invgamt33 - 2*dgamt323*exAtUU23*invgamt33 - 
  2*dgamt333*exAtUU33*invgamt33 + 2*dexAt313*invgamt13*invgamt33 + 
  2*dexAt323*invgamt23*invgamt33 + dexAt333*Power(invgamt33,2)
;


vreal ZtC1
=
trGt1 - trGtd1
;

vreal ZtC2
=
trGt2 - trGtd2
;

vreal ZtC3
=
trGt3 - trGtd3
;

vreal HC
=
exAt11*exAtUU11 + 2*exAt12*exAtUU12 + 2*exAt13*exAtUU13 + exAt22*exAtUU22 + 
  2*exAt23*exAtUU23 + exAt33*exAtUU33 - (2*Power(exKh,2))/3. - 16*cpi*rho - 
  (8*exKh*Theta)/3. - (8*Power(Theta,2))/3. + trR
;

vreal MtC1
=
(-2*(dchi1*exAtUU11 + dchi2*exAtUU12 + dchi3*exAtUU13) + 
    chi*(3*dexAtUU111 + 3*dexAtUU212 + 3*dexAtUU313 + 3*exAtUU11*Gt111 + 
       6*exAtUU12*Gt112 + 6*exAtUU13*Gt113 + 3*exAtUU22*Gt122 + 
       6*exAtUU23*Gt123 + 3*exAtUU33*Gt133 - 2*dexKh1*invgamt11 - 
       4*dTheta1*invgamt11 - 2*dexKh2*invgamt12 - 4*dTheta2*invgamt12 - 
       2*dexKh3*invgamt13 - 4*dTheta3*invgamt13 - 24*cpi*invgamt11*Sm1 - 
       24*cpi*invgamt12*Sm2 - 24*cpi*invgamt13*Sm3))/(3.*chi)
;

vreal MtC2
=
(-2*(dchi1*exAtUU12 + dchi2*exAtUU22 + dchi3*exAtUU23) + 
    chi*(3*dexAtUU112 + 3*dexAtUU222 + 3*dexAtUU323 + 3*exAtUU11*Gt211 + 
       6*exAtUU12*Gt212 + 6*exAtUU13*Gt213 + 3*exAtUU22*Gt222 + 
       6*exAtUU23*Gt223 + 3*exAtUU33*Gt233 - 2*dexKh1*invgamt12 - 
       4*dTheta1*invgamt12 - 2*dexKh2*invgamt22 - 4*dTheta2*invgamt22 - 
       2*dexKh3*invgamt23 - 4*dTheta3*invgamt23 - 24*cpi*invgamt12*Sm1 - 
       24*cpi*invgamt22*Sm2 - 24*cpi*invgamt23*Sm3))/(3.*chi)
;

vreal MtC3
=
(-2*(dchi1*exAtUU13 + dchi2*exAtUU23 + dchi3*exAtUU33) + 
    chi*(3*dexAtUU113 + 3*dexAtUU223 + 3*dexAtUU333 + 3*exAtUU11*Gt311 + 
       6*exAtUU12*Gt312 + 6*exAtUU13*Gt313 + 3*exAtUU22*Gt322 + 
       6*exAtUU23*Gt323 + 3*exAtUU33*Gt333 - 2*dexKh1*invgamt13 - 
       4*dTheta1*invgamt13 - 2*dexKh2*invgamt23 - 4*dTheta2*invgamt23 - 
       2*dexKh3*invgamt33 - 4*dTheta3*invgamt33 - 24*cpi*invgamt13*Sm1 - 
       24*cpi*invgamt23*Sm2 - 24*cpi*invgamt33*Sm3))/(3.*chi)
;

dtchi
=
(-2*chi*(dbeta11 + dbeta22 + dbeta33 - alpha*exKh - 2*alpha*Theta))/3.
;

dtgamt11
=
-2*alpha*exAt11 + 2*dbeta11*gamt11 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt11)/3. + 2*dbeta12*gamt12 + 
  2*dbeta13*gamt13
;

dtgamt12
=
-2*alpha*exAt12 + dbeta21*gamt11 + dbeta11*gamt12 + dbeta22*gamt12 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt12)/3. + dbeta23*gamt13 + 
  dbeta12*gamt22 + dbeta13*gamt23
;

dtgamt13
=
-2*alpha*exAt13 + dbeta31*gamt11 + dbeta32*gamt12 + dbeta11*gamt13 + 
  dbeta33*gamt13 - (2*(dbeta11 + dbeta22 + dbeta33)*gamt13)/3. + 
  dbeta12*gamt23 + dbeta13*gamt33
;

dtgamt22
=
-2*alpha*exAt22 + 2*dbeta21*gamt12 + 2*dbeta22*gamt22 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt22)/3. + 2*dbeta23*gamt23
;

dtgamt23
=
-2*alpha*exAt23 + dbeta31*gamt12 + dbeta21*gamt13 + dbeta32*gamt22 + 
  dbeta22*gamt23 + dbeta33*gamt23 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt23)/3. + dbeta23*gamt33
;

dtgamt33
=
-2*alpha*exAt33 + 2*dbeta31*gamt13 + 2*dbeta32*gamt23 + 2*dbeta33*gamt33 - 
  (2*(dbeta11 + dbeta22 + dbeta33)*gamt33)/3.
;

dtexKh
=
-(DDalpha11*invgam11) - 2*DDalpha12*invgam12 - 2*DDalpha13*invgam13 - 
  DDalpha22*invgam22 - 2*DDalpha23*invgam23 - DDalpha33*invgam33 + 
  alpha*(exAt11*exAtUU11 + 2*exAt12*exAtUU12 + 2*exAt13*exAtUU13 + 
     exAt22*exAtUU22 + 2*exAt23*exAtUU23 + exAt33*exAtUU33 + 
     Power(exKh,2)/3. + 4*cpi*rho + ckappa1*Theta - ckappa1*ckappa2*Theta + 
     (4*exKh*Theta)/3. + (4*Power(Theta,2))/3. + 4*cpi*trSs)
;

dtexAt11
=
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
;

dtexAt12
=
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
;

dtexAt13
=
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
;

dtexAt22
=
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
;

dtexAt23
=
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
;

dtexAt33
=
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
;

dttrGt1
=
(-6*dalpha1*exAtUU11 - 6*dalpha2*exAtUU12 - 6*dalpha3*exAtUU13 - 
    (9*alpha*(dchi1*exAtUU11 + dchi2*exAtUU12 + dchi3*exAtUU13))/chi + 
    4*ddbeta111*invgamt11 + ddbeta122*invgamt11 + ddbeta133*invgamt11 + 
    7*ddbeta121*invgamt12 + ddbeta222*invgamt12 + ddbeta233*invgamt12 + 
    7*ddbeta131*invgamt13 + ddbeta232*invgamt13 + ddbeta333*invgamt13 + 
    3*ddbeta221*invgamt22 + 6*ddbeta231*invgamt23 + 3*ddbeta331*invgamt33 - 
    2*dbeta11*trGtd1 - 2*dbeta22*trGtd1 - 2*dbeta33*trGtd1 + 
    2*alpha*(3*exAtUU11*Gt111 + 6*exAtUU12*Gt112 + 6*exAtUU13*Gt113 + 
       3*exAtUU22*Gt122 + 6*exAtUU23*Gt123 + 3*exAtUU33*Gt133 - 
       2*dexKh1*invgamt11 - dTheta1*invgamt11 - 2*dexKh2*invgamt12 - 
       dTheta2*invgamt12 - 2*dexKh3*invgamt13 - dTheta3*invgamt13 - 
       24*cpi*invgam11*Sm1 - 24*cpi*invgam12*Sm2 - 24*cpi*invgam13*Sm3 - 
       3*ckappa1*trGt1 + 3*ckappa1*trGtd1))/3.
;

dttrGt2
=
(-6*dalpha1*exAtUU12 - 6*dalpha2*exAtUU22 - 6*dalpha3*exAtUU23 - 
    (9*alpha*(dchi1*exAtUU12 + dchi2*exAtUU22 + dchi3*exAtUU23))/chi + 
    3*ddbeta112*invgamt11 + ddbeta111*invgamt12 + 7*ddbeta122*invgamt12 + 
    ddbeta133*invgamt12 + 6*ddbeta132*invgamt13 + ddbeta121*invgamt22 + 
    4*ddbeta222*invgamt22 + ddbeta233*invgamt22 + ddbeta131*invgamt23 + 
    7*ddbeta232*invgamt23 + ddbeta333*invgamt23 + 3*ddbeta332*invgamt33 - 
    2*dbeta11*trGtd2 - 2*dbeta22*trGtd2 - 2*dbeta33*trGtd2 + 
    2*alpha*(3*exAtUU11*Gt211 + 6*exAtUU12*Gt212 + 6*exAtUU13*Gt213 + 
       3*exAtUU22*Gt222 + 6*exAtUU23*Gt223 + 3*exAtUU33*Gt233 - 
       2*dexKh1*invgamt12 - dTheta1*invgamt12 - 2*dexKh2*invgamt22 - 
       dTheta2*invgamt22 - 2*dexKh3*invgamt23 - dTheta3*invgamt23 - 
       24*cpi*invgam12*Sm1 - 24*cpi*invgam22*Sm2 - 24*cpi*invgam23*Sm3 - 
       3*ckappa1*trGt2 + 3*ckappa1*trGtd2))/3.
;

dttrGt3
=
(-6*dalpha1*exAtUU13 - 6*dalpha2*exAtUU23 - 6*dalpha3*exAtUU33 - 
    (9*alpha*(dchi1*exAtUU13 + dchi2*exAtUU23 + dchi3*exAtUU33))/chi + 
    3*ddbeta113*invgamt11 + 6*ddbeta123*invgamt12 + ddbeta111*invgamt13 + 
    ddbeta122*invgamt13 + 7*ddbeta133*invgamt13 + 3*ddbeta223*invgamt22 + 
    ddbeta121*invgamt23 + ddbeta222*invgamt23 + 7*ddbeta233*invgamt23 + 
    ddbeta131*invgamt33 + ddbeta232*invgamt33 + 4*ddbeta333*invgamt33 - 
    2*dbeta11*trGtd3 - 2*dbeta22*trGtd3 - 2*dbeta33*trGtd3 + 
    2*alpha*(3*exAtUU11*Gt311 + 6*exAtUU12*Gt312 + 6*exAtUU13*Gt313 + 
       3*exAtUU22*Gt322 + 6*exAtUU23*Gt323 + 3*exAtUU33*Gt333 - 
       2*dexKh1*invgamt13 - dTheta1*invgamt13 - 2*dexKh2*invgamt23 - 
       dTheta2*invgamt23 - 2*dexKh3*invgamt33 - dTheta3*invgamt33 - 
       24*cpi*invgam13*Sm1 - 24*cpi*invgam23*Sm2 - 24*cpi*invgam33*Sm3 - 
       3*ckappa1*trGt3 + 3*ckappa1*trGtd3))/3.
;

dtTheta
=
-0.16666666666666666*(alpha*(3*exAt11*exAtUU11 + 6*exAt12*exAtUU12 + 
      6*exAt13*exAtUU13 + 3*exAt22*exAtUU22 + 6*exAt23*exAtUU23 + 
      3*exAt33*exAtUU33 - 2*Power(exKh,2) + 48*cpi*rho + 12*ckappa1*Theta + 
      6*ckappa1*ckappa2*Theta - 8*exKh*Theta - 8*Power(Theta,2) - 3*trR))
;

dtalpha
=
-(Power(alpha,2)*cmuL*exKh)
;

dtbeta1
=
-(beta1*ceta) + Power(alpha,2)*cmuS*trGt1
;

dtbeta2
=
-(beta2*ceta) + Power(alpha,2)*cmuS*trGt2
;

dtbeta3
=
-(beta3*ceta) + Power(alpha,2)*cmuS*trGt3
;


/* Z4co_set_rhs.hxx */
