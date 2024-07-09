/* Z4co_set_ADM.hxx */
/* Produced with Mathematica */

const GF3D2<CCTK_REAL> &local_ADMgam11 = gf_ADMgam(0,0);
const GF3D2<CCTK_REAL> &local_ADMgam12 = gf_ADMgam(0,1);
const GF3D2<CCTK_REAL> &local_ADMgam13 = gf_ADMgam(0,2);
const GF3D2<CCTK_REAL> &local_ADMgam22 = gf_ADMgam(1,1);
const GF3D2<CCTK_REAL> &local_ADMgam23 = gf_ADMgam(1,2);
const GF3D2<CCTK_REAL> &local_ADMgam33 = gf_ADMgam(2,2);
const GF3D2<CCTK_REAL> &local_ADMK11 = gf_ADMK(0,0);
const GF3D2<CCTK_REAL> &local_ADMK12 = gf_ADMK(0,1);
const GF3D2<CCTK_REAL> &local_ADMK13 = gf_ADMK(0,2);
const GF3D2<CCTK_REAL> &local_ADMK22 = gf_ADMK(1,1);
const GF3D2<CCTK_REAL> &local_ADMK23 = gf_ADMK(1,2);
const GF3D2<CCTK_REAL> &local_ADMK33 = gf_ADMK(2,2);
const GF3D2<CCTK_REAL> &local_ADMalpha = gf_ADMalpha;
const GF3D2<CCTK_REAL> &local_ADMbeta1 = gf_ADMbeta(0);
const GF3D2<CCTK_REAL> &local_ADMbeta2 = gf_ADMbeta(1);
const GF3D2<CCTK_REAL> &local_ADMbeta3 = gf_ADMbeta(2);

grid.loop_all_device<0, 0, 0, vsize>(
  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
  const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
  const GF3D2index index2(layout2, p.I);

const auto &tmp_chi = gf_chi(mask, index2);
const auto &tmp_gamt = gf_gamt(mask, index2);
const auto &tmp_exKh = gf_exKh(mask, index2);
const auto &tmp_exAt = gf_exAt(mask, index2);
const auto &tmp_Theta = gf_Theta(mask, index2);
const auto &tmp_alpha = gf_alpha(mask, index2);
const auto &tmp_beta = gf_beta(mask, index2);

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
const vreal Theta = tmp_Theta;
const vreal alpha = tmp_alpha;
const vreal beta1 = tmp_beta(0);
const vreal beta2 = tmp_beta(1);
const vreal beta3 = tmp_beta(2);

local_ADMgam11.store(mask, index2, 
gamt11/chi
);

local_ADMgam12.store(mask, index2, 
gamt12/chi
);

local_ADMgam13.store(mask, index2, 
gamt13/chi
);

local_ADMgam22.store(mask, index2, 
gamt22/chi
);

local_ADMgam23.store(mask, index2, 
gamt23/chi
);

local_ADMgam33.store(mask, index2, 
gamt33/chi
);

local_ADMK11.store(mask, index2, 
(exAt11 + (gamt11*(exKh + 2*Theta))/3.)/chi
);

local_ADMK12.store(mask, index2, 
(exAt12 + (gamt12*(exKh + 2*Theta))/3.)/chi
);

local_ADMK13.store(mask, index2, 
(exAt13 + (gamt13*(exKh + 2*Theta))/3.)/chi
);

local_ADMK22.store(mask, index2, 
(exAt22 + (gamt22*(exKh + 2*Theta))/3.)/chi
);

local_ADMK23.store(mask, index2, 
(exAt23 + (gamt23*(exKh + 2*Theta))/3.)/chi
);

local_ADMK33.store(mask, index2, 
(exAt33 + (gamt33*(exKh + 2*Theta))/3.)/chi
);

local_ADMalpha.store(mask, index2, 
alpha
);

local_ADMbeta1.store(mask, index2, 
beta1
);

local_ADMbeta2.store(mask, index2, 
beta2
);

local_ADMbeta3.store(mask, index2, 
beta3
);


});

/* Z4co_set_ADM.hxx */
