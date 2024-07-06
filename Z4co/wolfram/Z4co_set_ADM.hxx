/* Z4co_set_ADM.hxx */
/* Produced with Mathematica */

const vreal &chi = gf_chi(mask, index2, 1);
const vreal &gamt11 = gf_gamt(mask, index2, 1)(0,0);
const vreal &gamt12 = gf_gamt(mask, index2, 0)(0,1);
const vreal &gamt13 = gf_gamt(mask, index2, 0)(0,2);
const vreal &gamt22 = gf_gamt(mask, index2, 1)(1,1);
const vreal &gamt23 = gf_gamt(mask, index2, 0)(1,2);
const vreal &gamt33 = gf_gamt(mask, index2, 1)(2,2);
const vreal &exKh = gf_exKh(mask, index2);
const vreal &exAt11 = gf_exAt(mask, index2)(0,0);
const vreal &exAt12 = gf_exAt(mask, index2)(0,1);
const vreal &exAt13 = gf_exAt(mask, index2)(0,2);
const vreal &exAt22 = gf_exAt(mask, index2)(1,1);
const vreal &exAt23 = gf_exAt(mask, index2)(1,2);
const vreal &exAt33 = gf_exAt(mask, index2)(2,2);
const vreal &Theta = gf_Theta(mask, index2);
const vreal &alpha = gf_alpha(mask, index2, 1);
const vreal &beta1 = gf_beta(mask, index2)(0);
const vreal &beta2 = gf_beta(mask, index2)(1);
const vreal &beta3 = gf_beta(mask, index2)(2);

ADMgam11.store(mask, index2, 
gamt11/chi
);

ADMgam12.store(mask, index2, 
gamt12/chi
);

ADMgam13.store(mask, index2, 
gamt13/chi
);

ADMgam22.store(mask, index2, 
gamt22/chi
);

ADMgam23.store(mask, index2, 
gamt23/chi
);

ADMgam33.store(mask, index2, 
gamt33/chi
);

ADMK11.store(mask, index2, 
(exAt11 + (gamt11*(exKh + 2*Theta))/3.)/chi
);

ADMK12.store(mask, index2, 
(exAt12 + (gamt12*(exKh + 2*Theta))/3.)/chi
);

ADMK13.store(mask, index2, 
(exAt13 + (gamt13*(exKh + 2*Theta))/3.)/chi
);

ADMK22.store(mask, index2, 
(exAt22 + (gamt22*(exKh + 2*Theta))/3.)/chi
);

ADMK23.store(mask, index2, 
(exAt23 + (gamt23*(exKh + 2*Theta))/3.)/chi
);

ADMK33.store(mask, index2, 
(exAt33 + (gamt33*(exKh + 2*Theta))/3.)/chi
);

ADMalpha.store(mask, index2, 
alpha
);

ADMbeta1.store(mask, index2, 
beta1
);

ADMbeta2.store(mask, index2, 
beta2
);

ADMbeta3.store(mask, index2, 
beta3
);


/* Z4co_set_ADM.hxx */
