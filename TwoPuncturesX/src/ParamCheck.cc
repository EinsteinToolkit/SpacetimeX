/* TwoPuncturesX:  File  "TwoPuncturesX.c"*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

/* -------------------------------------------------------------------*/
extern "C"
void TwoPuncturesX_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (use_sources) {
    if (!CCTK_IsFunctionAliased("Set_Rho_ADM"))
      CCTK_WARN(0, "Matter sources have been enabled, "
                   "but there is no aliased function for matter sources.");
  }

  if (CCTK_EQUALS(initial_lapse, "psi^n")) {
    CCTK_VInfo(CCTK_THORNSTRING, "Setting initial lapse to psi^%f profile,",
               (double)initial_lapse_psi_exponent);
    CCTK_VInfo(CCTK_THORNSTRING, "where psi is the *uncorrected* conformal");
    CCTK_VInfo(CCTK_THORNSTRING, "factor, not psi+u!");
  }
  if (CCTK_EQUALS(initial_lapse, "W"))
    CCTK_VInfo(CCTK_THORNSTRING,
               "Setting initial lapse to W=(psi+u)^{-2}, needed for slow-start lapse.");
  if (CCTK_EQUALS(initial_lapse, "brownsville"))
    CCTK_VINFO("Setting initial lapse to a Brownsville-style profile "
               "with exp %f.",
               (double)initial_lapse_psi_exponent);
}
