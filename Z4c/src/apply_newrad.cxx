// Apply NewRadX boundary conditions to Z4c variables
//=============================================================================

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <newradx.hxx>

using namespace NewRadX;

namespace Z4c {

extern "C" void apply_boundary_conditions(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_apply_boundary_conditions;
  DECLARE_CCTK_PARAMETERS;

  NewRadX_Apply(cctkGH, chi, chi_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, gammatxx, gammatxx_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, gammatxy, gammatxy_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, gammatxz, gammatxz_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, gammatyy, gammatyy_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, gammatyz, gammatyz_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, gammatzz, gammatzz_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, Kh, Kh_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, Atxx, Atxx_rhs, 0, 1, 2);
  NewRadX_Apply(cctkGH, Atxy, Atxy_rhs, 0, 1, 2);
  NewRadX_Apply(cctkGH, Atxz, Atxz_rhs, 0, 1, 2);
  NewRadX_Apply(cctkGH, Atyy, Atyy_rhs, 0, 1, 2);
  NewRadX_Apply(cctkGH, Atyz, Atyz_rhs, 0, 1, 2);
  NewRadX_Apply(cctkGH, Atzz, Atzz_rhs, 0, 1, 2);
  NewRadX_Apply(cctkGH, Gamtx, Gamtx_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, Gamty, Gamty_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, Gamtz, Gamtz_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, Theta, Theta_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, alphaG, alphaG_rhs, 1, 1, 1);
  NewRadX_Apply(cctkGH, betaGx, betaGx_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, betaGy, betaGy_rhs, 0, 1, 1);
  NewRadX_Apply(cctkGH, betaGz, betaGz_rhs, 0, 1, 1);
}

} // namespace Z4c
