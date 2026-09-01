// Apply NewRadX boundary conditions to Z4c variables
//=============================================================================

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <newradx.hxx>

using namespace NewRadX;

namespace Z4c {

extern "C" void Z4c_apply_newradx_boundary_conditions(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Z4c_apply_newradx_boundary_conditions;
  DECLARE_CCTK_PARAMETERS;

  // The fourth argument is the value at infinity. Every Z4c variable is
  // stored as a deviation that tends to zero in flat space (see
  // interface.ccl), so it is zero for all of them -- including chi, the
  // diagonal of gamma_tilde and alphaG, which used to pass one. Passing one
  // turned the -v0 (var - var0) / r term of the radiative condition into a
  // spurious +v0/r source that drove those variables towards one at the
  // outer boundary.
  NewRadX_Apply(cctkGH, chi, chi_rhs, 0, 1, n_chi);
  NewRadX_Apply(cctkGH, gammatxx, gammatxx_rhs, 0, 1, n_gammat);
  NewRadX_Apply(cctkGH, gammatxy, gammatxy_rhs, 0, 1, n_gammat);
  NewRadX_Apply(cctkGH, gammatxz, gammatxz_rhs, 0, 1, n_gammat);
  NewRadX_Apply(cctkGH, gammatyy, gammatyy_rhs, 0, 1, n_gammat);
  NewRadX_Apply(cctkGH, gammatyz, gammatyz_rhs, 0, 1, n_gammat);
  NewRadX_Apply(cctkGH, gammatzz, gammatzz_rhs, 0, 1, n_gammat);
  NewRadX_Apply(cctkGH, Kh, Kh_rhs, 0, 1, n_Kh);
  NewRadX_Apply(cctkGH, Atxx, Atxx_rhs, 0, 1, n_At);
  NewRadX_Apply(cctkGH, Atxy, Atxy_rhs, 0, 1, n_At);
  NewRadX_Apply(cctkGH, Atxz, Atxz_rhs, 0, 1, n_At);
  NewRadX_Apply(cctkGH, Atyy, Atyy_rhs, 0, 1, n_At);
  NewRadX_Apply(cctkGH, Atyz, Atyz_rhs, 0, 1, n_At);
  NewRadX_Apply(cctkGH, Atzz, Atzz_rhs, 0, 1, n_At);
  NewRadX_Apply(cctkGH, Gamtx, Gamtx_rhs, 0, 1, n_Gamt);
  NewRadX_Apply(cctkGH, Gamty, Gamty_rhs, 0, 1, n_Gamt);
  NewRadX_Apply(cctkGH, Gamtz, Gamtz_rhs, 0, 1, n_Gamt);
  NewRadX_Apply(cctkGH, Theta, Theta_rhs, 0, 1, n_Theta);
  NewRadX_Apply(cctkGH, alphaG, alphaG_rhs, 0, 1, n_alphaG);
  NewRadX_Apply(cctkGH, betaGx, betaGx_rhs, 0, 1, n_betaG);
  NewRadX_Apply(cctkGH, betaGy, betaGy_rhs, 0, 1, n_betaG);
  NewRadX_Apply(cctkGH, betaGz, betaGz_rhs, 0, 1, n_betaG);
  if (evolveA)
    NewRadX_Apply(cctkGH, A, A_rhs, 0, 1, n_A);
  if (evolveB) {
    NewRadX_Apply(cctkGH, Bx, Bx_rhs, 0, 1, n_B);
    NewRadX_Apply(cctkGH, By, By_rhs, 0, 1, n_B);
    NewRadX_Apply(cctkGH, Bz, Bz_rhs, 0, 1, n_B);
  }
}

} // namespace Z4c
