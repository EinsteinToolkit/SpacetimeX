#include <cctk.h>

#ifdef __CUDACC__
// Disable CCTK_DEBUG since the debug information takes too much
// parameter space to launch the kernels
#ifdef CCTK_DEBUG
#undef CCTK_DEBUG
#endif
#endif

// #define Power(x, y) (Arith::pown((x), (y)))

#include "diss.hxx"

#include <derivs.hxx>
#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef __CUDACC__
#include <nvToolsExt.h>
#endif

#include <cmath>

namespace Z4co {
using namespace Arith;
using namespace Loop;
using namespace std;

template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T Power(T x, int y) {
  return (y == 2) ? Arith::pow2(x) : Arith::pown(x, y);
}

extern "C" void Z4co_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4co_RHS;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2 + 1)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2 + 1);

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};
  vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);

  // Suffix 1: with ghost zones, suffix 0: without ghost zones
  const GF3D2layout layout2(cctkGH, indextype);
  const GF3D5layout layout5(imin, imax);

  // Input grid functions
  const GF3D2<const CCTK_REAL> gf_chi(layout2, chi);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_gamt{
      GF3D2<const CCTK_REAL>(layout2, gammatxx),
      GF3D2<const CCTK_REAL>(layout2, gammatxy),
      GF3D2<const CCTK_REAL>(layout2, gammatxz),
      GF3D2<const CCTK_REAL>(layout2, gammatyy),
      GF3D2<const CCTK_REAL>(layout2, gammatyz),
      GF3D2<const CCTK_REAL>(layout2, gammatzz)};

  const GF3D2<const CCTK_REAL> gf_exKh(layout2, Kh);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_exAt{
      GF3D2<const CCTK_REAL>(layout2, Atxx),
      GF3D2<const CCTK_REAL>(layout2, Atxy),
      GF3D2<const CCTK_REAL>(layout2, Atxz),
      GF3D2<const CCTK_REAL>(layout2, Atyy),
      GF3D2<const CCTK_REAL>(layout2, Atyz),
      GF3D2<const CCTK_REAL>(layout2, Atzz)};

  const vec<GF3D2<const CCTK_REAL>, 3> gf_trGt{
      GF3D2<const CCTK_REAL>(layout2, Gamtx),
      GF3D2<const CCTK_REAL>(layout2, Gamty),
      GF3D2<const CCTK_REAL>(layout2, Gamtz)};

  const GF3D2<const CCTK_REAL> gf_Theta(layout2, Theta);

  const GF3D2<const CCTK_REAL> gf_alpha(layout2, alphaG);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{
      GF3D2<const CCTK_REAL>(layout2, betaGx),
      GF3D2<const CCTK_REAL>(layout2, betaGy),
      GF3D2<const CCTK_REAL>(layout2, betaGz)};

  // Tile variables for derivatives and so on
  const int ntmps = 136;
  GF3D5vector<CCTK_REAL> tmps(layout5, ntmps);
  int itmp = 0;

  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(tmps(itmp++)); };
  const auto make_vec = [&](const auto &f) {
    return vec<result_of_t<decltype(f)()>, 3>([&](int) { return f(); });
  };
  const auto make_mat = [&](const auto &f) {
    return smat<result_of_t<decltype(f)()>, 3>([&](int, int) { return f(); });
  };
  const auto make_vec_gf = [&]() { return make_vec(make_gf); };
  const auto make_mat_gf = [&]() { return make_mat(make_gf); };
  const auto make_vec_vec_gf = [&]() { return make_vec(make_vec_gf); };
  const auto make_vec_mat_gf = [&]() { return make_vec(make_mat_gf); };
  const auto make_mat_vec_gf = [&]() { return make_mat(make_vec_gf); };
  const auto make_mat_mat_gf = [&]() { return make_mat(make_mat_gf); };

  const GF3D5<CCTK_REAL> tl_chi(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dchi(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3> tl_ddchi(make_mat_gf());

  const Loop::GridDescBaseDevice grid(cctkGH);

  const vect<CCTK_REAL, dim> dx(std::array<CCTK_REAL, dim>{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  });

  const auto calccopy = [&](const auto &gf, const auto &gf0) {
    Derivs::calc_copy<0, 0, 0>(gf, layout5, grid, gf0);
  };

  const auto calcderivs = [&](const auto &gf, const auto &dgf,
                              const auto &gf0) {
    Derivs::calc_derivs<0, 0, 0>(gf, dgf, layout5, grid, gf0, dx, deriv_order);
  };

  const auto calcderivs2 = [&](const auto &gf, const auto &dgf,
                               const auto &ddgf, const auto &gf0) {
    Derivs::calc_derivs2<0, 0, 0>(gf, dgf, ddgf, layout5, grid, gf0, dx,
                                  deriv_order);
  };

  calcderivs2(tl_chi, tl_dchi, tl_ddchi, gf_chi);

  const smat<GF3D5<CCTK_REAL>, 3> tl_gamt(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> tl_dgamt(make_mat_vec_gf());
  const smat<smat<GF3D5<CCTK_REAL>, 3>, 3> tl_ddgamt(make_mat_mat_gf());
  calcderivs2(tl_gamt, tl_dgamt, tl_ddgamt, gf_gamt);

  const GF3D5<CCTK_REAL> tl_exKh(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dexKh(make_vec_gf());
  calcderivs(tl_exKh, tl_dexKh, gf_exKh);

  const smat<GF3D5<CCTK_REAL>, 3> tl_exAt(make_mat_gf());
  calccopy(tl_exAt, gf_exAt);

  const vec<GF3D5<CCTK_REAL>, 3> tl_trGt(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> tl_dtrGt(make_vec_vec_gf());
  calcderivs(tl_trGt, tl_dtrGt, gf_trGt);

  const GF3D5<CCTK_REAL> tl_Theta(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dTheta(make_vec_gf());
  calcderivs(tl_Theta, tl_dTheta, gf_Theta);

  const GF3D5<CCTK_REAL> tl_alpha(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dalpha(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3> tl_ddalpha(make_mat_gf());
  calcderivs2(tl_alpha, tl_dalpha, tl_ddalpha, gf_alpha);

  const vec<GF3D5<CCTK_REAL>, 3> tl_beta(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> tl_dbeta(make_vec_vec_gf());
  const vec<smat<GF3D5<CCTK_REAL>, 3>, 3> tl_ddbeta(make_vec_mat_gf());
  calcderivs2(tl_beta, tl_dbeta, tl_ddbeta, gf_beta);

  if (itmp != ntmps)
    CCTK_VERROR("Wrong number of temporary variables: ntmps=%d itmp=%d", ntmps,
                itmp);
  itmp = -1;

  // More input grid functions
  const GF3D2<const CCTK_REAL> gf_eTtt(layout2, eTtt);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_eTt{
      GF3D2<const CCTK_REAL>(layout2, eTtx),
      GF3D2<const CCTK_REAL>(layout2, eTty),
      GF3D2<const CCTK_REAL>(layout2, eTtz)};

  const smat<GF3D2<const CCTK_REAL>, 3> gf_eT{
      GF3D2<const CCTK_REAL>(layout2, eTxx),
      GF3D2<const CCTK_REAL>(layout2, eTxy),
      GF3D2<const CCTK_REAL>(layout2, eTxz),
      GF3D2<const CCTK_REAL>(layout2, eTyy),
      GF3D2<const CCTK_REAL>(layout2, eTyz),
      GF3D2<const CCTK_REAL>(layout2, eTzz)};

  // Output grid functions
  const GF3D2<CCTK_REAL> gf_dtchi(layout2, chi_rhs);

  const smat<GF3D2<CCTK_REAL>, 3> gf_dtgamt{
      GF3D2<CCTK_REAL>(layout2, gammatxx_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatxy_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatxz_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatyy_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatyz_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatzz_rhs)};

  const GF3D2<CCTK_REAL> gf_dtexKh(layout2, Kh_rhs);

  const smat<GF3D2<CCTK_REAL>, 3> gf_dtexAt{
      GF3D2<CCTK_REAL>(layout2, Atxx_rhs), GF3D2<CCTK_REAL>(layout2, Atxy_rhs),
      GF3D2<CCTK_REAL>(layout2, Atxz_rhs), GF3D2<CCTK_REAL>(layout2, Atyy_rhs),
      GF3D2<CCTK_REAL>(layout2, Atyz_rhs), GF3D2<CCTK_REAL>(layout2, Atzz_rhs)};

  const vec<GF3D2<CCTK_REAL>, 3> gf_dttrGt{
      GF3D2<CCTK_REAL>(layout2, Gamtx_rhs),
      GF3D2<CCTK_REAL>(layout2, Gamty_rhs),
      GF3D2<CCTK_REAL>(layout2, Gamtz_rhs)};

  const GF3D2<CCTK_REAL> gf_dtTheta(layout2, Theta_rhs);

  const GF3D2<CCTK_REAL> gf_dtalpha(layout2, alphaG_rhs);

  const vec<GF3D2<CCTK_REAL>, 3> gf_dtbeta{
      GF3D2<CCTK_REAL>(layout2, betaGx_rhs),
      GF3D2<CCTK_REAL>(layout2, betaGy_rhs),
      GF3D2<CCTK_REAL>(layout2, betaGz_rhs)};

  // simd types
  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  // parameters
  const CCTK_REAL cpi = acos(-1.0);
  const CCTK_REAL ckappa1 = kappa1;
  const CCTK_REAL ckappa2 = kappa2;
  const CCTK_REAL cmuL = f_mu_L;
  const CCTK_REAL cmuS = f_mu_S;
  const CCTK_REAL ceta = eta;

  // Loop
#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4co_RHS::rhs");
#endif

#include "../wolfram/Z4co_set_rhs.hxx"

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif

  // Upwind and dissipation terms

  // TODO: Consider fusing the loops to reduce memory bandwidth

  apply_upwind_diss(cctkGH, gf_chi, gf_beta, gf_dtchi);

  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      apply_upwind_diss(cctkGH, gf_gamt(a, b), gf_beta, gf_dtgamt(a, b));

  apply_upwind_diss(cctkGH, gf_exKh, gf_beta, gf_dtexKh);

  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      apply_upwind_diss(cctkGH, gf_exAt(a, b), gf_beta, gf_dtexAt(a, b));

  for (int a = 0; a < 3; ++a)
    apply_upwind_diss(cctkGH, gf_trGt(a), gf_beta, gf_dttrGt(a));

  if (!set_Theta_zero)
    apply_upwind_diss(cctkGH, gf_Theta, gf_beta, gf_dtTheta);

  apply_upwind_diss(cctkGH, gf_alpha, gf_beta, gf_dtalpha);

  for (int a = 0; a < 3; ++a)
    apply_upwind_diss(cctkGH, gf_beta(a), gf_beta, gf_dtbeta(a));
}

} // namespace Z4co
