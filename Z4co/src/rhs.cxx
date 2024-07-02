#include <cctk.h>

#ifdef __CUDACC__
// Disable CCTK_DEBUG since the debug information takes too much
// parameter space to launch the kernels
#ifdef CCTK_DEBUG
#undef CCTK_DEBUG
#endif
#endif

#include "derivs.hxx"
#include "physics.hxx"

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

extern "C" void Z4co_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4co_RHS;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2 + 1)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2 + 1);

  //

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};
  vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  // Suffix 1: with ghost zones, suffix 0: without ghost zones
  const GF3D2layout layout2(cctkGH, indextype);
  const GF3D5layout layout5(imin, imax);

  const GF3D2<const CCTK_REAL> gf_chi1(layout2, chi);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_gammat1{
      GF3D2<const CCTK_REAL>(layout2, gammatxx),
      GF3D2<const CCTK_REAL>(layout2, gammatxy),
      GF3D2<const CCTK_REAL>(layout2, gammatxz),
      GF3D2<const CCTK_REAL>(layout2, gammatyy),
      GF3D2<const CCTK_REAL>(layout2, gammatyz),
      GF3D2<const CCTK_REAL>(layout2, gammatzz)};

  const GF3D2<const CCTK_REAL> gf_Kh1(layout2, Kh);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_At1{
      GF3D2<const CCTK_REAL>(layout2, Atxx),
      GF3D2<const CCTK_REAL>(layout2, Atxy),
      GF3D2<const CCTK_REAL>(layout2, Atxz),
      GF3D2<const CCTK_REAL>(layout2, Atyy),
      GF3D2<const CCTK_REAL>(layout2, Atyz),
      GF3D2<const CCTK_REAL>(layout2, Atzz)};

  const vec<GF3D2<const CCTK_REAL>, 3> gf_Gamt1{
      GF3D2<const CCTK_REAL>(layout2, Gamtx),
      GF3D2<const CCTK_REAL>(layout2, Gamty),
      GF3D2<const CCTK_REAL>(layout2, Gamtz)};

  const GF3D2<const CCTK_REAL> gf_Theta1(layout2, Theta);

  const GF3D2<const CCTK_REAL> gf_alphaG1(layout2, alphaG);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_betaG1{
      GF3D2<const CCTK_REAL>(layout2, betaGx),
      GF3D2<const CCTK_REAL>(layout2, betaGy),
      GF3D2<const CCTK_REAL>(layout2, betaGz)};

  //

  // Ideas:
  //
  // - Outline certain functions, e.g. `det` or `raise_index`. Ensure
  //   they are called with floating-point arguments, not tensor
  //   indices.

  const int ntmps = 154;
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

  const GF3D5<CCTK_REAL> gf_chi(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_dchi(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3> gf_ddchi(make_mat_gf());
  calc_derivs2(cctkGH, gf_chi1, gf_chi, gf_dchi, gf_ddchi, layout5);

  const smat<GF3D5<CCTK_REAL>, 3> gf_gamt(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dgamt(make_mat_vec_gf());
  const smat<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddgamt(make_mat_mat_gf());
  calc_derivs2(cctkGH, gf_gammat1, gf_gamt, gf_dgamt, gf_ddgamt,
               layout5);

  const GF3D5<CCTK_REAL> gf_exKh(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_dexKh(make_vec_gf());
  calc_derivs(cctkGH, gf_Kh1, gf_exKh, gf_dexKh, layout5);

  const smat<GF3D5<CCTK_REAL>, 3> gf_exAt(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dexAt(make_mat_vec_gf());
  calc_derivs(cctkGH, gf_At1, gf_exAt, gf_dexAt, layout5);

  const vec<GF3D5<CCTK_REAL>, 3> gf_trGt(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dtrGt(make_vec_vec_gf());
  calc_derivs(cctkGH, gf_Gamt1, gf_trGt, gf_dtrGt, layout5);

  const GF3D5<CCTK_REAL> gf_Theta(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_dTheta0(make_vec_gf());
  calc_derivs(cctkGH, gf_Theta1, gf_Theta, gf_dTheta0, layout5);

  const GF3D5<CCTK_REAL> gf_alpha(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_dalpha(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3> gf_ddalpha(make_mat_gf());
  calc_derivs2(cctkGH, gf_alphaG1, gf_alpha, gf_dalpha, gf_ddalpha,
               layout5);

  const vec<GF3D5<CCTK_REAL>, 3> gf_beta(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dbeta(make_vec_vec_gf());
  const vec<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddbeta(make_vec_mat_gf());
  calc_derivs2(cctkGH, gf_betaG1, gf_beta, gf_dbeta, gf_ddbeta, layout5);

  if (itmp != ntmps)
    CCTK_VERROR("Wrong number of temporary variables: ntmps=%d itmp=%d", ntmps,
                itmp);
  itmp = -1;

  //

  const GF3D2<const CCTK_REAL> gf_eTtt1(layout2, eTtt);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_eTti1{
      GF3D2<const CCTK_REAL>(layout2, eTtx),
      GF3D2<const CCTK_REAL>(layout2, eTty),
      GF3D2<const CCTK_REAL>(layout2, eTtz)};

  const smat<GF3D2<const CCTK_REAL>, 3> gf_eTij1{
      GF3D2<const CCTK_REAL>(layout2, eTxx),
      GF3D2<const CCTK_REAL>(layout2, eTxy),
      GF3D2<const CCTK_REAL>(layout2, eTxz),
      GF3D2<const CCTK_REAL>(layout2, eTyy),
      GF3D2<const CCTK_REAL>(layout2, eTyz),
      GF3D2<const CCTK_REAL>(layout2, eTzz)};

  //

  const GF3D2<CCTK_REAL> gf_chi_rhs1(layout2, chi_rhs);

  const smat<GF3D2<CCTK_REAL>, 3> gf_gammat_rhs1{
      GF3D2<CCTK_REAL>(layout2, gammatxx_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatxy_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatxz_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatyy_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatyz_rhs),
      GF3D2<CCTK_REAL>(layout2, gammatzz_rhs)};

  const GF3D2<CCTK_REAL> gf_Kh_rhs1(layout2, Kh_rhs);

  const smat<GF3D2<CCTK_REAL>, 3> gf_At_rhs1{
      GF3D2<CCTK_REAL>(layout2, Atxx_rhs), GF3D2<CCTK_REAL>(layout2, Atxy_rhs),
      GF3D2<CCTK_REAL>(layout2, Atxz_rhs), GF3D2<CCTK_REAL>(layout2, Atyy_rhs),
      GF3D2<CCTK_REAL>(layout2, Atyz_rhs), GF3D2<CCTK_REAL>(layout2, Atzz_rhs)};

  const vec<GF3D2<CCTK_REAL>, 3> gf_Gamt_rhs1{
      GF3D2<CCTK_REAL>(layout2, Gamtx_rhs),
      GF3D2<CCTK_REAL>(layout2, Gamty_rhs),
      GF3D2<CCTK_REAL>(layout2, Gamtz_rhs)};

  const GF3D2<CCTK_REAL> gf_Theta_rhs1(layout2, Theta_rhs);

  const GF3D2<CCTK_REAL> gf_alphaG_rhs1(layout2, alphaG_rhs);

  const vec<GF3D2<CCTK_REAL>, 3> gf_betaG_rhs1{
      GF3D2<CCTK_REAL>(layout2, betaGx_rhs),
      GF3D2<CCTK_REAL>(layout2, betaGy_rhs),
      GF3D2<CCTK_REAL>(layout2, betaGz_rhs)};

  //

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);

#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4co_RHS::rhs");
#endif
  noinline([&]() __attribute__((__flatten__, __hot__)) {
    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          const GF3D2index index2(layout2, p.I);
          const GF3D5index index5(layout5, p.I);

#include "../wolfram/Z4co_set_rhs.hxx"
        });
  });
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif

  // Upwind and dissipation terms

  // TODO: Consider fusing the loops to reduce memory bandwidth

  apply_upwind_diss(cctkGH, gf_chi1, gf_betaG1, gf_chi_rhs1);

  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      apply_upwind_diss(cctkGH, gf_gammat1(a, b), gf_betaG1,
                        gf_gammat_rhs1(a, b));

  apply_upwind_diss(cctkGH, gf_Kh1, gf_betaG1, gf_Kh_rhs1);

  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      apply_upwind_diss(cctkGH, gf_At1(a, b), gf_betaG1, gf_At_rhs1(a, b));

  for (int a = 0; a < 3; ++a)
    apply_upwind_diss(cctkGH, gf_Gamt1(a), gf_betaG1, gf_Gamt_rhs1(a));

  if (!set_Theta_zero)
    apply_upwind_diss(cctkGH, gf_Theta1, gf_betaG1, gf_Theta_rhs1);

  apply_upwind_diss(cctkGH, gf_alphaG1, gf_betaG1, gf_alphaG_rhs1);

  for (int a = 0; a < 3; ++a)
    apply_upwind_diss(cctkGH, gf_betaG1(a), gf_betaG1, gf_betaG_rhs1(a));
}

} // namespace Z4co
