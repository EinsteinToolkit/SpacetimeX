#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef __CUDACC__
// Disable CCTK_DEBUG since the debug information takes too much
// parameter space to launch the kernels
#ifdef CCTK_DEBUG
#undef CCTK_DEBUG
#endif
#endif

#include <derivs.hxx>
#include <loop_device.hxx>
#include <simd.hxx>

#ifdef __CUDACC__
#include <nvtx3/nvToolsExt.h>
#endif

#include <cmath>

namespace ADMconstraints {
using namespace Arith;
using namespace Loop;

template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T Power(T x, int y) {
  return (y == 2) ? Arith::pow2(x) : Arith::pown(x, y);
}

extern "C" void ADMconstraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMconstraints;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2);

  const vect<CCTK_REAL, dim> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};
  vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  // suffix 2: with ghost zones, suffix 5: without ghost zones
  const GF3D2layout layout2(cctkGH, indextype);
  const GF3D5layout layout5(imin, imax);

  // Input grid functions
  const GF3D2<const CCTK_REAL> &gf_W = W;
  const smat<GF3D2<const CCTK_REAL>, 3> gf_gamt{gammatxx, gammatxy, gammatxz,
                                                gammatyy, gammatyz, gammatzz};
  const GF3D2<const CCTK_REAL> &gf_exKh = Kh;
  const smat<GF3D2<const CCTK_REAL>, 3> gf_exAt{Atxx, Atxy, Atxz,
                                                Atyy, Atyz, Atzz};
  const vec<GF3D2<const CCTK_REAL>, 3> gf_trGt{Gamtx, Gamty, Gamtz};
  const GF3D2<const CCTK_REAL> &gf_Theta = Theta;
  const GF3D2<const CCTK_REAL> &gf_alpha = alphaG;
  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{betaGx, betaGy, betaGz};

  // More input grid functions
  const GF3D2<const CCTK_REAL> &gf_eTtt = eTtt;
  const vec<GF3D2<const CCTK_REAL>, 3> gf_eTt{eTtx, eTty, eTtz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_eT{eTxx, eTxy, eTxz,
                                              eTyy, eTyz, eTzz};

  // Output grid functions
  const vec<GF3D2<CCTK_REAL>, 3> gf_ZtC{ZtCx, ZtCy, ZtCz};
  const GF3D2<CCTK_REAL> &gf_HC = HC;
  const vec<GF3D2<CCTK_REAL>, 3> gf_MtC{MtCx, MtCy, MtCz};

  // Define derivs lambdas
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

  // Tile variables for derivatives and so on
  const int ntmps = 118;
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
  const auto make_mat_vec_gf = [&]() { return make_mat(make_vec_gf); };
  const auto make_mat_mat_gf = [&]() { return make_mat(make_mat_gf); };

  const GF3D5<CCTK_REAL> tl_W(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dW(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3> tl_ddW(make_mat_gf());
  calcderivs2(tl_W, tl_dW, tl_ddW, gf_W);

  const smat<GF3D5<CCTK_REAL>, 3> tl_gamt(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> tl_dgamt(make_mat_vec_gf());
  const smat<smat<GF3D5<CCTK_REAL>, 3>, 3> tl_ddgamt(make_mat_mat_gf());
  calcderivs2(tl_gamt, tl_dgamt, tl_ddgamt, gf_gamt);

  const GF3D5<CCTK_REAL> tl_exKh(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dexKh(make_vec_gf());
  calcderivs(tl_exKh, tl_dexKh, gf_exKh);

  const smat<GF3D5<CCTK_REAL>, 3> tl_exAt(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> tl_dexAt(make_mat_vec_gf());
  calcderivs(tl_exAt, tl_dexAt, gf_exAt);

  const vec<GF3D5<CCTK_REAL>, 3> tl_trGt(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> tl_dtrGt(make_vec_vec_gf());
  calcderivs(tl_trGt, tl_dtrGt, gf_trGt);

  const GF3D5<CCTK_REAL> tl_Theta(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dTheta(make_vec_gf());
  calcderivs(tl_Theta, tl_dTheta, gf_Theta);

  const GF3D5<CCTK_REAL> tl_alpha(make_gf());
  calccopy(tl_alpha, gf_alpha);

  const vec<GF3D5<CCTK_REAL>, 3> tl_beta(make_vec_gf());
  calccopy(tl_beta, gf_beta);

  if (itmp != ntmps)
    CCTK_VERROR("Wrong number of temporary variables: ntmps=%d itmp=%d", ntmps,
                itmp);
  itmp = -1;

  // simd types
  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  // parameters
  const CCTK_REAL cpi = M_PI;

#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4cow_Constraints::constraints");
#endif

#include "../wolfram/Z4cow_set_constraint.hxx"

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace ADMconstraints
