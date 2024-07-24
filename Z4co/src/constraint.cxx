#include <cctk.h>

#ifdef __CUDACC__
// Disable CCTK_DEBUG since the debug information takes too much
// parameter space to launch the kernels
#ifdef CCTK_DEBUG
#undef CCTK_DEBUG
#endif
#endif

#include "derivs.hxx"

#include <loop_device.hxx>
#include <simd.hxx>

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

extern "C" void Z4co_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4co_Constraints;
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

  const GF3D5<CCTK_REAL> tl_chi(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dchi(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3> tl_ddchi(make_mat_gf());
  calc_derivs2(cctkGH, gf_chi, tl_chi, tl_dchi, tl_ddchi, layout5);

  const smat<GF3D5<CCTK_REAL>, 3> tl_gamt(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> tl_dgamt(make_mat_vec_gf());
  const smat<smat<GF3D5<CCTK_REAL>, 3>, 3> tl_ddgamt(make_mat_mat_gf());
  calc_derivs2(cctkGH, gf_gamt, tl_gamt, tl_dgamt, tl_ddgamt, layout5);

  const GF3D5<CCTK_REAL> tl_exKh(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dexKh(make_vec_gf());
  calc_derivs(cctkGH, gf_exKh, tl_exKh, tl_dexKh, layout5);

  const smat<GF3D5<CCTK_REAL>, 3> tl_exAt(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> tl_dexAt(make_mat_vec_gf());
  calc_derivs(cctkGH, gf_exAt, tl_exAt, tl_dexAt, layout5);

  const vec<GF3D5<CCTK_REAL>, 3> tl_trGt(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> tl_dtrGt(make_vec_vec_gf());
  calc_derivs(cctkGH, gf_trGt, tl_trGt, tl_dtrGt, layout5);

  const GF3D5<CCTK_REAL> tl_Theta(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> tl_dTheta(make_vec_gf());
  calc_derivs(cctkGH, gf_Theta, tl_Theta, tl_dTheta, layout5);

  const GF3D5<CCTK_REAL> tl_alpha(make_gf());
  calc_copy(cctkGH, gf_alpha, tl_alpha, layout5);

  const vec<GF3D5<CCTK_REAL>, 3> tl_beta(make_vec_gf());
  calc_copy(cctkGH, gf_beta, tl_beta, layout5);

  if (itmp != ntmps)
    CCTK_VERROR("Wrong number of temporary variables: ntmps=%d itmp=%d", ntmps,
                itmp);
  itmp = -1;

  //

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

  //

  const vec<GF3D2<CCTK_REAL>, 3> gf_ZtC{GF3D2<CCTK_REAL>(layout2, ZtCx),
                                        GF3D2<CCTK_REAL>(layout2, ZtCy),
                                        GF3D2<CCTK_REAL>(layout2, ZtCz)};

  const GF3D2<CCTK_REAL> gf_HC(layout2, HC);

  const vec<GF3D2<CCTK_REAL>, 3> gf_MtC{GF3D2<CCTK_REAL>(layout2, MtCx),
                                        GF3D2<CCTK_REAL>(layout2, MtCy),
                                        GF3D2<CCTK_REAL>(layout2, MtCz)};

  //

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const CCTK_REAL cpi = acos(-1.0);

  const Loop::GridDescBaseDevice grid(cctkGH);
#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4co_Constraints::constraints");
#endif

#include "../wolfram/Z4co_set_constraint.hxx"

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace Z4co
