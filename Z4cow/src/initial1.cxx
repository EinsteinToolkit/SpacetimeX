#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#ifdef __CUDACC__
#include <nvToolsExt.h>
#endif

#include <cmath>

namespace Z4cow {
using namespace Arith;
using namespace Loop;

extern "C" void Z4cow_Initial1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Z4cow_Initial1;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  // Input grid functions
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g1{gxx, gxy, gxz, gyy, gyz, gzz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_K1{kxx, kxy, kxz, kyy, kyz, kzz};
  const GF3D2<const CCTK_REAL> &gf_alp1 = alp;
  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta1{betax, betay, betaz};

  // Output grid functions
  const GF3D2<CCTK_REAL> &gf_W1 = W;
  const smat<GF3D2<CCTK_REAL>, 3> gf_gammat1{gammatxx, gammatxy, gammatxz,
                                             gammatyy, gammatyz, gammatzz};
  const GF3D2<CCTK_REAL> &gf_Kh1 = Kh;
  const smat<GF3D2<CCTK_REAL>, 3> gf_At1{Atxx, Atxy, Atxz, Atyy, Atyz, Atzz};
  const GF3D2<CCTK_REAL> &gf_Theta1 = Theta;
  const GF3D2<CCTK_REAL> &gf_alphaG1 = alphaG;
  const vec<GF3D2<CCTK_REAL>, 3> gf_betaG1{betaGx, betaGy, betaGz};

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4cow_Initial1::initial1");
#endif
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);

        // Load (initilize those masked grid with more reasonable value)
        const smat<vreal, 3> g = gf_g1(mask, index1, one<smat<int, 3>>()());
        const smat<vreal, 3> K = gf_K1(mask, index1);
        const vreal alp = gf_alp1(mask, index1, 1);
        const vec<vreal, 3> beta = gf_beta1(mask, index1);

        // Calculate Z4c variables (all except Gammat)
        const vreal detg = calc_det(g);
        const smat<vreal, 3> gu = calc_inv(g, detg);

        const vreal W = 1 / cbrt(sqrt(detg));
        const smat<vreal, 3> gammat(
            [&](int a, int b) ARITH_INLINE { return W * W * g(a, b); });
        const vreal trK = sum_symm<3>(
            [&](int x, int y) ARITH_INLINE { return gu(x, y) * K(x, y); });
        const vreal Theta = 0;
        const vreal Kh = trK - 2 * Theta;
        const smat<vreal, 3> At([&](int a, int b) ARITH_INLINE {
          return W * W * (K(a, b) - trK / 3 * g(a, b));
        });
        const vreal alphaG = alp;
        const vec<vreal, 3> betaG([&](int a) ARITH_INLINE { return beta(a); });

        // Store
        gf_W1.store(mask, index1, W);
        gf_gammat1.store(mask, index1, gammat);
        gf_Kh1.store(mask, index1, Kh);
        gf_At1.store(mask, index1, At);
        gf_Theta1.store(mask, index1, Theta);
        gf_alphaG1.store(mask, index1, alphaG);
        gf_betaG1.store(mask, index1, betaG);
      });
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace Z4cow
