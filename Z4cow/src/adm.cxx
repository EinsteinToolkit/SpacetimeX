#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef __CUDACC__
#include <nvtx3/nvToolsExt.h>
#endif

#include <cmath>

namespace Z4cow {
using namespace Arith;
using namespace Loop;

template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T Power(T x, int y) {
  return (y == 2) ? Arith::pow2(x) : Arith::pown(x, y);
}

extern "C" void Z4cow_ADM(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Z4cow_ADM;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout2(cctkGH, indextype);

  // Input Z4c grid functions
  const GF3D2<const CCTK_REAL> &gf_W = W;
  const smat<GF3D2<const CCTK_REAL>, 3> gf_gamt{gammatxx, gammatxy, gammatxz,
                                                gammatyy, gammatyz, gammatzz};
  const GF3D2<const CCTK_REAL> &gf_exKh = Kh;
  const smat<GF3D2<const CCTK_REAL>, 3> gf_exAt{Atxx, Atxy, Atxz,
                                                Atyy, Atyz, Atzz};
  const GF3D2<const CCTK_REAL> &gf_Theta = Theta;
  const GF3D2<const CCTK_REAL> &gf_alpha = alphaG;
  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{betaGx, betaGy, betaGz};

  // Output ADM grid functions
  const smat<GF3D2<CCTK_REAL>, 3> gf_ADMgam{gxx, gxy, gxz, gyy, gyz, gzz};
  const smat<GF3D2<CCTK_REAL>, 3> gf_ADMK{kxx, kxy, kxz, kyy, kyz, kzz};
  const GF3D2<CCTK_REAL> &gf_ADMalpha = alp;
  const vec<GF3D2<CCTK_REAL>, 3> gf_ADMbeta{betax, betay, betaz};

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4cow_ADM::adm");
#endif

#include "../wolfram/Z4cow_set_ADM.hxx"

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace Z4cow
