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

namespace Z4coo {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4coo_ADM(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4coo_ADM;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout2(cctkGH, indextype);

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

  const GF3D2<const CCTK_REAL> gf_Theta(layout2, Theta);

  const GF3D2<const CCTK_REAL> gf_alpha(layout2, alphaG);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{
      GF3D2<const CCTK_REAL>(layout2, betaGx),
      GF3D2<const CCTK_REAL>(layout2, betaGy),
      GF3D2<const CCTK_REAL>(layout2, betaGz)};

  const smat<GF3D2<CCTK_REAL>, 3> gf_ADMgam{
      GF3D2<CCTK_REAL>(layout2, gxx), GF3D2<CCTK_REAL>(layout2, gxy),
      GF3D2<CCTK_REAL>(layout2, gxz), GF3D2<CCTK_REAL>(layout2, gyy),
      GF3D2<CCTK_REAL>(layout2, gyz), GF3D2<CCTK_REAL>(layout2, gzz)};

  const smat<GF3D2<CCTK_REAL>, 3> gf_ADMK{
      GF3D2<CCTK_REAL>(layout2, kxx), GF3D2<CCTK_REAL>(layout2, kxy),
      GF3D2<CCTK_REAL>(layout2, kxz), GF3D2<CCTK_REAL>(layout2, kyy),
      GF3D2<CCTK_REAL>(layout2, kyz), GF3D2<CCTK_REAL>(layout2, kzz)};

  const GF3D2<CCTK_REAL> gf_ADMalpha(layout2, alp);

  const vec<GF3D2<CCTK_REAL>, 3> gf_ADMbeta{GF3D2<CCTK_REAL>(layout2, betax),
                                            GF3D2<CCTK_REAL>(layout2, betay),
                                            GF3D2<CCTK_REAL>(layout2, betaz)};

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);

#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4coo_ADM::adm");
#endif

#include "../wolfram/Z4coo_set_ADM.hxx"

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace Z4coo
