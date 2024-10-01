#include <loop_device.hxx>
#include <simd.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

namespace Z4cow {
using namespace Arith;
using namespace Loop;

extern "C" void Z4cow_InitialLapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Z4cow_InitialLapse;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);
  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);
        // Store
        alphaG.store(mask, index1, W(mask, index1));
      });
}

} // namespace Z4cow
