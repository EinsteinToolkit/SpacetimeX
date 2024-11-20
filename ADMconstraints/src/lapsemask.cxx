#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

namespace ADMconstraints {
using namespace Arith;
using namespace Loop;

extern "C" void ADMconstraints_LapseMask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMconstraints_LapseMask;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL local_cutoff = lapse_mask_cutoff;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        if (alp(p.I) < local_cutoff) {
          HC(p.I) = 0.0;
          MCx(p.I) = 0.0;
          MCy(p.I) = 0.0;
          MCz(p.I) = 0.0;
        }
      });
}

} // namespace ADMconstraints
