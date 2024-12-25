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
  const CCTK_REAL local_outer_radius = lapse_mask_outer_radius;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const CCTK_REAL rad = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        if (alp(p.I) < local_cutoff && rad > local_outer_radius) {
          HC(p.I) = 0.0;
          MCx(p.I) = 0.0;
          MCy(p.I) = 0.0;
          MCz(p.I) = 0.0;
        }
      });
}

} // namespace ADMconstraints
