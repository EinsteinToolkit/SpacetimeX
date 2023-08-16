#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

namespace ADMBaseX {
using namespace std;
using namespace Loop;

extern "C" void ADMBaseX_initial_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMBaseX_initial_data;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      gxx(p.I) = 1;
                                      gxy(p.I) = 0;
                                      gxz(p.I) = 0;
                                      gyy(p.I) = 1;
                                      gyz(p.I) = 0;
                                      gzz(p.I) = 1;

                                      kxx(p.I) = 0;
                                      kxy(p.I) = 0;
                                      kxz(p.I) = 0;
                                      kyy(p.I) = 0;
                                      kyz(p.I) = 0;
                                      kzz(p.I) = 0;
                                    });
}

extern "C" void ADMBaseX_initial_lapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMBaseX_initial_lapse;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { alp(p.I) = 1; });
}

extern "C" void ADMBaseX_initial_dtlapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMBaseX_initial_dtlapse;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { dtalp(p.I) = 0; });
}

extern "C" void ADMBaseX_initial_shift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMBaseX_initial_shift;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      betax(p.I) = 0;
                                      betay(p.I) = 0;
                                      betaz(p.I) = 0;
                                    });
}

extern "C" void ADMBaseX_initial_dtshift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMBaseX_initial_dtshift;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      dtbetax(p.I) = 0;
                                      dtbetay(p.I) = 0;
                                      dtbetaz(p.I) = 0;
                                    });
}

} // namespace ADMBaseX
