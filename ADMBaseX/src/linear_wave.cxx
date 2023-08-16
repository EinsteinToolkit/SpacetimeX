#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace ADMBaseX {
using namespace Loop;
using namespace std;

extern "C" void ADMBaseX_linear_wave(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMBaseX_linear_wave;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  // See arXiv:1111.2177 [gr-qc], (74-75)

  const auto b = [&](const PointDesc &p) {
    return linear_wave_amplitude *
           sin(2 * CCTK_REAL(M_PI) * (p.x - t) / linear_wave_wavelength);
  };
  const auto bt = [&](const PointDesc &p) {
    return -2 * CCTK_REAL(M_PI) * linear_wave_amplitude /
           linear_wave_wavelength *
           cos(2 * CCTK_REAL(M_PI) * (p.x - t) / linear_wave_wavelength);
  };

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      gxx(p.I) = 1;
                                      gxy(p.I) = 0;
                                      gxz(p.I) = 0;
                                      gyy(p.I) = 1 + b(p);
                                      gyz(p.I) = 0;
                                      gzz(p.I) = 1 - b(p);

                                      kxx(p.I) = 0;
                                      kxy(p.I) = 0;
                                      kxz(p.I) = 0;
                                      kyy(p.I) = bt(p) / 2;
                                      kyz(p.I) = 0;
                                      kzz(p.I) = -bt(p) / 2;
                                    });
}

} // namespace ADMBaseX
