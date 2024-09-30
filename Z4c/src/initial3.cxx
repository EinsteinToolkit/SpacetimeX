#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include <cmath>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_Initial3(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Z4c_Initial3;

  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      chi_pre(p.I) = 0.0;

                                      gammatxx_pre(p.I) = 0.0;
                                      gammatxy_pre(p.I) = 0.0;
                                      gammatxz_pre(p.I) = 0.0;
                                      gammatyy_pre(p.I) = 0.0;
                                      gammatyz_pre(p.I) = 0.0;
                                      gammatzz_pre(p.I) = 0.0;

                                      Kh_pre(p.I) = 0.0;

                                      Atxx_pre(p.I) = 0.0;
                                      Atxy_pre(p.I) = 0.0;
                                      Atxz_pre(p.I) = 0.0;
                                      Atyy_pre(p.I) = 0.0;
                                      Atyz_pre(p.I) = 0.0;
                                      Atzz_pre(p.I) = 0.0;

                                      Gamtx_pre(p.I) = 0.0;
                                      Gamty_pre(p.I) = 0.0;
                                      Gamtz_pre(p.I) = 0.0;

                                      Theta_pre(p.I) = 0.0;

                                      alphaG_pre(p.I) = 0.0;

                                      betaGx_pre(p.I) = 0.0;
                                      betaGy_pre(p.I) = 0.0;
                                      betaGz_pre(p.I) = 0.0;
                                    });
}

} // namespace Z4c
