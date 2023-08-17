#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <random>

namespace ADMBaseX {
using namespace Loop;
using namespace std;

extern "C" void ADMBaseX_add_noise(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMBaseX_add_noise;
  DECLARE_CCTK_PARAMETERS;

  // Hardware random device
  random_device device;
  // Create and seed software random number engine from hardware random number
  default_random_engine engine(device());
  // Random number distribution
  uniform_real_distribution<CCTK_REAL> distribution(-noise_amplitude,
                                                    noise_amplitude);
  const auto add_noise = [&](CCTK_REAL &restrict var) {
    var += distribution(engine);
  };

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      add_noise(gxx(p.I));
                                      add_noise(gxy(p.I));
                                      add_noise(gxz(p.I));
                                      add_noise(gyy(p.I));
                                      add_noise(gyz(p.I));
                                      add_noise(gzz(p.I));

                                      add_noise(kxx(p.I));
                                      add_noise(kxy(p.I));
                                      add_noise(kxz(p.I));
                                      add_noise(kyy(p.I));
                                      add_noise(kyz(p.I));
                                      add_noise(kzz(p.I));

                                      add_noise(alp(p.I));

                                      add_noise(dtalp(p.I));

                                      add_noise(betax(p.I));
                                      add_noise(betay(p.I));
                                      add_noise(betaz(p.I));

                                      add_noise(dtbetax(p.I));
                                      add_noise(dtbetay(p.I));
                                      add_noise(dtbetaz(p.I));
                                    });
}

} // namespace ADMBaseX
