#ifndef MULTIPOLE_SPHERICALHARMONIC_HXX
#define MULTIPOLE_SPHERICALHARMONIC_HXX

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#include <vector>

namespace Multipole {

void SphericalHarmonic(int s, int l, int m, CCTK_REAL th, CCTK_REAL ph,
                       CCTK_REAL *reY, CCTK_REAL *imY);

void HarmonicSetup(int s, int l, int m, int array_size,
                   const std::vector<CCTK_REAL> th,
                   const std::vector<CCTK_REAL> ph, std::vector<CCTK_REAL> &reY,
                   std::vector<CCTK_REAL> &imY);

} // namespace Multipole

#endif // #ifndef MULTIPOLE_SPHERICALHARMONIC_HXX
