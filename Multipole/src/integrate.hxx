#ifndef MULTIPOLE_INTEGRATE_HXX
#define MULTIPOLE_INTEGRATE_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <vector>

namespace Multipole {

CCTK_REAL Midpoint2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                             CCTK_REAL hy);

CCTK_REAL Trapezoidal2DIntegral(CCTK_REAL const *f, int nx, int ny,
                                CCTK_REAL hx, CCTK_REAL hy);

CCTK_REAL Simpson2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                            CCTK_REAL hy);

CCTK_REAL DriscollHealy2DIntegral(CCTK_REAL const *f, int nx, int ny,
                                  CCTK_REAL hx, CCTK_REAL hy);

void Integrate(int array_size, int nthetap,
               const std::vector<CCTK_REAL> array1r,
               const std::vector<CCTK_REAL> array1i,
               const std::vector<CCTK_REAL> array2r,
               const std::vector<CCTK_REAL> array2i, CCTK_REAL const th[],
               CCTK_REAL const ph[], CCTK_REAL *outre, CCTK_REAL *outim);

} // namespace Multipole

#endif // #ifndef MULTIPOLE_INTEGRATE_HXX
