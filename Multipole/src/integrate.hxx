#ifndef MULTIPOLE_INTEGRATE_HXX
#define MULTIPOLE_INTEGRATE_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Multipole {

CCTK_REAL Midpoint2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                             CCTK_REAL hy);

CCTK_REAL Trapezoidal2DIntegral(CCTK_REAL const *f, int nx, int ny,
                                CCTK_REAL hx, CCTK_REAL hy);

CCTK_REAL Simpson2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                            CCTK_REAL hy);

CCTK_REAL DriscollHealy2DIntegral(CCTK_REAL const *f, int nx, int ny,
                                  CCTK_REAL hx, CCTK_REAL hy);

void Integrate(int array_size, int ntheta, CCTK_REAL const array1r[],
               CCTK_REAL const array1i[], CCTK_REAL const array2r[],
               CCTK_REAL const array2i[], CCTK_REAL const th[],
               CCTK_REAL const pph[], CCTK_REAL out_arrayr[],
               CCTK_REAL out_arrayi[]);

} // namespace Multipole

#endif // #ifndef MULTIPOLE_INTEGRATE_HXX
