#ifndef MULTIPOLE_UTILS_HXX
#define MULTIPOLE_UTILS_HXX

#include "multipole.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Multipole {
using namespace std;

enum mp_coord { mp_theta, mp_phi };

static inline int Index_2d(int it, int ip, int ntheta) {
  return it + (ntheta + 1) * ip;
}

void Output1D(CCTK_ARGUMENTS, const string &name, int array_size,
              CCTK_REAL const th[], CCTK_REAL const ph[], mp_coord coord,
              CCTK_REAL const data[]);

void OutputComplexToFile(CCTK_ARGUMENTS, const string &name, CCTK_REAL redata,
                         CCTK_REAL imdata);

void OutputComplexToH5File(CCTK_ARGUMENTS, const variable_desc vars[],
                           const CCTK_REAL radii[], const ModeArray &modes);

void CoordSetup(CCTK_REAL xhat[], CCTK_REAL yhat[], CCTK_REAL zhat[],
                CCTK_REAL th[], CCTK_REAL ph[]);

void ScaleCartesian(int ntheta, int nphi, CCTK_REAL r, CCTK_REAL const xhat[],
                    CCTK_REAL const yhat[], CCTK_REAL const zhat[],
                    CCTK_REAL x[], CCTK_REAL y[], CCTK_REAL z[]);

} // namespace Multipole

#endif // #ifndef MULTIPOLE_UTILS_HXX
