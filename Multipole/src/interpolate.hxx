#ifndef MULTIPOLE_INTERPOLATE_HXX
#define MULTIPOLE_INTERPOLATE_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>
#include <util_Table.h>

namespace Multipole {

// This function interpolates psi4 onto the sphere in cartesian coordinates as
// created by CoordSetup.
void Interp(CCTK_ARGUMENTS, CCTK_REAL x[], CCTK_REAL y[], CCTK_REAL z[],
            int real_idx, int imag_idx, CCTK_REAL psi4r[], CCTK_REAL psi4i[]);

} // namespace Multipole

#endif // #ifndef MULTIPOLE_INTERPOLATE_HXX
