/* surface.cxx */
/* (c) Liwei Ji 06/2024 */

#include "surface.hxx"

#include <loop_device.hxx>

#include <util_Table.h>

namespace Multipole {

void Surface::interp(CCTK_ARGUMENTS, int realIdx, int imagIdx) {
  DECLARE_CCTK_PARAMETERS;

  const CCTK_INT nPoints =
      CCTK_MyProc(cctkGH) == 0 ? (nTheta_ + 1) * (nPhi_ + 1) : 0;

  const void *interpCoords[Loop::dim] = {(const void *)x_.data(),
                                         (const void *)y_.data(),
                                         (const void *)z_.data()};

  CCTK_INT nInputArrays = imagIdx == -1 ? 1 : 2;
  CCTK_INT nOutputArrays = imagIdx == -1 ? 1 : 2;

  const CCTK_INT inputArrayIndices[2] = {realIdx, imagIdx};

  // Interpolation result
  CCTK_POINTER outputArrays[2];
  outputArrays[0] = real_.data();
  outputArrays[1] = imag_.data();

  /* DriverInterpolate arguments that aren't currently used */
  const int coordSystemHandle = 0;
  CCTK_INT const interpCoords_type_code = 0;
  CCTK_INT const outputArrayTypes[1] = {0};

  int interpHandle = CCTK_InterpHandle("CarpetX");
  if (interpHandle < 0) {
    CCTK_VERROR("Could not obtain inteprolator handle for built-in 'CarpetX' "
                "interpolator: %d",
                interpHandle);
  }

  // Interpolation parameter table
  int paramTableHandle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);

  int ierr = Util_TableSetFromString(paramTableHandle, interpolator_pars);

  ierr = DriverInterpolate(cctkGH, Loop::dim, interpHandle, paramTableHandle,
                           coordSystemHandle, nPoints, interpCoords_type_code,
                           interpCoords, nInputArrays, inputArrayIndices,
                           nOutputArrays, outputArrayTypes, outputArrays);

  if (ierr < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "CCTK_InterpGridArrays returned error code %d", ierr);
  }

  if (imagIdx == -1) {
    for (int i = 0; i < (nTheta_ + 1) * (nPhi_ + 1); i++) {
      imag_[i] = 0;
    }
  }

  Util_TableDestroy(paramTableHandle);
}

} // namespace Multipole
