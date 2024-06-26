/* surface.cxx */
/* (c) Liwei Ji 06/2024 */

#include "surface.hxx"

#include <loop_device.hxx>

#include <util_Table.h>

namespace Multipole {

void Surface::interpolate(CCTK_ARGUMENTS, int realFieldIndex,
                          int imagFieldIndex) {
  DECLARE_CCTK_PARAMETERS;

  const CCTK_INT nPoints =
      CCTK_MyProc(cctkGH) == 0 ? (nTheta_ + 1) * (nPhi_ + 1) : 0;

  const void *interpCoords[Loop::dim] = {x_.data(), y_.data(), z_.data()};

  CCTK_INT nInputArrays = imagFieldIndex == -1 ? 1 : 2;
  CCTK_INT nOutputArrays = imagFieldIndex == -1 ? 1 : 2;

  const CCTK_INT inputArrayIndices[2] = {realFieldIndex, imagFieldIndex};

  // Interpolation result
  CCTK_POINTER outputArrays[2] = {realF_.data(), imagF_.data()};

  /* DriverInterpolate arguments that aren't currently used */
  const int coordSystemHandle = 0;
  const CCTK_INT interpCoordsTypeCode = 0;
  const CCTK_INT outputArrayTypes[1] = {0};

  int interpHandle = CCTK_InterpHandle("CarpetX");
  if (interpHandle < 0) {
    CCTK_VERROR("Could not obtain interpolator handle for built-in 'CarpetX' "
                "interpolator: %d",
                interpHandle);
  }

  // Interpolation parameter table
  int paramTableHandle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);

  int ierr = Util_TableSetFromString(paramTableHandle, interpolator_pars);

  ierr = DriverInterpolate(cctkGH, Loop::dim, interpHandle, paramTableHandle,
                           coordSystemHandle, nPoints, interpCoordsTypeCode,
                           interpCoords, nInputArrays, inputArrayIndices,
                           nOutputArrays, outputArrayTypes, outputArrays);

  if (ierr < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "CCTK_InterpGridArrays returned error code %d", ierr);
  }

  if (imagFieldIndex == -1) {
    std::fill(imagF_.begin(), imagF_.end(), 0);
  }

  Util_TableDestroy(paramTableHandle);
}

// Take the integral of conj(array1)*array2*sin(th)
void Surface::integrate(const std::vector<CCTK_REAL> &array1r,
                        const std::vector<CCTK_REAL> &array1i,
                        const std::vector<CCTK_REAL> &array2r,
                        const std::vector<CCTK_REAL> &array2i, CCTK_REAL *outRe,
                        CCTK_REAL *outIm) {
  DECLARE_CCTK_PARAMETERS;

  std::vector<CCTK_REAL> fReal(theta_.size());
  std::vector<CCTK_REAL> fImag(theta_.size());

  // integrand: conj(array1)*array2*sin(th)
  for (size_t i = 0; i < theta_.size(); ++i) {
    fReal[i] = (array1r[i] * array2r[i] + array1i[i] * array2i[i]) *
               std::sin(theta_[i]);
    fImag[i] = (array1r[i] * array2i[i] - array1i[i] * array2r[i]) *
               std::sin(theta_[i]);
  }

  if (CCTK_Equals(integration_method, "midpoint")) {
    *outRe = Midpoint2DIntegral(fReal.data(), nTheta_, nPhi_, dTheta_, dPhi_);
    *outIm = Midpoint2DIntegral(fImag.data(), nTheta_, nPhi_, dTheta_, dPhi_);
  } else if (CCTK_Equals(integration_method, "trapezoidal")) {
    *outRe =
        Trapezoidal2DIntegral(fReal.data(), nTheta_, nPhi_, dTheta_, dPhi_);
    *outIm =
        Trapezoidal2DIntegral(fImag.data(), nTheta_, nPhi_, dTheta_, dPhi_);
  } else if (CCTK_Equals(integration_method, "Simpson")) {
    if (nPhi_ % 2 != 0 || nTheta_ % 2 != 0) {
      CCTK_WARN(CCTK_WARN_ABORT, "The Simpson integration method requires even "
                                 "nTheta_ and even nPhi_");
    }
    *outRe = Simpson2DIntegral(fReal.data(), nTheta_, nPhi_, dTheta_, dPhi_);
    *outIm = Simpson2DIntegral(fImag.data(), nTheta_, nPhi_, dTheta_, dPhi_);
  } else if (CCTK_Equals(integration_method, "DriscollHealy")) {
    if (nTheta_ % 2 != 0) {
      CCTK_WARN(CCTK_WARN_ABORT,
                "The Driscoll&Healy integration method requires even nTheta_");
    }
    *outRe =
        DriscollHealy2DIntegral(fReal.data(), nTheta_, nPhi_, dTheta_, dPhi_);
    *outIm =
        DriscollHealy2DIntegral(fImag.data(), nTheta_, nPhi_, dTheta_, dPhi_);
  } else {
    CCTK_WARN(CCTK_WARN_ABORT, "internal error");
  }
}

} // namespace Multipole
