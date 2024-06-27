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

void Surface::output1DSingle(CCTK_ARGUMENTS, const std::string &fileName,
                             MpCoord coord,
                             const std::vector<CCTK_REAL> &data) const {
  DECLARE_CCTK_ARGUMENTS;

  const int n = (coord == MpTheta) ? nTheta_ : nPhi_;
  const std::vector<CCTK_REAL> &x = (coord == MpTheta) ? theta_ : phi_;

  if (FILE *f = OpenOutputFile(CCTK_PASS_CTOC, fileName)) {
    fprintf(f, "\"Time = %.19g\n", cctk_time);

    for (int i = 0; i <= n; ++i) {
      int idx = (coord == MpTheta) ? index2D(i, 0) : index2D(nTheta_ / 4, i);
      fprintf(f, "%f %.19g\n", x[idx], data[idx]);
    }

    fprintf(f, "\n\n");
    fclose(f);
  }
}

void Surface::output1D(CCTK_ARGUMENTS, const VariableParse &var,
                       CCTK_REAL rad) const {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_MyProc(cctkGH) == 0 && output_tsv) {
    if (out_1d_every != 0 && cctk_iteration % out_1d_every == 0) {
      std::ostringstream realBase;
      realBase << "mp_" << std::string(CCTK_VarName(var.realIndex)) << "_r"
               << std::fixed << std::setprecision(2) << rad;

      // Output real part data
      output1DSingle(CCTK_PASS_CTOC, realBase.str() + ".th.tsv", MpTheta,
                     realF_);
      output1DSingle(CCTK_PASS_CTOC, realBase.str() + ".ph.tsv", MpPhi, realF_);

      // Output imaginary part data if available
      if (var.imagIndex != -1) {
        std::ostringstream imagBase;
        imagBase << "mp_" << std::string(CCTK_VarName(var.imagIndex)) << "_r"
                 << std::fixed << std::setprecision(2) << rad;
        output1DSingle(CCTK_PASS_CTOC, imagBase.str() + ".th.tsv", MpTheta,
                       imagF_);
        output1DSingle(CCTK_PASS_CTOC, imagBase.str() + ".ph.tsv", MpPhi,
                       imagF_);
      }
    }
  }
}

} // namespace Multipole
