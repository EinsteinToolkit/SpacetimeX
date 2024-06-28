#include "puncture.hxx"

#include <util_Table.h>

namespace PunctureTracker {

void PunctureContainer::updatePreviousTime(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  for (size_t n = 0; n < time_.size(); ++n) {
    previousTime_[n] = time_[n];
    time_[n] = cctk_time;
  }
}

void PunctureContainer::interpolate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  // Only processor 0 interpolates
  const CCTK_INT nPoints = CCTK_MyProc(cctkGH) == 0 ? numPunctures_ : 0;

  // Interpolation coordinates
  const void *interpCoords[Loop::dim] = {
      location_[0].data(), location_[1].data(), location_[2].data()};

  // Interpolated variables
  const CCTK_INT nInputArrays = 3;
  const CCTK_INT inputArrayIndices[3] = {CCTK_VarIndex("ADMBaseX::betax"),
                                         CCTK_VarIndex("ADMBaseX::betay"),
                                         CCTK_VarIndex("ADMBaseX::betaz")};

  CCTK_POINTER outputArrays[3] = {beta_[0].data(), beta_[1].data(),
                                  beta_[2].data()};

  /* DriverInterpolate arguments that aren't currently used */
  const int coordSystemHandle = 0;
  const CCTK_INT interpCoordsTypeCode = 0;
  const CCTK_INT outputArrayTypes[1] = {0};

  const int interpHandle = CCTK_InterpHandle("CarpetX");
  if (interpHandle < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Can't get interpolation handle");
    return;
  }

  int ierr;

  int paramTableHandle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (paramTableHandle < 0) {
    CCTK_VERROR("Can't create parameter table: %d", paramTableHandle);
  }

  if ((ierr = Util_TableSetInt(paramTableHandle, interp_order, "order")) < 0) {
    CCTK_VERROR("Can't set order in parameter table: %d", ierr);
  }

  // Interpolate
  ierr = DriverInterpolate(cctkGH, Loop::dim, interpHandle, paramTableHandle,
                           coordSystemHandle, nPoints, interpCoordsTypeCode,
                           interpCoords, nInputArrays, inputArrayIndices,
                           nInputArrays, outputArrayTypes, outputArrays);

  if (ierr < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Interpolation error");
  }

  Util_TableDestroy(paramTableHandle);
}

void PunctureContainer::evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  if (CCTK_MyProc(cctkGH) == 0) {
    // First order time integrator
    // Michael Koppitz says this works...
    // if it doesn't, we can make it second order accurate
    for (size_t n = 0; n < time_.size(); ++n) {
      const CCTK_REAL dt = time_[n] - previousTime_[n];
      for (int i = 0; i < Loop::dim; ++i) {
        location_[i][n] += dt * (-beta_[i][n]);
        velocity_[i][n] = -beta_[i][n];
      }
    }
  }
}

void PunctureContainer::broadcast(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  // 3 components for location, 3 components for velocity
  const CCTK_INT bufferSize = 6 * numPunctures_;

  CCTK_REAL buffer[bufferSize];

  if (CCTK_MyProc(cctkGH) == 0) {
    for (int i = 0; i < Loop::dim; ++i) {
      for (int n = 0; n < numPunctures_; ++n) {
        buffer[i * numPunctures_ + n] = location_[i][n];
        buffer[(i + Loop::dim) * numPunctures_ + n] = velocity_[i][n];
      }
    }
  }

  MPI_Bcast(buffer, bufferSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (int i = 0; i < Loop::dim; ++i) {
    for (int n = 0; n < numPunctures_; ++n) {
      location_[i][n] = buffer[i * numPunctures_ + n];
      velocity_[i][n] = buffer[(i + Loop::dim) * numPunctures_ + n];
    }
  }
}

} // namespace PunctureTracker
