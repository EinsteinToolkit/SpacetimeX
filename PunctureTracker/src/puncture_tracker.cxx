#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>
#include <loop_device.hxx>
#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <array>
#include <ctype.h>

namespace PunctureTracker {

const int max_num_tracked = 10;

extern "C" void PunctureTracker_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PunctureTracker_Init;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("Initializing PunctureTracker");
  }

  for (int n = 0; n < max_num_tracked; ++n) {
    if (track[n]) {
      pt_loc_t[n] = cctk_time;
      pt_loc_x[n] = initial_x[n];
      pt_loc_y[n] = initial_y[n];
      pt_loc_z[n] = initial_z[n];
      pt_vel_t[n] = cctk_time;
      pt_vel_x[n] = 0.0;
      pt_vel_y[n] = 0.0;
      pt_vel_z[n] = 0.0;
    } else {
      // Initialise to some sensible but unimportant values
      pt_loc_t[n] = 0.0;
      pt_loc_x[n] = 0.0;
      pt_loc_y[n] = 0.0;
      pt_loc_z[n] = 0.0;
      pt_vel_t[n] = 0.0;
      pt_vel_x[n] = 0.0;
      pt_vel_y[n] = 0.0;
      pt_vel_z[n] = 0.0;
    }
  }

  // enabled if refinement regions should follow the punctures
  if (track_boxes) {
    const int max_num_regions = 2;
    for (int i = 0; i < max_num_regions; i++) {
      CCTK_VINFO("Writing punc coords to box %d.", i);
      position_x[i] = pt_loc_x[i];
      position_y[i] = pt_loc_y[i];
      position_z[i] = pt_loc_z[i];
    }
  }
}

extern "C" void PunctureTracker_Track(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PunctureTracker_Track;
  DECLARE_CCTK_PARAMETERS;

  // Do not track while setting up initial data;
  // time interpolation may fail

  if (cctk_iteration == 0) {
    return;
  }

  // Some output

  if (verbose) {
    CCTK_INFO("Tracking punctures...");
  }

  if (verbose) {
    for (int n = 0; n < max_num_tracked; ++n) {
      if (track[n]) {
        CCTK_VINFO("Puncture #%d is at (%g,%g,%g)", n, double(pt_loc_x[n]),
                   double(pt_loc_y[n]), double(pt_loc_z[n]));
      }
    }
  }

  // Manual time level cycling
  CCTK_REAL pt_t_prev[max_num_tracked];

  for (int n = 0; n < max_num_tracked; ++n) {
    if (track[n]) {
      pt_t_prev[n] = pt_loc_t[n];
      pt_loc_t[n] = cctk_time;
      pt_vel_t[n] = cctk_time;
    } else {
      pt_t_prev[n] = 0.0;
    }
  }

  // Interpolate
  {
    // Only processor 0 interpolates
    const CCTK_INT nPoints = CCTK_MyProc(cctkGH) == 0 ? max_num_tracked : 0;

    // Interpolation coordinates
    const void *interpCoords[Loop::dim] = {pt_loc_x, pt_loc_y, pt_loc_z};

    // Interpolated variables
    const CCTK_INT nInputArrays = 3;
    const CCTK_INT inputArrayIndices[3] = {CCTK_VarIndex("ADMBaseX::betax"),
                                           CCTK_VarIndex("ADMBaseX::betay"),
                                           CCTK_VarIndex("ADMBaseX::betaz")};

    // Interpolation result
    CCTK_REAL pt_betax[max_num_tracked];
    CCTK_REAL pt_betay[max_num_tracked];
    CCTK_REAL pt_betaz[max_num_tracked];
    CCTK_POINTER outputArrays[3] = {pt_betax, pt_betay, pt_betaz};

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

    if ((ierr = Util_TableSetInt(paramTableHandle, interp_order, "order")) <
        0) {
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

    if (CCTK_MyProc(cctkGH) == 0) {

      // Some more output

      if (verbose) {
        for (int n = 0; n < max_num_tracked; ++n) {
          if (track[n]) {
            CCTK_VINFO("Shift at puncture #%d is at (%g,%g,%g)", n,
                       double(pt_betax[n]), double(pt_betay[n]),
                       double(pt_betaz[n]));
          }
        }
      }

      // Check for NaNs and large shift components
      for (int n = 0; n < max_num_tracked; ++n) {
        if (track[n]) {
          CCTK_REAL norm = sqrt(pow(pt_betax[n], 2) + pow(pt_betay[n], 2) +
                                pow(pt_betaz[n], 2));

          if (!CCTK_isfinite(norm) || norm > shift_limit) {
            CCTK_VERROR("Shift at puncture #%d is (%g,%g,%g).  This likely "
                        "indicates an error in the simulation.",
                        n, double(pt_betax[n]), double(pt_betay[n]),
                        double(pt_betaz[n]));
          }
        }
      }

      // Time evolution
      for (int n = 0; n < max_num_tracked; ++n) {
        if (track[n]) {
          const CCTK_REAL dt = pt_loc_t[n] - pt_t_prev[n];
          // First order time integrator
          // Michael Koppitz says this works...
          // if it doesn't, we can make it second order accurate
          pt_loc_x[n] += dt * (-pt_betax[n]);
          pt_loc_y[n] += dt * (-pt_betay[n]);
          pt_loc_z[n] += dt * (-pt_betaz[n]);
          pt_vel_x[n] = -pt_betax[n];
          pt_vel_y[n] = -pt_betay[n];
          pt_vel_z[n] = -pt_betaz[n];
        }
      }
    }

    // Broadcast result: 3 components for location, 3 components for velocity
    CCTK_REAL loc_global[6 * max_num_tracked];
    if (CCTK_MyProc(cctkGH) == 0) {
      for (int n = 0; n < max_num_tracked; ++n) {
        loc_global[n] = pt_loc_x[n];
        loc_global[max_num_tracked + n] = pt_loc_y[n];
        loc_global[2 * max_num_tracked + n] = pt_loc_z[n];
        loc_global[3 * max_num_tracked + n] = pt_vel_x[n];
        loc_global[4 * max_num_tracked + n] = pt_vel_y[n];
        loc_global[5 * max_num_tracked + n] = pt_vel_z[n];
      }
    }

    // CarpetX doesn't register reduction handle, here's a quick fix.
    MPI_Bcast(loc_global, 6 * max_num_tracked, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int n = 0; n < max_num_tracked; ++n) {
      pt_loc_x[n] = loc_global[n];
      pt_loc_y[n] = loc_global[max_num_tracked + n];
      pt_loc_z[n] = loc_global[2 * max_num_tracked + n];
      pt_vel_x[n] = loc_global[3 * max_num_tracked + n];
      pt_vel_y[n] = loc_global[4 * max_num_tracked + n];
      pt_vel_z[n] = loc_global[5 * max_num_tracked + n];
    }
  }

  if (track_boxes) {
    const int max_num_regions = 2;
    for (int i = 0; i < max_num_regions; i++) {
      position_x[i] = pt_loc_x[i];
      position_y[i] = pt_loc_y[i];
      position_z[i] = pt_loc_z[i];
    }
  }
}

extern "C" void PunctureTracker_CheckShift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_PunctureTracker_CheckShift;
  DECLARE_CCTK_PARAMETERS;

  const int level = ilogb(CCTK_REAL(cctk_levfac[0]));
  const int finest_lvl = num_levels[0] - 1;

  if (level == finest_lvl) {
    for (int n = 0; n < max_num_tracked; ++n) {
      if (track[n]) {
        const Arith::vect<CCTK_REAL, Loop::dim> loc_vec = {
            pt_loc_x[n], pt_loc_y[n], pt_loc_z[n]};

        grid.loop_all_device<0, 0, 0>(
            grid.nghostzones,
            [=] CCTK_DEVICE(
                const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
              if (maximum(abs(p.X - loc_vec)) <= (1 / pow(2, level))) {
                printf("Shift at level %d near puncture #%d is {%g, %g, %g} at "
                       "coords {%g, %g, %g}.",
                       level, n, betax(p.I), betay(p.I), betaz(p.I), p.x, p.y,
                       p.z);
              }
            });
      }
    }
  }
}

} // namespace PunctureTracker
