#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <cassert>
#include <cmath>
#include <cstdio>

namespace PunctureTracker {
using namespace std;

const int max_num_tracked = 10;

extern "C" void PunctureTracker_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
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
    pt_loc_t_p[n] = 0.0;
    pt_loc_x_p[n] = 0.0;
    pt_loc_y_p[n] = 0.0;
    pt_loc_z_p[n] = 0.0;
  }
}

extern "C" void PunctureTracker_Track(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
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

  for (int n = 0; n < max_num_tracked; ++n) {
    if (track[n]) {
      pt_loc_t_p[n] = pt_loc_t[n];
      pt_loc_x_p[n] = pt_loc_x[n];
      pt_loc_y_p[n] = pt_loc_y[n];
      pt_loc_z_p[n] = pt_loc_z[n];

      pt_loc_t[n] = cctk_time;
      pt_vel_t[n] = cctk_time;
    }
  }

  // Interpolate

  // Dimensions
  const int dim = 3;

  // Interpolation operator
  const int operator_handle =
      CCTK_InterpHandle("Lagrange polynomial interpolation");
  if (operator_handle < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Can't get interpolation handle");
    return;
  }

  // Interpolation parameter table
  const int order = 4;
  const int param_table_handle = Util_TableCreateFromString("order=4");
  if (param_table_handle < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Can't create parameter table");
    return;
  }

  {

    // Interpolation coordinate system
    const int coordsys_handle = CCTK_CoordSystemHandle("cart3d");
    if (coordsys_handle < 0) {
      CCTK_WARN(CCTK_WARN_ALERT, "Can't get coordinate system handle");
      goto label_free_param_table;
    }

    // Only processor 0 interpolates
    const int num_points = CCTK_MyProc(cctkGH) == 0 ? max_num_tracked : 0;

    // Interpolation coordinates
    assert(dim == 3);
    CCTK_POINTER_TO_CONST interp_coords[3];
    interp_coords[0] = pt_loc_x_p;
    interp_coords[1] = pt_loc_y_p;
    interp_coords[2] = pt_loc_z_p;

    // Number of interpolation variables
    int const num_vars = 3;

    // Interpolated variables
    assert(num_vars == 3);
    int input_array_indices[3];
    input_array_indices[0] = CCTK_VarIndex("ADMBase::betax");
    input_array_indices[1] = CCTK_VarIndex("ADMBase::betay");
    input_array_indices[2] = CCTK_VarIndex("ADMBase::betaz");

    // Interpolation result types
    assert(num_vars == 3);
    CCTK_INT output_array_type_codes[3];
    output_array_type_codes[0] = CCTK_VARIABLE_REAL;
    output_array_type_codes[1] = CCTK_VARIABLE_REAL;
    output_array_type_codes[2] = CCTK_VARIABLE_REAL;

    // Interpolation result
    CCTK_REAL pt_betax[max_num_tracked];
    CCTK_REAL pt_betay[max_num_tracked];
    CCTK_REAL pt_betaz[max_num_tracked];

    assert(num_vars == 3);
    CCTK_POINTER output_arrays[3];
    output_arrays[0] = pt_betax;
    output_arrays[1] = pt_betay;
    output_arrays[2] = pt_betaz;

    // Interpolate
    int ierr;
    if (CCTK_IsFunctionAliased("InterpGridArrays")) {
      // TODO: use correct array types
      // (CCTK_POINTER[] vs. CCTK_REAL[])
      ierr = InterpGridArrays(cctkGH, dim, order, num_points, interp_coords,
                              num_vars, input_array_indices, num_vars,
                              output_arrays);
    } else {
      ierr = CCTK_InterpGridArrays(
          cctkGH, dim, operator_handle, param_table_handle, coordsys_handle,
          num_points, CCTK_VARIABLE_REAL, interp_coords, num_vars,
          input_array_indices, num_vars, output_array_type_codes,
          output_arrays);
    }
    if (ierr < 0) {
      CCTK_WARN(CCTK_WARN_ALERT, "Interpolation error");
      goto label_free_param_table;
    }

    if (CCTK_MyProc(cctkGH) == 0) {

      // Some more output

      if (verbose && CCTK_MyProc(cctkGH) == 0) {
        for (int n = 0; n < max_num_tracked; ++n) {
          if (track[n]) {
            CCTK_VINFO("Shift at puncture #%d is at (%g,%g,%g)", n,
                       double(pt_betax[n]), double(pt_betay[n]),
                       double(pt_betaz[n]));
          }
        }
      }

      // Check for NaNs and large shift components
      if (CCTK_MyProc(cctkGH) == 0) {
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
      }

      // Time evolution

      for (int n = 0; n < max_num_tracked; ++n) {
        if (track[n]) {
          const CCTK_REAL dt = pt_loc_t[n] - pt_loc_t_p[n];
          // First order time integrator
          // Michael Koppitz says this works...
          // if it doesn't, we can make it second order accurate
          pt_loc_x[n] = pt_loc_x_p[n] + dt * (-pt_betax[n]);
          pt_loc_y[n] = pt_loc_y_p[n] + dt * (-pt_betay[n]);
          pt_loc_z[n] = pt_loc_z_p[n] + dt * (-pt_betaz[n]);
          pt_vel_x[n] = -pt_betax[n];
          pt_vel_y[n] = -pt_betay[n];
          pt_vel_z[n] = -pt_betaz[n];
        }
      }
    }

    // Broadcast result

    CCTK_REAL loc_local[6 * max_num_tracked]; /* 3 components for location, 3
                                                 components for velocity */
    if (CCTK_MyProc(cctkGH) == 0) {
      for (int n = 0; n < max_num_tracked; ++n) {
        loc_local[n] = pt_loc_x[n];
        loc_local[max_num_tracked + n] = pt_loc_y[n];
        loc_local[2 * max_num_tracked + n] = pt_loc_z[n];
        loc_local[3 * max_num_tracked + n] = pt_vel_x[n];
        loc_local[4 * max_num_tracked + n] = pt_vel_y[n];
        loc_local[5 * max_num_tracked + n] = pt_vel_z[n];
      }
    } else {
      for (int n = 0; n < max_num_tracked; ++n) {
        loc_local[n] = 0.0;
        loc_local[max_num_tracked + n] = 0.0;
        loc_local[2 * max_num_tracked + n] = 0.0;
        loc_local[3 * max_num_tracked + n] = 0.0;
        loc_local[4 * max_num_tracked + n] = 0.0;
        loc_local[5 * max_num_tracked + n] = 0.0;
      }
    }

    CCTK_REAL loc_global[6 * max_num_tracked]; /* 3 components for location, 3
                                                  components for velocity */

    const int handle_sum = CCTK_ReductionArrayHandle("sum");
    if (handle_sum < 0) {
      CCTK_WARN(CCTK_WARN_ALERT, "Can't get redunction handle");
      goto label_free_param_table;
    }

    const int ierr2 = CCTK_ReduceLocArrayToArray1D(
        cctkGH, -1, handle_sum, loc_local, loc_global, 6 * max_num_tracked,
        CCTK_VARIABLE_REAL);
    if (ierr2 < 0) {
      CCTK_WARN(CCTK_WARN_ALERT, "Reduction error");
      goto label_free_param_table;
    }

    for (int n = 0; n < max_num_tracked; ++n) {
      pt_loc_x[n] = loc_global[n];
      pt_loc_y[n] = loc_global[max_num_tracked + n];
      pt_loc_z[n] = loc_global[2 * max_num_tracked + n];
      pt_vel_x[n] = loc_global[3 * max_num_tracked + n];
      pt_vel_y[n] = loc_global[4 * max_num_tracked + n];
      pt_vel_z[n] = loc_global[5 * max_num_tracked + n];
    }
  }

// Done

// Poor man's exception handling
label_free_param_table:
  Util_TableDestroy(param_table_handle);
}

extern "C" void PunctureTracker_SetPositions(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL dist;

  for (int n = 0; n < max_num_tracked; ++n) {
    if (track[n]) {
      // store puncture location in spherical surface
      if (which_surface_to_store_info[n] != -1) {
        int sn = which_surface_to_store_info[n];

        sf_centroid_x[sn] = pt_loc_x[n];
        sf_centroid_y[sn] = pt_loc_y[n];
        sf_centroid_z[sn] = pt_loc_z[n];

        sf_active[sn] = 1;
        sf_valid[sn] = 1;

        if (verbose) {
          CCTK_VINFO("Setting spherical surface %d centroid "
                     "from puncture #%d to (%g,%g,%g)",
                     sn, n, double(pt_loc_x[n]), double(pt_loc_y[n]),
                     double(pt_loc_z[n]));
        }
      }
    }
  }

  if (modify_puncture[0] >= 0 && modify_puncture[0] < max_num_tracked &&
      modify_puncture[1] >= 0 && modify_puncture[1] < max_num_tracked &&
      modify_puncture[0] != modify_puncture[1]) {

    if (track[modify_puncture[0]] && track[modify_puncture[1]]) {

      dist = sqrt(
          pow(pt_loc_x[modify_puncture[0]] - pt_loc_x[modify_puncture[1]], 2) +
          pow(pt_loc_y[modify_puncture[0]] - pt_loc_y[modify_puncture[1]], 2) +
          pow(pt_loc_z[modify_puncture[0]] - pt_loc_z[modify_puncture[1]], 2));

      if (dist < modify_distance) {

        if (new_reflevel_number[0] > -1) {
          if (verbose) {
            CCTK_VINFO("Setting the number of refinement levels to %d for "
                       "refinement region #%d",
                       new_reflevel_number[0], modify_puncture[0]);
          }
          num_levels[modify_puncture[0]] = new_reflevel_number[0];
        }

        if (new_reflevel_number[1] > -1) {
          if (verbose) {
            CCTK_VINFO("Setting the number of refinement levels to %d for "
                       "refinement region #%d",
                       new_reflevel_number[1], modify_puncture[1]);
          }
          num_levels[modify_puncture[1]] = new_reflevel_number[1];
        }
      }
    }
  }
}
}
