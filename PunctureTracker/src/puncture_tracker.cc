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
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

namespace PunctureTracker {
using namespace std;
using namespace Loop;

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
    }
		else {
			pt_t_prev[n] = 0.0;
		}
  }

  // Interpolate

  // Dimensions
  const int dim = 3;

	// Number of interpolation variables
	int const num_vars = 3;

  const int operator_handle = 0;

  // Interpolation parameter table
  CCTK_INT operations[1][dim];
	for (int var = 0 ; var < num_vars; var++) {
	  operations[0][var] = 0;
	}

  int operands[1][dim];
	for (int var = 0 ; var < num_vars; var++) {
		operands[0][var] = var;
	}

	int ierr;
  // const int param_table_handle = Util_TableCreateFromString("order=4");
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (param_table_handle < 0)
    CCTK_VERROR("Can't create parameter table: %d", param_table_handle);
  if ((ierr = Util_TableSetInt(param_table_handle, interp_order, "order")) < 0)
    CCTK_VERROR("Can't set order in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, num_vars, (int const*const)operands,
                            "operand_indices")) < 0)
    CCTK_VERROR("Can't set operand_indices array in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, num_vars, (int const*const)operations,
                            "operation_codes")) < 0)
    CCTK_VERROR("Can't set operation_codes array in parameter table: %d", ierr);

  {

    // Interpolation coordinate system: Not used in CarpetX_DriverInterpolate
    const int coordsys_handle = 0;
		CCTK_INT const interp_coords_type_code = 0;

    // Only processor 0 interpolates
    const int num_points = CCTK_MyProc(cctkGH) == 0 ? max_num_tracked : 0;

    // Interpolation coordinates
    assert(dim == 3);
    CCTK_POINTER_TO_CONST interp_coords[dim];
		interp_coords[0] = pt_loc_x;
		interp_coords[1] = pt_loc_y;
		interp_coords[2] = pt_loc_z;

		// const CCTK_REAL interp_coords[dim][num_points] = {*pt_loc_x_p, *pt_loc_y_p, *pt_loc_z_p};
		// const void* interp_coords[dim] = {*pt_loc_x_p, *pt_loc_y_p, *pt_loc_z_p};

    // Interpolated variables
    assert(num_vars == 3);
    int input_array_indices[3];
    input_array_indices[0] = CCTK_VarIndex("ADMBaseX::betax");
    input_array_indices[1] = CCTK_VarIndex("ADMBaseX::betay");
    input_array_indices[2] = CCTK_VarIndex("ADMBaseX::betaz");

    // Interpolation result types: Not used by CarpetX DriverInterp
		CCTK_INT const output_array_type_codes[1] = {0};

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
    // Use CarpetX Funtion:
		ierr = DriverInterpolate(
		cctkGH, dim, operator_handle, param_table_handle, coordsys_handle,
		num_points, interp_coords_type_code, interp_coords, num_vars, (int const * const)input_array_indices,
		num_vars, output_array_type_codes, output_arrays);

		// Interpolate(cctkGH, num_points, interp_coords[0], interp_coords[1], interp_coords[2],
		// 	num_vars, (CCTK_INT const * const)input_array_indices, (CCTK_INT const * const)operations,
		// 	(CCTK_REAL **)output_arrays); 

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

			// CCTK_VINFO("Not updating puncture locations!");
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

    // Broadcast result

    CCTK_REAL loc_global[6 * max_num_tracked]; /* 3 components for location, 3
                                                 components for velocity */
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
// Done

// Poor man's exception handling
label_free_param_table:
  Util_TableDestroy(param_table_handle);
}

using namespace Arith;

extern "C" void PunctureTracker_CheckShift(CCTK_ARGUMENTS) {
	DECLARE_CCTK_ARGUMENTS_PunctureTracker_CheckShift;
	DECLARE_CCTK_PARAMETERS;

  const int dim = 3;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> betax_(layout, betax);
  const GF3D2<const CCTK_REAL> betay_(layout, betay);
  const GF3D2<const CCTK_REAL> betaz_(layout, betaz);

  const GridDescBaseDevice grid(cctkGH);

	const int level = ilogb(CCTK_REAL(cctk_levfac[0]));
	const int finest_lvl = num_levels[0] - 1;	

	if (level == finest_lvl) {
		for (int n = 0; n < max_num_tracked; ++n) {
			if (track[n]) {

				const vect<CCTK_REAL, dim> loc_vec = {pt_loc_x[n], pt_loc_y[n], pt_loc_z[n]};

				grid.loop_all_device<0, 0, 0>(grid.nghostzones,
																			[=] CCTK_DEVICE(const PointDesc &p)
																					CCTK_ATTRIBUTE_ALWAYS_INLINE {
						if (maximum(abs(p.X - loc_vec)) <= (1 / pow(2, level))) {
							CCTK_VINFO("Shift at level %d near puncture #%d is {%g, %g, %g} at coords {%g, %g, %g}.", level, n,
								betax_(p.I), betay_(p.I), betaz_(p.I), p.x, p.y, p.z); 
						}
					});
			}
		}
	}
}

extern "C" void CheckInterpolate(CCTK_ARGUMENTS) {
	DECLARE_CCTK_ARGUMENTS_CheckInterpolate;
	DECLARE_CCTK_PARAMETERS;

	const int level = ilogb(CCTK_REAL(cctk_levfac[0]));
	const int finest_lvl = num_levels[0] - 1;	

	const CCTK_REAL dgrid = 1 / pow(2, level);

  // Interpolate

  // Dimensions
  const int dim = 3;

	// Number of interpolation variables
	int const num_vars = 3;

  const int operator_handle = 0;

  // Interpolation parameter table
  CCTK_INT operations[1][dim];
	for (int var = 0 ; var < num_vars; var++) {
	  operations[0][var] = 0;
	}

  int operands[1][dim];
	for (int var = 0 ; var < num_vars; var++) {
		operands[0][var] = var;
	}

	int ierr;
  // const int param_table_handle = Util_TableCreateFromString("order=4");
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (param_table_handle < 0)
    CCTK_VERROR("Can't create parameter table: %d", param_table_handle);
  if ((ierr = Util_TableSetInt(param_table_handle, interp_order, "order")) < 0)
    CCTK_VERROR("Can't set order in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, num_vars, (int const*const)operands,
                            "operand_indices")) < 0)
    CCTK_VERROR("Can't set operand_indices array in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, num_vars, (int const*const)operations,
                            "operation_codes")) < 0)
    CCTK_VERROR("Can't set operation_codes array in parameter table: %d", ierr);

  {

    // Interpolation coordinate system: Not used in CarpetX_DriverInterpolate
    const int coordsys_handle = 0;
		CCTK_INT const interp_coords_type_code = 0;

    // Only processor 0 interpolates
    const int num_points = CCTK_MyProc(cctkGH) == 0 ? 2 : 0;

    // Interpolation coordinates
    assert(dim == 3);
    const CCTK_REAL coords[dim][2] = {
			{4, -4}, {dgrid, dgrid}, {dgrid, dgrid}
		};
		const void* interp_coords[dim] = {
			coords[0], coords[1], coords[2]
		};


    // Interpolated variables
    assert(num_vars == 3);
    int input_array_indices[3];
    input_array_indices[0] = CCTK_VarIndex("ADMBaseX::betax");
    input_array_indices[1] = CCTK_VarIndex("ADMBaseX::betay");
    input_array_indices[2] = CCTK_VarIndex("ADMBaseX::betaz");

    // Interpolation result types: Not used by CarpetX DriverInterp
		CCTK_INT const output_array_type_codes[1] = {0};

    // Interpolation result
    CCTK_REAL pt_betax[2];
    CCTK_REAL pt_betay[2];
    CCTK_REAL pt_betaz[2];

    assert(num_vars == 3);
    CCTK_POINTER output_arrays[3];
    output_arrays[0] = pt_betax;
    output_arrays[1] = pt_betay;
    output_arrays[2] = pt_betaz;

    // Interpolate
    int ierr;
    // Use CarpetX Funtion:
		ierr = DriverInterpolate(
		cctkGH, dim, operator_handle, param_table_handle, coordsys_handle,
		num_points, interp_coords_type_code, interp_coords, num_vars, (int const * const)input_array_indices,
		num_vars, output_array_type_codes, output_arrays);

    if (CCTK_MyProc(cctkGH) == 0) {

      // Some more output

      if (verbose && CCTK_MyProc(cctkGH) == 0) {
      	CCTK_VINFO("Shift at x=+4 interpolated to be (%g,%g,%g)", 
                   double(pt_betax[0]), double(pt_betay[0]),
                   double(pt_betaz[0]));
      	CCTK_VINFO("Shift at x=-4 interpolated to be (%g,%g,%g)", 
                   double(pt_betax[1]), double(pt_betay[1]),
                   double(pt_betaz[1]));
      }
		}
	}

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> betax_(layout, betax);
  const GF3D2<const CCTK_REAL> betay_(layout, betay);
  const GF3D2<const CCTK_REAL> betaz_(layout, betaz);

  const GridDescBaseDevice grid(cctkGH);

	if (level == finest_lvl) {
		const vect<CCTK_REAL, dim> loc_vec1 = {4, dgrid, dgrid};
		const vect<CCTK_REAL, dim> loc_vec2 = {-4, dgrid, dgrid};

		grid.loop_all_device<0, 0, 0>(grid.nghostzones,
																	[=] CCTK_DEVICE(const PointDesc &p)
																			CCTK_ATTRIBUTE_ALWAYS_INLINE {
				if (maximum(abs(p.X - loc_vec1)) <= (1 / pow(2, level))) {
					CCTK_VINFO("Shift near x = +4 is {%g, %g, %g} at coords {%g, %g, %g}.", 
						betax_(p.I), betay_(p.I), betaz_(p.I), p.x, p.y, p.z); 
				}
				else if (maximum(abs(p.X - loc_vec2)) <= (1 / pow(2, level))) {
					CCTK_VINFO("Shift near x = -4 is {%g, %g, %g} at coords {%g, %g, %g}.", 
						betax_(p.I), betay_(p.I), betaz_(p.I), p.x, p.y, p.z); 
				}

			});
	}
}

} //namespace PunctureTracker
