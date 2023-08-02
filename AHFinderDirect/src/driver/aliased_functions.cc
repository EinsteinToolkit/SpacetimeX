// horizon_radius.cc -- provide horizon radius & other info for other thorns
// $Header$
//
// <<<access to persistent data>>>
// AHFinderDirect_local_coordinate_origin - provide our local coordinate origin
// AHFinderDirect_horizon_was_found - query if a given horizon was found
// AHFinderDirect_horizon_centroid - query horizon centroid & if it was found
// AHFinderDirect_radius_in_direction - provide r(angle) function
//

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <vector>

#include "util_Table.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "config.h"
#include "stdc.h"
#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"
using jtutil::error_exit;

#include "../patch/coords.hh"
#include "../patch/grid.hh"
#include "../patch/fd_grid.hh"
#include "../patch/patch.hh"
#include "../patch/patch_edge.hh"
#include "../patch/patch_interp.hh"
#include "../patch/ghost_zone.hh"
#include "../patch/patch_system.hh"

#include "../elliptic/Jacobian.hh"

#include "../gr/gfns.hh"
#include "../gr/gr.hh"

#include "horizon_sequence.hh"
#include "BH_diagnostics.hh"
#include "driver.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// ***** access to persistent data *****
//
extern struct state state;

//******************************************************************************

//
// This function is called (via the magic of function aliasing) by
// other thorns to find out our local coordinate origin for a given AH.
//
// Results:
// This function returns 0 for ok, or -1 if the horizon number is invalid.
//
extern "C"
  CCTK_INT AHFinderDirect_local_coordinate_origin
    (CCTK_INT horizon_number,
     CCTK_REAL* p_origin_x, CCTK_REAL* p_origin_y, CCTK_REAL* p_origin_z)
{
const struct verbose_info& verbose_info = state.verbose_info;

if (!  ((horizon_number >= 1) && (horizon_number <= state.N_horizons))  )
   then {
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_local_coordinate_origin():\n"
"        horizon_number=%d must be in the range [1,N_horizons=%d]!\n"
		   ,
		   int(horizon_number), state.N_horizons);
	return -1;					// *** ERROR RETURN ***
	}

assert(state.AH_data_array[horizon_number] != NULL);
const struct AH_data& AH_data = *state.AH_data_array[horizon_number];

assert(AH_data.ps_ptr != NULL);
const patch_system& ps = *AH_data.ps_ptr;

assert(p_origin_x != NULL);
assert(p_origin_y != NULL);
assert(p_origin_z != NULL);
*p_origin_x = ps.origin_x();
*p_origin_y = ps.origin_y();
*p_origin_z = ps.origin_z();

return 0;						// *** NORMAL RETURN ***
}

//******************************************************************************

//
// This function is called (via the Cactus flesh function-aliasing mechanism)
// by other thorns to query whether or not the specified horizon was found
// the last time we searched for it.
//
// Results:
// This function returns
//  1 if the horizon was found
//  0 if the horizon was not found
//  -1 if the horizon number is invalid.
//
extern "C"
  CCTK_INT AHFinderDirect_horizon_was_found(CCTK_INT horizon_number)
{
const struct verbose_info& verbose_info = state.verbose_info;

if (!  ((horizon_number >= 1) && (horizon_number <= state.N_horizons))  )
   then {
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_horizon_was_found():\n"
"        horizon_number=%d must be in the range [1,N_horizons=%d]!\n"
		   ,
		   int(horizon_number), state.N_horizons);
	return -1;					// *** ERROR RETURN ***
	}

assert(state.AH_data_array[horizon_number] != NULL);
const struct AH_data& AH_data = *state.AH_data_array[horizon_number];

return AH_data.found_flag ? 1 : 0;
}

//******************************************************************************

//
// This function is called (via the Cactus flesh function-aliasing mechanism)
// by other thorns to query whether or not the specified horizon was found
// the last time we searched for it, and if so, to determine the horizon
// centroid.
//
// Results:
// This function returns:
//  1 if the horizon was found; in this case  *centroid_[xyz]_ptr
//    are set to the centroid position 
//  0 if the horizon was not found; in this case  *centroid_[xyz]_ptr
//    set to zeros
//  negative for an error
//
extern "C"
  CCTK_INT AHFinderDirect_horizon_centroid
    (CCTK_INT horizon_number,
     CCTK_REAL* p_centroid_x, CCTK_REAL* p_centroid_y, CCTK_REAL* p_centroid_z)
{
const struct verbose_info& verbose_info = state.verbose_info;

if (!  ((horizon_number >= 1) && (horizon_number <= state.N_horizons))  )
   then {
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_horizon_centroid():\n"
"        horizon_number=%d must be in the range [1,N_horizons=%d]!\n"
		   ,
		   int(horizon_number), state.N_horizons);
	return -1;					// *** ERROR RETURN ***
	}

assert(state.AH_data_array[horizon_number] != NULL);
const struct AH_data& AH_data = *state.AH_data_array[horizon_number];
const struct BH_diagnostics& BH_diagnostics = AH_data.BH_diagnostics;

assert(p_centroid_x != NULL);
assert(p_centroid_y != NULL);
assert(p_centroid_z != NULL);
if (AH_data.found_flag)
   then {
	*p_centroid_x = BH_diagnostics.centroid_x;
	*p_centroid_y = BH_diagnostics.centroid_y;
	*p_centroid_z = BH_diagnostics.centroid_z;
	return 1;
	}
   else {
	*p_centroid_x = 0.0;
	*p_centroid_y = 0.0;
	*p_centroid_z = 0.0;
	return 0;
	}
}

//******************************************************************************

//
// This function is called (via the Cactus flesh function-aliasing mechanism)
// by other thorns to find out a given AH's radius in the direction from
// its local coordinate origin to a given (x,y,z) coordinate or coordinates.
//
// Arguments:
// horizon_number = must be in the range 1 to N_horizons
// N_points = should be >= 0
// x[], y[], z[] = these give the (x,y,z) coordinates
// radius[] = this is set to the horizon radius values (Euclidean distance
//	      from the local coordinate origin), or to all -1.0 if we didn't
//	      find this horizon the last time we looked for it
//
// Results:
// This function returns 0 for ok, or -1 if the horizon number is invalid.
//
extern "C"
  CCTK_INT AHFinderDirect_radius_in_direction
    (CCTK_INT horizon_number,
     CCTK_INT N_points,
     const CCTK_REAL* const x, const CCTK_REAL* const y, const CCTK_REAL* const z,
     CCTK_REAL* const radius)
{
  const struct verbose_info& verbose_info = state.verbose_info;
  
  if (! ((horizon_number >= 1) && (horizon_number <= state.N_horizons)) ) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "AHFinderDirect_distance_outside_thorn():\n"
               "        horizon_number=%d must be in the range [1,N_horizons=%d]!\n"
               ,
               int(horizon_number), state.N_horizons);
    return -1;					// *** ERROR RETURN ***
  }
  
  assert(state.AH_data_array[horizon_number] != NULL);
  const struct AH_data& AH_data = *state.AH_data_array[horizon_number];
  
  if (AH_data.found_flag) {
    
    assert(AH_data.ps_ptr != NULL);
    const patch_system& ps = *AH_data.ps_ptr;
    
    std::vector<fp> local_xs(N_points);
    std::vector<fp> local_ys(N_points);
    std::vector<fp> local_zs(N_points);
    
    for (int point = 0 ; point < N_points ; ++point) {
      
      local_xs.at(point) = x[point] - ps.origin_x();
      local_ys.at(point) = y[point] - ps.origin_y();
      local_zs.at(point) = z[point] - ps.origin_z();
      
    }
    
    ps.radii_in_local_xyz_directions (gfns::gfn__h,
                                      N_points,
                                      & local_xs.front(),
                                      & local_ys.front(),
                                      & local_zs.front(),
                                      radius);
    
  } else {
    // if not found
    
    for (int point = 0 ; point < N_points ; ++point) {
      radius[point] = -1.0;
    }
    
  } // if not found

return 0;						// *** NORMAL RETURN ***
}

//******************************************************************************

	  }	// namespace AHFinderDirect
