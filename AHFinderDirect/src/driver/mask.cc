// mask.cc -- set a mask gridfn based on each horizon shape
// $Header$
//
// <<<access to persistent data>>>
/// <<<data structures local to this file>>>
/// <<<prototypes for functions local to this file>>>
// AHFinderDirect_do_masks - top-level driver for all mask stuff
/// setup_mask_grid_info - setup mask grid origin/delta etc
/// setup_mask_dataptrs_and_bitfields - map gridfn/bitfield names to ptr/bitmask
/// set_mask_gridfn - set mask gridfn(s) based on each horizon's shape
/// set_mask_gridfn_to_outside_value - ... "outside" value
/// set_mask_gridfn_to_inside_and_buffer_values - "inside"/"buffer" values
/// inner_mask_radius - compute inner mask radius r_inner given r_horizon
/// outer_mask_radius - compute outer mask radius r_outer given r_inner
//

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "util_Table.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "SpaceMask.h"			// from thorn SpaceMask

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

// define this to get extra debugging on the mask grid origin/delta etc
#undef DEBUG_MASK_GRID

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
// ***** data structures local to this file *****
//

//
// This structure holds all the information we need about the (Cactus)
// grid where the mask gridfn(s) is/are stored.
//
namespace {
struct	mask_grid_info
	{
	const cGH *GH;			// --> Cactus grid hierarchy

	// Cactus coordinate system
	// (for the *current* grid if we are doing mesh refinement)
	fp proc_coord_origin[N_GRID_DIMS];	// global (x,y,z) of
						// *this processor's*
						// (i,j,k) = (0,0,0)
						// grid point
	fp coord_delta[N_GRID_DIMS];		// (x,y,z) grid spacing

	// maximum of x,y,z grid spacings
	fp max_coord_delta;

	// geometric mean of x,y,z grid spacings,
	// on the *base* grid if we are doing mesh refinement
	// (we need the base-grid semantics to make excision consistent
	//  across refinement levels)
	fp base_grid_mean_coord_delta;

	// dimensions of gridfn data on this processor, viewed as a 3-D array
	// n.b. storage ordering is Fortran,
	//	i.e. i is contiguous, j has stride Ni, k has stride Ni*Nj
	CCTK_INT proc_gridfn_dims[N_GRID_DIMS];


	//
	// coordinate conversion functions
	//
	// here ijk = gridfn array indices on this processor
	//      xyz = floating-point global coordinates
	//

	// convert integer ijk --> floating-point global xyz coordinates
	fp global_xyz_of_ijk(int axis, int ijk) const
		{ return proc_coord_origin[axis] + ijk*coord_delta[axis]; }
	// convert floating-point global xyz --> integer ijk coordinates
	// ... but as a floating-point number
	fp fp_ijk_of_global_xyz(int axis, fp xyz) const
		{ return (xyz - proc_coord_origin[axis]) / coord_delta[axis]; }
	// convert floating-point global_xyz --> integer ijk coordinates
	// ... rounding down/up to the next integer (= to the next grid point)
	int ijk_floor_of_global_xyz(int axis, fp xyz) const
		{ return jtutil::ifloor(fp_ijk_of_global_xyz(axis, xyz)); }
	int ijk_ceil_of_global_xyz (int axis, fp xyz) const
		{ return jtutil::iceil(fp_ijk_of_global_xyz(axis, xyz)); }
	};
	  }

//******************************************************************************

//
// ***** prototypes for functions local to this file *****
//

namespace {
void setup_mask_grid_info(CCTK_ARGUMENTS, struct mask_grid_info& mgi);
void setup_mask_dataptrs_and_bitfields(const cGH *GH,
				       struct mask_info& mask_info);
void set_mask_gridfn(int N_horizons,
		     const struct AH_data* const AH_data_array[],
		     const struct mask_grid_info& mgi,
		     const struct mask_info& mask_info,
		     const struct verbose_info& verbose_info);

void set_mask_gridfn_to_outside_value(const struct mask_grid_info& mgi,
				      const struct mask_info& mask_info,
				      const struct verbose_info& verbose_info);
void set_mask_gridfn_to_inside_and_buffer_values
	(const struct mask_grid_info& mgi,
	 int use_min_i, int use_max_i,
	 int use_min_j, int use_max_j,
	 int use_min_k, int use_max_k,
	 const struct mask_info& mask_info,
	 const patch_system& ps,
	 const struct verbose_info& verbose_info);
fp inner_mask_radius(const struct mask_grid_info& mgi,
		     const struct mask_info& mask_info,
		     fp r_horizon);
fp outer_mask_radius(const struct mask_grid_info& mgi,
		     const struct mask_info& mask_info,
		     fp r_inner);
	  }

//******************************************************************************

//
// This function is called by the Cactus scheduler after we find any
// apparent horizons, to do all this thorn's mask processing.
//
extern "C"
  void AHFinderDirect_maybe_do_masks(CCTK_ARGUMENTS)
{
DECLARE_CCTK_ARGUMENTS_AHFinderDirect_maybe_do_masks
DECLARE_CCTK_PARAMETERS

const struct verbose_info& verbose_info = state.verbose_info;
      struct    mask_info&    mask_info = state.mask_info;

// optionally set the mask gridfn based on each horizon's shape
if (mask_info.set_mask_for_any_horizon)
   then {
	//
	// this setup has to be done each time we're called, and
	// the  mask_grid_info  structure can't live in  struct state ,
	// because if mesh refinement is in effect we'll be called
	// separately for each separate locally-uniform grid patch
	//
	mask_grid_info mgi;
	setup_mask_grid_info(CCTK_PASS_CTOC, mgi);

	setup_mask_dataptrs_and_bitfields(cctkGH, mask_info);
	set_mask_gridfn(state.N_horizons, state.AH_data_array,
			mgi, mask_info,
			verbose_info);
	}
}

//******************************************************************************

//
// This function sets up the  mask_grid_info  structure to describe
// the (current) Cactus grid where the mask gridfn(s) live.
//
namespace {
void setup_mask_grid_info(CCTK_ARGUMENTS, struct mask_grid_info& mgi)
{
DECLARE_CCTK_ARGUMENTS

mgi.GH = cctkGH;

// Cactus grid spacing on *this* grid
mgi.coord_delta[X_AXIS] = CCTK_DELTA_SPACE(X_AXIS);
mgi.coord_delta[Y_AXIS] = CCTK_DELTA_SPACE(Y_AXIS);
mgi.coord_delta[Z_AXIS] = CCTK_DELTA_SPACE(Z_AXIS);
mgi.max_coord_delta = jtutil::max(mgi.coord_delta[X_AXIS],
		      jtutil::max(mgi.coord_delta[Y_AXIS],
				  mgi.coord_delta[Z_AXIS]));

// Cactus grid spacings on the *base* grid
const fp base_grid_delta_product =   cctk_delta_space[X_AXIS]
				   * cctk_delta_space[Y_AXIS]
				   * cctk_delta_space[Z_AXIS];
mgi.base_grid_mean_coord_delta = pow(base_grid_delta_product, 1.0/3.0);

// get global/local Cactus grid origin
// KLUDGE -- is this the right way to get this??
mgi.proc_coord_origin[X_AXIS] = CCTK_ORIGIN_SPACE(X_AXIS)
				 + cctk_lbnd[X_AXIS] * mgi.coord_delta[X_AXIS];
mgi.proc_coord_origin[Y_AXIS] = CCTK_ORIGIN_SPACE(Y_AXIS)
				 + cctk_lbnd[Y_AXIS] * mgi.coord_delta[Y_AXIS];
mgi.proc_coord_origin[Z_AXIS] = CCTK_ORIGIN_SPACE(Z_AXIS)
				 + cctk_lbnd[Z_AXIS] * mgi.coord_delta[Z_AXIS];
mgi.proc_gridfn_dims[X_AXIS] = cctk_lsh[X_AXIS];
mgi.proc_gridfn_dims[Y_AXIS] = cctk_lsh[Y_AXIS];
mgi.proc_gridfn_dims[Z_AXIS] = cctk_lsh[Z_AXIS];

#ifdef DEBUG_MASK_GRID
printf("mask.cc:: cctk_lsh[] = [%d,%d,%d]\n",
	int(cctk_lsh[X_AXIS]),
	int(cctk_lsh[Y_AXIS]),
	int(cctk_lsh[Z_AXIS]));
printf("mask.cc:: mgi.coord_delta[] = [%g,%g,%g]\n",
	double(mgi.coord_delta[X_AXIS]),
	double(mgi.coord_delta[Y_AXIS]),
	double(mgi.coord_delta[Z_AXIS]));
printf("mask.cc:: mgi.proc_coord_origin[] = [%g,%g,%g]\n",
	double(mgi.proc_coord_origin[X_AXIS]),
	double(mgi.proc_coord_origin[Y_AXIS]),
	double(mgi.proc_coord_origin[Z_AXIS]));
#endif
}
	  }

//******************************************************************************

//
// This function maps the character-string names of the mask gridfn(s)
// and/or bitfield(s) into internal data pointers and/or bit masks/values:
//
//	mask_info.old_style_mask_info
//		 .gridfn_varindex --> .gridfn_dataptr
//	mask_info.new_style_mask_info
//		 .gridfn_varindex --> .gridfn_dataptr
//		 .bitfield_name   --> .bitfield_bitmask
//		 .inside_value    --> .inside_bitvalue
//		 .buffer_value    --> .buffer_bitvalue
//		 .outside_value   --> .outside_bitvalue
//
namespace {
void setup_mask_dataptrs_and_bitfields(const cGH *GH,
				       struct mask_info& mask_info)
{
const int mask_time_level = 0;
const int SPACEMASK_ERROR = -1;

struct mask_info::old_style_mask_info& osmi = mask_info.old_style_mask_info;
if (mask_info.set_old_style_mask)
   then {
	osmi.gridfn_dataptr = static_cast<CCTK_REAL*>(
				CCTK_VarDataPtrI(GH, mask_time_level,
						 osmi.gridfn_varindex)
						     );
	if (osmi.gridfn_dataptr == NULL)
	   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"setup_mask_dataptrs_and_bitfields():\n"
"        couldn't get data pointer for old-style mask \"%s\"!\n"
			   ,
			   osmi.gridfn_name);			/*NOTREACHED*/
	}

if (mask_info.set_new_style_mask)
   then {
	struct mask_info::new_style_mask_info& nsmi
		= mask_info.new_style_mask_info;
	nsmi.gridfn_dataptr = static_cast<CCTK_INT*>(
				CCTK_VarDataPtrI(GH, mask_time_level,
						 nsmi.gridfn_varindex)
						    );
	if (nsmi.gridfn_dataptr == NULL)
	   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"setup_mask_dataptrs_and_bitfields():\n"
"        couldn't get data pointer for new-style mask \"%s\"!\n"
			   ,
			   nsmi.gridfn_name);			/*NOTREACHED*/

	nsmi.bitfield_bitmask = SpaceMask_GetTypeBits(nsmi.bitfield_name);
	if (nsmi.bitfield_bitmask == SPACEMASK_ERROR)
	   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"setup_mask_dataptrs_and_bitfields():\n"
"        couldn't get bit mask for new-style mask bit field \"%s\"!\n"
"        (this almost certainly means this bit field was never registered\n"
"         with SpaceMask; remember that AHFinderDirect doesn't do this\n"
"         registration, rather you need to arrange for some other thorn(s)\n"
"         to do it)\n"
		   ,
		   nsmi.bitfield_name);				/*NOTREACHED*/

	nsmi.inside_bitvalue  = SpaceMask_GetStateBits(nsmi.bitfield_name,
						       nsmi.inside_value);
	nsmi.buffer_bitvalue  = SpaceMask_GetStateBits(nsmi.bitfield_name,
						       nsmi.buffer_value);
	nsmi.outside_bitvalue = SpaceMask_GetStateBits(nsmi.bitfield_name,
						       nsmi.outside_value);
	if (    (nsmi.inside_bitvalue == SPACEMASK_ERROR)
	     || (nsmi.buffer_bitvalue == SPACEMASK_ERROR)
	     || (nsmi.outside_bitvalue == SPACEMASK_ERROR)    )
	   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"setup_mask_dataptrs_and_bitfields():\n"
"        couldn't get bit value(s) for one or more of the new-style mask\n"
"        inside/buffer/outside values \"%s\"/\"%s\"/\"%s\"!\n"
"        (this almost certainly means one or more of the values was/were\n"
"         never registered with SpaceMask; remember that AHFinderDirect\n"
"         doesn't do this registration, rather you need to arrange for\n"
"         some other thorn(s) to do it)\n"
			   ,
			   nsmi.inside_value,
			   nsmi.buffer_value,
			   nsmi.outside_value);			/*NOTREACHED*/
	}
}
	  }

//******************************************************************************

//
// This function sets the mask grid function specified by  mask_info
// based on each horizon's shape.
//
namespace {
void set_mask_gridfn(int N_horizons,
		     const struct AH_data* const AH_data_array[],
		     const struct mask_grid_info& mgi,
		     const struct mask_info& mask_info,
		     const struct verbose_info& verbose_info)
{
const bool set_old_style_mask = mask_info.set_old_style_mask;
const bool set_new_style_mask = mask_info.set_new_style_mask;
const struct mask_info::old_style_mask_info& osmi = mask_info.old_style_mask_info;
const struct mask_info::new_style_mask_info& nsmi = mask_info.new_style_mask_info;

if (verbose_info.print_algorithm_highlights)
   then {
	if (set_old_style_mask)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "setting old-style (CCTK_REAL) mask grid function %s",
			   osmi.gridfn_name);
	if (set_new_style_mask)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "setting new-style (CCTK_INT) mask grid function %s",
			   nsmi.gridfn_name);
		CCTK_VInfo(CCTK_THORNSTRING,
			   "                                  bit field %s",
			   nsmi.bitfield_name);
		}
	}

if (verbose_info.print_algorithm_debug)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "      grid on this processor has x=[%g,%g]",
		   double(mgi.global_xyz_of_ijk(X_AXIS, 0)),
		   double(mgi.global_xyz_of_ijk(X_AXIS, mgi.proc_gridfn_dims[X_AXIS]-1)));
	CCTK_VInfo(CCTK_THORNSTRING,
		   "                                 y=[%g,%g]",
		   double(mgi.global_xyz_of_ijk(Y_AXIS, 0)),
		   double(mgi.global_xyz_of_ijk(Y_AXIS, mgi.proc_gridfn_dims[Y_AXIS]-1)));
	CCTK_VInfo(CCTK_THORNSTRING,
		   "                                 z=[%g,%g]",
		   double(mgi.global_xyz_of_ijk(Z_AXIS, 0)),
		   double(mgi.global_xyz_of_ijk(Z_AXIS, mgi.proc_gridfn_dims[Z_AXIS]-1)));
	}


//
// set the mask to the outside value everywhere in
// (this processor's chunk of) the grid
//
set_mask_gridfn_to_outside_value(mgi,
				 mask_info,
				 verbose_info);


//
// loop over each horizon's xyz bounding box
// (intersected with (this processor's chunk of) the grid)
// to set the mask accurately (or skip this horizon if it's too small)
//
	for (int hn = 1 ; hn <= N_horizons ; ++hn)
	{
	if (! mask_info.set_mask_for_this_horizon[hn])
	   then continue;				// *** LOOP CONTROL ***

	const struct AH_data& AH_data = *AH_data_array[hn];
	if (! AH_data.found_flag)
	   then continue;				// *** LOOP CONTROL ***

	const patch_system& ps                      = *AH_data.ps_ptr;
	const struct BH_diagnostics& BH_diagnostics =  AH_data.BH_diagnostics;

	//
	// skip this horizon if it's too small
	//
	// ... The condition is documented as using the minimum over angle
	//     of r_inner(r_horizon(angle)).  But this is trivially equivalent
	//     to r_inner(the minimum over angle of r_horizon), i.e. to
	//     r_inner(BH_diagnostics.min_radius); this latter form is what
	//     we actually compute.
	//
	const fp r_inner_min
		= inner_mask_radius(mgi, mask_info, BH_diagnostics.min_radius);
	if (r_inner_min
	    < mask_info.min_horizon_radius_points_for_mask*mgi.max_coord_delta)
	   then {
		if (verbose_info.print_algorithm_details)
		   then CCTK_VInfo(CCTK_THORNSTRING,
   "   r_inner_min=%g < %g grid points ==> skipping mask for horizon %d",
				   r_inner_min,
				   mask_info.min_horizon_radius_points_for_mask,
				   hn);
		continue;				// *** LOOP CONTROL ***
		}

	//
	// get to here ==> normal mask processing for this horizon
	//

	if (verbose_info.print_algorithm_details)
	   then CCTK_VInfo(CCTK_THORNSTRING,
   "   setting mask grid function to \"buffer\"/\"inside\" for horizon %d",
			   hn);


	// horizon bounding box, rounded "out" to the next grid point
	const int AH_min_i = mgi.ijk_floor_of_global_xyz(X_AXIS, BH_diagnostics.min_x);
	const int AH_max_i = mgi.ijk_ceil_of_global_xyz (X_AXIS, BH_diagnostics.max_x);
	const int AH_min_j = mgi.ijk_floor_of_global_xyz(Y_AXIS, BH_diagnostics.min_y);
	const int AH_max_j = mgi.ijk_ceil_of_global_xyz (Y_AXIS, BH_diagnostics.max_y);
	const int AH_min_k = mgi.ijk_floor_of_global_xyz(Z_AXIS, BH_diagnostics.min_z);
	const int AH_max_k = mgi.ijk_ceil_of_global_xyz (Z_AXIS, BH_diagnostics.max_z);
	if (verbose_info.print_algorithm_debug)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "      horizon bounding box is x=[%g,%g]",
			   double(BH_diagnostics.min_x),
			   double(BH_diagnostics.max_x));
		CCTK_VInfo(CCTK_THORNSTRING,
			   "                              y=[%g,%g]",
			   double(BH_diagnostics.min_y),
			   double(BH_diagnostics.max_y));
		CCTK_VInfo(CCTK_THORNSTRING,
			   "                              z=[%g,%g]",
			   double(BH_diagnostics.min_z),
			   double(BH_diagnostics.max_z));
		}

	//
	// compute intersection of horizon bounding box with
	// (this processor's chunk of) the grid
	//
	// as can be seen from the following picture, the intersection
	// has the max/min of the min/max coordinate as its min/max:
	//
	//	           +--------------------+
	//	           |                    |
	//	+----------+--------------+     |
	//	|          |XXXXXXXXXXXXXX|     |
	//	|          |XXXXXXXXXXXXXX|     |
	//	|          |XXXXXXXXXXXXXX|     |
	//	|          +--------------+-----+
	//	|                         |
	//	|                         |
	//	+-------------------------+
	//
	const int use_min_i = jtutil::max(AH_min_i, 0);
	const int use_max_i = jtutil::min(AH_max_i, int(mgi.proc_gridfn_dims[X_AXIS]-1));
	const int use_min_j = jtutil::max(AH_min_j, 0);
	const int use_max_j = jtutil::min(AH_max_j, int(mgi.proc_gridfn_dims[Y_AXIS]-1));
	const int use_min_k = jtutil::max(AH_min_k, 0);
	const int use_max_k = jtutil::min(AH_max_k, int(mgi.proc_gridfn_dims[Z_AXIS]-1));

	if (verbose_info.print_algorithm_debug)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "      use bounding box is x=[%g,%g]",
			   double(mgi.global_xyz_of_ijk(X_AXIS, use_min_i)),
			   double(mgi.global_xyz_of_ijk(X_AXIS, use_max_i)));
		CCTK_VInfo(CCTK_THORNSTRING,
			   "                          y=[%g,%g]",
			   double(mgi.global_xyz_of_ijk(Y_AXIS, use_min_j)),
			   double(mgi.global_xyz_of_ijk(Y_AXIS, use_max_j)));
		CCTK_VInfo(CCTK_THORNSTRING,
			   "                          z=[%g,%g]",
			   double(mgi.global_xyz_of_ijk(Z_AXIS, use_min_k)),
			   double(mgi.global_xyz_of_ijk(Z_AXIS, use_max_k)));
		}

	set_mask_gridfn_to_inside_and_buffer_values(mgi,
						    use_min_i, use_max_i,
						    use_min_j, use_max_j,
						    use_min_k, use_max_k,
						    mask_info,
						    ps,
						    verbose_info);
	}
}
	  }

//******************************************************************************

//
// For each Cactus grid point in (this processor's chunk of) the Cactus
// grid, this function sets the mask grid function specified by  mask_info
// to the appropriate value for the "outside" region (possibly leaving
// the mask unchanged at its previous value if the  mask_is_noshrink
// option is set).
//
namespace {
void set_mask_gridfn_to_outside_value(const struct mask_grid_info& mgi,
				      const struct mask_info& mask_info,
				      const struct verbose_info& verbose_info)
{
const bool set_old_style_mask = mask_info.set_old_style_mask;
const bool set_new_style_mask = mask_info.set_new_style_mask;
const struct mask_info::old_style_mask_info& osmi = mask_info.old_style_mask_info;
const struct mask_info::new_style_mask_info& nsmi = mask_info.new_style_mask_info;


if (verbose_info.print_algorithm_details)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   setting mask grid function to \"outside\"");

	for (int k = 0 ; k < mgi.proc_gridfn_dims[Z_AXIS] ; ++k)
	{
	for (int j = 0 ; j < mgi.proc_gridfn_dims[Y_AXIS] ; ++j)
	{
	for (int i = 0 ; i < mgi.proc_gridfn_dims[X_AXIS] ; ++i)
	{
	const int posn = CCTK_GFINDEX3D(mgi.GH, i,j,k);
	if (set_old_style_mask)
	   then {
		const CCTK_REAL old_value = osmi.gridfn_dataptr[posn];
		if ( mask_info.mask_is_noshrink
		     && (    (old_value == osmi.inside_value)
			  || (old_value == osmi.buffer_value)    ) )
		   then {
			// we're in "noshrink mode"
			// and this point was previously "inside" or "buffer"
			// ==> no-op here
			}
		   else osmi.gridfn_dataptr[posn] = osmi.outside_value;
		}
	if (set_new_style_mask)
	   then {
		if ( mask_info.mask_is_noshrink
		     && (    SpaceMask_CheckStateBits(nsmi.gridfn_dataptr, posn,
						      nsmi.bitfield_bitmask,
						      nsmi.inside_bitvalue)
			  || SpaceMask_CheckStateBits(nsmi.gridfn_dataptr, posn,
						      nsmi.bitfield_bitmask,
						      nsmi.buffer_bitvalue)    ) )
		   then {
			// we're in "noshrink mode"
			// and this point was previously "inside" or "buffer"
			// ==> no-op here
			}
		   else SpaceMask_SetStateBits(nsmi.gridfn_dataptr, posn,
					       nsmi.bitfield_bitmask,
					       nsmi.outside_bitvalue);
		}
	}
	}
	}
}
	  }

//******************************************************************************

//
// For each Cactus grid point in a specified region of the Cactus grid,
// and for a specified (single) horizon, this function either sets the
// mask grid function specified by  mask_info to the appropriate value
// if the point is in the "inside" or "buffer" region for a specified
// horizon (possibly leaving the mask unchanged at its previous value
// if the  mask_is_noshrink  option is set), or leaves the mask grid
// function unchanged if the point is actually outside this horizon.
//
// Arguments:
// ps = The patch system for this horizon.
// use_{min,max}_{i,j,k} = Specifies the region of the Cactus grid.
//
// FIXME:
// The current implementation is quite inefficient: it does a full 2-D
// interpolation (to find the horizon radius) for each xyz grid point
// in the specified region.  It would be more efficient to batch the
// interpolations, and maybe also use r_min/r_max tests for early-out.
//
namespace {
void set_mask_gridfn_to_inside_and_buffer_values
	(const struct mask_grid_info& mgi,
	 int use_min_i, int use_max_i,
	 int use_min_j, int use_max_j,
	 int use_min_k, int use_max_k,
	 const struct mask_info& mask_info,
	 const patch_system& ps,
	 const struct verbose_info& verbose_info)
{
const bool set_old_style_mask = mask_info.set_old_style_mask;
const bool set_new_style_mask = mask_info.set_new_style_mask;
const struct mask_info::old_style_mask_info& osmi = mask_info.old_style_mask_info;
const struct mask_info::new_style_mask_info& nsmi = mask_info.new_style_mask_info;

long inside_count = 0;
long buffer_count = 0;

	for (int k = use_min_k ; k <= use_max_k ; ++k)
	{
	for (int j = use_min_j ; j <= use_max_j ; ++j)
	{
	for (int i = use_min_i ; i <= use_max_i ; ++i)
	{
	const int posn = CCTK_GFINDEX3D(mgi.GH, i,j,k);

	const fp global_x = mgi.global_xyz_of_ijk(X_AXIS, i);
	const fp global_y = mgi.global_xyz_of_ijk(Y_AXIS, j);
	const fp global_z = mgi.global_xyz_of_ijk(Z_AXIS, k);
	const fp local_x = global_x - ps.origin_x();
	const fp local_y = global_y - ps.origin_y();
	const fp local_z = global_z - ps.origin_z();

	const fp r = jtutil::hypot3(local_x, local_y, local_z);

	// interpolate to find r_horizon in this direction
	//
	// FIXME:
	// it would be more efficient here to compute the
	// radii of a whole batch of points at once
	const fp r_horizon = ps.radius_in_local_xyz_direction
				   (gfns::gfn__h,
				    local_x, local_y, local_z);

	const fp r_inner = inner_mask_radius(mgi, mask_info, r_horizon);
	const fp r_outer = outer_mask_radius(mgi, mask_info, r_inner);

	if	(r <= r_inner)
	   then {
		// set the mask to the "inside" value
		++inside_count;
		if (set_old_style_mask)
		   then osmi.gridfn_dataptr[posn] = osmi.inside_value;
		if (set_new_style_mask)
		   then SpaceMask_SetStateBits(nsmi.gridfn_dataptr, posn,
					       nsmi.bitfield_bitmask,
					       nsmi.inside_bitvalue);
		}
	else if (r <= r_outer)
	   then {
		// set the mask to the "buffer" value
		// ... except that if the  mask_is_noshrink  option is set
		//     and the mask is already set to the "inside" value
		//     then we need to leave it unchanged
		//     (at the "inside" value)
		++buffer_count;
		if (set_old_style_mask)
		   then {
			if ( mask_info.mask_is_noshrink
			     && (osmi.gridfn_dataptr[posn] == osmi.inside_value) )
			   then {
				// we're in "noshrink mode"
				// and this point was previously "inside"
				// ==> no-op here
				}
			   else osmi.gridfn_dataptr[posn] = osmi.buffer_value;
			}
		if (set_new_style_mask)
		   then {
			if ( mask_info.mask_is_noshrink
			     && SpaceMask_CheckStateBits(nsmi.gridfn_dataptr, posn,
							 nsmi.bitfield_bitmask,
							 nsmi.inside_bitvalue) )
			   then {
				// we're in "noshrink mode"
				// and this point was previously "inside"
				// ==> no-op here
				}
			   else SpaceMask_SetStateBits(nsmi.gridfn_dataptr, posn,
						       nsmi.bitfield_bitmask,
						       nsmi.buffer_bitvalue);
			}
		}
	else	{
		// it's "outside" ==> no-op here
		}
	}
	}
	}

if (verbose_info.print_algorithm_details)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      %ld \"inside\" points, %ld \"buffer\" points on this processor",
		   inside_count, buffer_count);
}
	  }

//******************************************************************************

//
// This function computes the inner mask radius  r_inner  for a given
// horizon radius  r_horizon .
//
namespace {
fp inner_mask_radius(const struct mask_grid_info& mgi,
		     const struct mask_info& mask_info,
		     fp r_horizon)
{
const fp Cactus_dx = mgi.base_grid_mean_coord_delta;
return mask_info.radius_multiplier*r_horizon
       + mask_info.radius_offset*Cactus_dx;
}
	  }

//******************************************************************************

//
// This function computes the outer mask radius  r_outer  for a given
// inner mask radius  r_inner .
//
namespace {
fp outer_mask_radius(const struct mask_grid_info& mgi,
		     const struct mask_info& mask_info,
		     fp r_inner)
{
const fp Cactus_dx = mgi.base_grid_mean_coord_delta;
return r_inner + mask_info.buffer_thickness*Cactus_dx;
}
	  }

//******************************************************************************

	  }	// namespace AHFinderDirect
