// setup.cc -- top level driver to setup our persistent data structures
// $Header$
//
// <<<prototypes for functions local to this file>>>
// <<<access to persistent data>>>
//
// AHFinderDirect_setup - top-level driver to setup persistent data structures
///
/// decode_method - decode the  method  parameter
/// decode_verbose_level - decode the  verbose_level  parameter
///
/// allocate_horizons_to_processor - choose which horizons this proc will find
///
/// choose_patch_system_type - choose patch system type
///

#include <assert.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include <vector>

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

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// ***** prototypes for functions local to this file *****
//

namespace {
enum method
  decode_method(const char method_string[]);
enum verbose_level
  decode_verbose_level(const char verbose_level_string[]);

int allocate_horizons_to_processor(int N_procs, int my_proc,
				   int N_horizons, bool multiproc_flag,
                                   const CCTK_INT depends_on[],
				   horizon_sequence& my_hs,
				   const struct verbose_info& verbose_info);

enum patch_system::patch_system_type
  choose_patch_system_type(const char grid_domain[],
			   const char grid_bitant_plane[],
			   const char grid_quadrant_direction[],
			   const char grid_rotation_axis[],
			   fp origin_x, fp origin_y, fp origin_z);
	  }

void set_initial_guess_parameters(struct AH_data& AH_data, const int hn,
                                  const fp ini_origin_x, const fp ini_origin_y, const fp ini_origin_z);





//******************************************************************************

//
// This function is used to setup initial guesses either from 
// the corresponding parameters, or from a custom ini_origin_*.
// The latter is necessary when tracking with a grid scalar.
//
//

void set_initial_guess_parameters(struct AH_data& AH_data, const int hn,
                                  const fp ini_origin_x, const fp ini_origin_y, const fp ini_origin_z)
{
DECLARE_CCTK_PARAMETERS;

                AH_data.initial_guess_info.method
			= decode_initial_guess_method(initial_guess_method[hn]);
                AH_data.initial_guess_info.reset_horizon_after_not_finding
                        = reset_horizon_after_not_finding[hn];
		// ... read from named file
		AH_data.initial_guess_info.read_from_named_file_info.file_name
			= initial_guess__read_from_named_file__file_name[hn];

if (!track_origin_from_grid_scalar[hn]) {

		// ... Kerr/Kerr
		AH_data.initial_guess_info.Kerr_Kerr_info.x_posn
			= initial_guess__Kerr_Kerr__x_posn[hn];
		AH_data.initial_guess_info.Kerr_Kerr_info.y_posn
			= initial_guess__Kerr_Kerr__y_posn[hn];
		AH_data.initial_guess_info.Kerr_Kerr_info.z_posn
			= initial_guess__Kerr_Kerr__z_posn[hn];
		AH_data.initial_guess_info.Kerr_Kerr_info.mass
			= initial_guess__Kerr_Kerr__mass[hn];
		AH_data.initial_guess_info.Kerr_Kerr_info.spin
			= initial_guess__Kerr_Kerr__spin[hn];
		// ... Kerr/Kerr-Schild
		AH_data.initial_guess_info.Kerr_KerrSchild_info.x_posn
			= initial_guess__Kerr_KerrSchild__x_posn[hn];
		AH_data.initial_guess_info.Kerr_KerrSchild_info.y_posn
			= initial_guess__Kerr_KerrSchild__y_posn[hn];
		AH_data.initial_guess_info.Kerr_KerrSchild_info.z_posn
			= initial_guess__Kerr_KerrSchild__z_posn[hn];
		AH_data.initial_guess_info.Kerr_KerrSchild_info.mass
			= initial_guess__Kerr_KerrSchild__mass[hn];
		AH_data.initial_guess_info.Kerr_KerrSchild_info.spin
			= initial_guess__Kerr_KerrSchild__spin[hn];
		// ... coordinate sphere
		AH_data.initial_guess_info.coord_sphere_info.x_center
			= initial_guess__coord_sphere__x_center[hn];
		AH_data.initial_guess_info.coord_sphere_info.y_center
			= initial_guess__coord_sphere__y_center[hn];
		AH_data.initial_guess_info.coord_sphere_info.z_center
			= initial_guess__coord_sphere__z_center[hn];
		AH_data.initial_guess_info.coord_sphere_info.radius
			= initial_guess__coord_sphere__radius[hn];
		// ... coordinate ellipsoid
		AH_data.initial_guess_info.coord_ellipsoid_info.x_center
			= initial_guess__coord_ellipsoid__x_center[hn];
		AH_data.initial_guess_info.coord_ellipsoid_info.y_center
			= initial_guess__coord_ellipsoid__y_center[hn];
		AH_data.initial_guess_info.coord_ellipsoid_info.z_center
			= initial_guess__coord_ellipsoid__z_center[hn];
		AH_data.initial_guess_info.coord_ellipsoid_info.x_radius
			= initial_guess__coord_ellipsoid__x_radius[hn];
		AH_data.initial_guess_info.coord_ellipsoid_info.y_radius
			= initial_guess__coord_ellipsoid__y_radius[hn];
		AH_data.initial_guess_info.coord_ellipsoid_info.z_radius
			= initial_guess__coord_ellipsoid__z_radius[hn];

} else {

		// ... Kerr/Kerr
		AH_data.initial_guess_info.Kerr_Kerr_info.x_posn
			= ini_origin_x;
		AH_data.initial_guess_info.Kerr_Kerr_info.y_posn
			= ini_origin_y;
		AH_data.initial_guess_info.Kerr_Kerr_info.z_posn
			= ini_origin_z;
		AH_data.initial_guess_info.Kerr_Kerr_info.mass
			= initial_guess__Kerr_Kerr__mass[hn];
		AH_data.initial_guess_info.Kerr_Kerr_info.spin
			= initial_guess__Kerr_Kerr__spin[hn];
		// ... Kerr/Kerr-Schild
		AH_data.initial_guess_info.Kerr_KerrSchild_info.x_posn
			= ini_origin_x;
		AH_data.initial_guess_info.Kerr_KerrSchild_info.y_posn
			= ini_origin_y;
		AH_data.initial_guess_info.Kerr_KerrSchild_info.z_posn
			= ini_origin_z;
		AH_data.initial_guess_info.Kerr_KerrSchild_info.mass
			= initial_guess__Kerr_KerrSchild__mass[hn];
		AH_data.initial_guess_info.Kerr_KerrSchild_info.spin
			= initial_guess__Kerr_KerrSchild__spin[hn];
		// ... coordinate sphere
		AH_data.initial_guess_info.coord_sphere_info.x_center
			= ini_origin_x;
		AH_data.initial_guess_info.coord_sphere_info.y_center
			= ini_origin_y;
		AH_data.initial_guess_info.coord_sphere_info.z_center
			= ini_origin_z;
		AH_data.initial_guess_info.coord_sphere_info.radius
			= initial_guess__coord_sphere__radius[hn];
		// ... coordinate ellipsoid
		AH_data.initial_guess_info.coord_ellipsoid_info.x_center
			= ini_origin_x;
		AH_data.initial_guess_info.coord_ellipsoid_info.y_center
			= ini_origin_y;
		AH_data.initial_guess_info.coord_ellipsoid_info.z_center
			= ini_origin_z;
		AH_data.initial_guess_info.coord_ellipsoid_info.x_radius
			= initial_guess__coord_ellipsoid__x_radius[hn];
		AH_data.initial_guess_info.coord_ellipsoid_info.y_radius
			= initial_guess__coord_ellipsoid__y_radius[hn];
		AH_data.initial_guess_info.coord_ellipsoid_info.z_radius
			= initial_guess__coord_ellipsoid__z_radius[hn];

}

}

//******************************************************************************

//
// ***** access to persistent data *****
//
extern struct state state;

//******************************************************************************

//
// This function is called by the Cactus scheduler to set up all our
// persistent data structures.  (These are stored in  struct state .)
//
extern "C"
  void AHFinderDirect_setup(CCTK_ARGUMENTS)
{
DECLARE_CCTK_ARGUMENTS_AHFinderDirect_setup
DECLARE_CCTK_PARAMETERS

CCTK_VInfo(CCTK_THORNSTRING,
	   "setting up AHFinderDirect data structures");


//
// check parameters
//
int need_zones = 0;
for (int n=1; n<N_horizons; ++n) {
  need_zones = jtutil::max(need_zones, int(N_zones_per_right_angle[n]));
}
if (need_zones > max_N_zones_per_right_angle) {
  CCTK_VWarn (FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
              "AHFinderDirect_setup(): "
              "The parameter max_N_zones_per_right_angle must be at least the maximum of all N_zones_per_right_angle[] parameters.  "
              "Set max_N_zones_per_right_angle to %d or higher to continue.",
              need_zones);      /*NOTREACHED*/
 }
for (int n=1; n<N_horizons; ++n) {
  if (depends_on[n] != 0) {
    assert (depends_on[n] >= 1 && depends_on[n] < n);
    if (N_zones_per_right_angle[n] != N_zones_per_right_angle[depends_on[n]]) {
      CCTK_VWarn (FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "AHFinderDirect_setup(): "
                  "The parameter N_zones_per_right_angle must be the same for a horizon as for the horizon on which it depends.  "
                  "Horizon %d depends on horizon %d, but they have different resolutions.",
                  n, int(depends_on[n])); /*NOTREACHED*/
    }
  }
}

bool find_individual_is_set = false;
for (int n = 1 ; n < N_horizons ; ++n)
{
if (find_every_individual[n] > 0)
   then find_individual_is_set = true;
}
if (move_origins && find_individual_is_set)
   then {
	CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
"Find_every_individual is currently not compatible with moving origins via the move_origins parameter.  "
"Make sure simultaneous moving horizons are not being calculated at different frequencies." );
	}

//
// basic setup
//
state.method = decode_method(method);

state.error_info.warn_level__point_outside__initial
   = warn_level__point_outside__initial;
state.error_info.warn_level__point_outside__subsequent
   = warn_level__point_outside__subsequent;
state.error_info.warn_level__skipping_finite_check
   = warn_level__skipping_finite_check;
state.error_info.warn_level__nonfinite_geometry
   = warn_level__nonfinite_geometry;
state.error_info.warn_level__gij_not_positive_definite__initial
   = warn_level__gij_not_positive_definite__initial;
state.error_info.warn_level__gij_not_positive_definite__subsequent
   = warn_level__gij_not_positive_definite__subsequent;

struct verbose_info& verbose_info = state.verbose_info;
verbose_info.verbose_level = decode_verbose_level(verbose_level);
verbose_info.print_physics_highlights
   = (state.verbose_info.verbose_level >= verbose_level__physics_highlights);
verbose_info.print_physics_details
   = (state.verbose_info.verbose_level >= verbose_level__physics_details);
verbose_info.print_algorithm_highlights
   = (state.verbose_info.verbose_level >= verbose_level__algorithm_highlights);
verbose_info.print_algorithm_details
   = (state.verbose_info.verbose_level >= verbose_level__algorithm_details);
verbose_info.print_algorithm_debug
   = (state.verbose_info.verbose_level >= verbose_level__algorithm_debug);

state.timer_handle = (print_timing_stats != 0) ? CCTK_TimerCreate("finding apparent horizons") : -1;

state.N_procs = CCTK_nProcs(cctkGH);
state.my_proc = CCTK_MyProc(cctkGH);

state.N_horizons = N_horizons;
state.N_active_procs = 0;	// dummy value, will be set properly later
CCTK_VInfo(CCTK_THORNSTRING,
           "           to search for %d horizon%s on %d processor%s",
	   state.N_horizons, ((state.N_horizons == 1) ? "" : "s"),
	   state.N_procs, ((state.N_procs == 1) ? "" : "s"));


//
// Cactus grid info
//
if (verbose_info.print_algorithm_highlights)
   then CCTK_VInfo(CCTK_THORNSTRING, "   setting up Cactus grid info");
struct cactus_grid_info& cgi = state.cgi;
cgi.GH = cctkGH;
cgi.coord_system_handle = CCTK_CoordSystemHandle(coordinate_system_name);
if (cgi.coord_system_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_setup(): can't find Cactus coordinate system \"%s\"!",
		   coordinate_system_name);			/*NOTREACHED*/
cgi.use_Cactus_conformal_metric = false;	// dummy value, may change later

cgi.mask_varindex    = Cactus_gridfn_varindex("AHFinderDirect::ahmask");
cgi.g_dd_11_varindex = Cactus_gridfn_varindex("ADMBase::gxx");
cgi.g_dd_12_varindex = Cactus_gridfn_varindex("ADMBase::gxy");
cgi.g_dd_13_varindex = Cactus_gridfn_varindex("ADMBase::gxz");
cgi.g_dd_22_varindex = Cactus_gridfn_varindex("ADMBase::gyy");
cgi.g_dd_23_varindex = Cactus_gridfn_varindex("ADMBase::gyz");
cgi.g_dd_33_varindex = Cactus_gridfn_varindex("ADMBase::gzz");
cgi.K_dd_11_varindex = Cactus_gridfn_varindex("ADMBase::kxx");
cgi.K_dd_12_varindex = Cactus_gridfn_varindex("ADMBase::kxy");
cgi.K_dd_13_varindex = Cactus_gridfn_varindex("ADMBase::kxz");
cgi.K_dd_22_varindex = Cactus_gridfn_varindex("ADMBase::kyy");
cgi.K_dd_23_varindex = Cactus_gridfn_varindex("ADMBase::kyz");
cgi.K_dd_33_varindex = Cactus_gridfn_varindex("ADMBase::kzz");
cgi.psi_varindex     = Cactus_gridfn_varindex("StaticConformal::psi");


//
// geometry info
//
if (verbose_info.print_algorithm_highlights)
   then CCTK_VInfo(CCTK_THORNSTRING, "   setting up geometry interpolator");
struct geometry_info& gi = state.gi;
gi.hardwire_Schwarzschild_EF_geometry
	= (hardwire_Schwarzschild_EF_geometry != 0);

gi.operator_handle = CCTK_InterpHandle(geometry_interpolator_name);
if (gi.operator_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_setup(): couldn't find interpolator \"%s\"!",
		   geometry_interpolator_name);		/*NOTREACHED*/
gi.param_table_handle = Util_TableCreateFromString(geometry_interpolator_pars);
if (gi.param_table_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_setup(): bad geometry-interpolator parameter(s) \"%s\"!",
		   geometry_interpolator_pars);		/*NOTREACHED*/

gi.geometry__Schwarzschild_EF__mass     = geometry__Schwarzschild_EF__mass;
gi.geometry__Schwarzschild_EF__x_posn   = geometry__Schwarzschild_EF__x_posn;
gi.geometry__Schwarzschild_EF__y_posn   = geometry__Schwarzschild_EF__y_posn;
gi.geometry__Schwarzschild_EF__z_posn   = geometry__Schwarzschild_EF__z_posn;
gi.geometry__Schwarzschild_EF__epsilon  = geometry__Schwarzschild_EF__epsilon;
gi.geometry__Schwarzschild_EF__Delta_xyz= geometry__Schwarzschild_EF__Delta_xyz;

gi.check_that_h_is_finite        = (check_that_h_is_finite        != 0);
gi.check_that_geometry_is_finite = (check_that_geometry_is_finite != 0);
gi.mask_is_noshrink		 = (mask_is_noshrink != 0);

//
// Jacobian info
//
struct Jacobian_info& Jac_info = state.Jac_info;
Jac_info.Jacobian_compute_method
	= decode_Jacobian_compute_method(Jacobian_compute_method);
Jac_info.Jacobian_store_solve_method
	= decode_Jacobian_store_solve_method(Jacobian_store_solve_method);
Jac_info.perturbation_amplitude = Jacobian_perturbation_amplitude;


//
// solver info
//
struct solver_info& solver_info = state.solver_info;
solver_info.debugging_output_at_each_Newton_iteration
			= (debugging_output_at_each_Newton_iteration != 0);
solver_info.linear_solver_pars.ILUCG_pars.error_tolerance
	= ILUCG__error_tolerance;
solver_info.linear_solver_pars.ILUCG_pars.limit_CG_iterations
	= (ILUCG__limit_CG_iterations != 0);
solver_info.linear_solver_pars.UMFPACK_pars.N_II_iterations
	= UMFPACK__N_II_iterations;
solver_info.max_Newton_iterations__initial
	= max_Newton_iterations__initial;
solver_info.max_Newton_iterations__subsequent
	= max_Newton_iterations__subsequent;
solver_info.max_allowable_Delta_h_over_h = max_allowable_Delta_h_over_h;
solver_info.Theta_norm_for_convergence   = Theta_norm_for_convergence;
solver_info.max_allowable_Theta          = max_allowable_Theta;
solver_info.max_allowable_Theta_growth_iterations
                                       = max_allowable_Theta_growth_iterations;
solver_info.max_allowable_Theta_nonshrink_iterations
                                    = max_allowable_Theta_nonshrink_iterations;
// ... horizon numbers run from 1 to N_horizons inclusive
//     so the array size is N_horizons+1
solver_info.max_allowable_horizon_radius = new fp[state.N_horizons+1];
	  {
	for (int hn = 0 ; hn <= N_horizons ; ++hn)
	{
	solver_info.max_allowable_horizon_radius[hn]
		= max_allowable_horizon_radius[hn];
	}
	  }
solver_info.want_expansion_gradients = want_expansion_gradients;


//
// I/O info
//
struct IO_info& IO_info = state.IO_info;
IO_info.output_ASCII_files = (output_ASCII_files != 0);
IO_info.output_HDF5_files = (output_HDF5_files != 0);
IO_info.output_initial_guess = (output_initial_guess != 0);
IO_info.output_h_every     = output_h_every;
IO_info.output_Theta_every = output_Theta_every;
IO_info.output_mean_curvature_every = output_mean_curvature_every;
IO_info.output_h     = false;	// dummy value
IO_info.output_Theta = false;	// dummy value
IO_info.output_mean_curvature = false;	// dummy value

IO_info.output_BH_diagnostics              = (output_BH_diagnostics != 0);
IO_info.BH_diagnostics_directory
	= (strlen(BH_diagnostics_directory) == 0)
	  ? /* IO:: */ out_dir
	  : BH_diagnostics_directory;
IO_info.BH_diagnostics_base_file_name      = BH_diagnostics_base_file_name;
IO_info.BH_diagnostics_file_name_extension = BH_diagnostics_file_name_extension;

IO_info.output_ghost_zones_for_h  = (output_ghost_zones_for_h != 0);
IO_info.ASCII_gnuplot_file_name_extension = ASCII_gnuplot_file_name_extension;
IO_info.HDF5_file_name_extension          = HDF5_file_name_extension;
IO_info.h_directory
	= (strlen(h_directory) == 0)
	  ? /* IO:: */ out_dir
	  : h_directory;
IO_info.h_base_file_name         = h_base_file_name;
IO_info.Theta_base_file_name     = Theta_base_file_name;
IO_info.mean_curvature_base_file_name     = mean_curvature_base_file_name;
IO_info.Delta_h_base_file_name   = Delta_h_base_file_name;
IO_info.h_min_digits             = h_min_digits;
IO_info.Jacobian_base_file_name  = Jacobian_base_file_name;
IO_info.output_OpenDX_control_files  = (output_OpenDX_control_files != 0);
IO_info.OpenDX_control_file_name_extension = OpenDX_control_file_name_extension;
IO_info.time_iteration = 0;
IO_info.time           = 0.0;


//
// other misc setup
//
state.BH_diagnostics_info.integral_method
   = patch::decode_integration_method(integral_method);


//
// mask parameters
//
struct mask_info& mask_info = state.mask_info;
mask_info.set_mask_for_any_horizon = false;
// ... horizon numbers run from 1 to N_horizons inclusive
//     so the array size is N_horizons+1
mask_info.set_mask_for_this_horizon = new bool[N_horizons+1];
	  {
	for (int hn = 1 ; hn <= N_horizons ; ++hn)
	{
	mask_info.set_mask_for_this_horizon[hn]
		= (set_mask_for_all_horizons != 0)
		  || (set_mask_for_individual_horizon[hn] != 0);
	mask_info.set_mask_for_any_horizon
		|= mask_info.set_mask_for_this_horizon[hn];
	}
	  }
if (mask_info.set_mask_for_any_horizon)
   then {
	mask_info.radius_multiplier  = mask_radius_multiplier;
	mask_info.radius_offset      = mask_radius_offset;
	mask_info.buffer_thickness   = mask_buffer_thickness;
	mask_info.mask_is_noshrink   = mask_is_noshrink;
	mask_info.min_horizon_radius_points_for_mask
				     = min_horizon_radius_points_for_mask;
	mask_info.set_old_style_mask = (set_old_style_mask != 0);
	mask_info.set_new_style_mask = (set_new_style_mask != 0);
	if (mask_info.set_old_style_mask)
	   then {
		struct mask_info::old_style_mask_info& osmi
			= mask_info.old_style_mask_info;
		osmi.gridfn_name     = old_style_mask_gridfn_name;
		osmi.gridfn_varindex = Cactus_gridfn_varindex(osmi.gridfn_name);
		osmi.gridfn_dataptr  = NULL;	// dummy value; fixup later
		osmi.inside_value  = old_style_mask_inside_value;
		osmi.buffer_value  = old_style_mask_buffer_value;
		osmi.outside_value = old_style_mask_outside_value;
		}
	if (mask_info.set_new_style_mask)
	   then {
		struct mask_info::new_style_mask_info& nsmi
			= mask_info.new_style_mask_info;
		nsmi.gridfn_name     = new_style_mask_gridfn_name;
		nsmi.gridfn_varindex = Cactus_gridfn_varindex(nsmi.gridfn_name);
		nsmi.gridfn_dataptr  = NULL;	// dummy value; fixup later
		nsmi.bitfield_name   = new_style_mask_bitfield_name;
		nsmi.bitfield_bitmask = 0;	// dummy value; fixup later
		nsmi.inside_value    = new_style_mask_inside_value;
		nsmi.buffer_value    = new_style_mask_buffer_value;
		nsmi.outside_value   = new_style_mask_outside_value;
		nsmi.inside_bitvalue = 0;	// dummy value; fixup later
		nsmi.buffer_bitvalue = 0;	// dummy value; fixup later
		nsmi.outside_bitvalue = 0;	// dummy value; fixup later
		}
	}


//
// (genuine) horizon sequence for this processor
//
state.my_hs = new horizon_sequence(state.N_horizons);
horizon_sequence& hs = *state.my_hs;

//
// if we're going to actually find horizons
//    we spread the horizons over multiple processors for maximum efficiency,
// otherwise (we're just doing testing/debugging computations, so)
//    we allocate all the horizons to processor #0 for simplicity
//
const bool multiproc_flag = (state.method == method__find_horizons);
state.N_active_procs
   = allocate_horizons_to_processor(state.N_procs, state.my_proc,
				    state.N_horizons, multiproc_flag,
                                    depends_on,
				    hs,
				    verbose_info);

// ... horizon numbers run from 1 to N_horizons inclusive
//     so the array size is N_horizons+1
state.AH_data_array = new AH_data*[N_horizons+1];
	  {
	for (int hn = 0 ; hn <= N_horizons ; ++hn)
	{
	state.AH_data_array[hn] = NULL;
	}
	  }


//
// horizon-specific info for each horizon
//

// set up the interpatch interpolator
if (verbose_info.print_algorithm_highlights)
   then CCTK_VInfo(CCTK_THORNSTRING, "   setting up interpatch interpolator");
const int ip_interp_handle = CCTK_InterpHandle(interpatch_interpolator_name);
if (ip_interp_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_setup(): couldn't find interpatch interpolator \"%s\"!",
		   interpatch_interpolator_name);		/*NOTREACHED*/
const int ip_interp_param_table_handle
	= Util_TableCreateFromString(interpatch_interpolator_pars);
if (ip_interp_param_table_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_setup(): bad interpatch-interpolator parameter(s) \"%s\"!",
		   interpatch_interpolator_pars);		/*NOTREACHED*/

// set up the surface interpolator if it's going to be used
int surface_interp_handle = -1;
int surface_interp_param_table_handle = -1;
if (strlen(surface_interpolator_name) > 0)
   then {
	if (verbose_info.print_algorithm_highlights)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "   setting up surface interpolator");
	surface_interp_handle = CCTK_InterpHandle(surface_interpolator_name);
	if (surface_interp_handle < 0)
	   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_setup(): couldn't find surface interpolator \"%s\"!",
			   surface_interpolator_name);		/*NOTREACHED*/
	surface_interp_param_table_handle
		= Util_TableCreateFromString(surface_interpolator_pars);
	if (surface_interp_param_table_handle < 0)
	   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_setup(): bad surface-interpolator parameter(s) \"%s\"!",
			   surface_interpolator_pars);		/*NOTREACHED*/
	}

// setup all horizons on this processor,
// with full-fledged patch systems for genuine horizons
// and skeletal patch systems for others
	  {
	for (int hn = 1 ; hn <= hs.N_horizons() ; ++hn)
	{
	const bool genuine_flag = hs.is_hn_genuine(hn);
	state.AH_data_array[hn] = new AH_data;
	struct AH_data& AH_data = *state.AH_data_array[hn];

	if (verbose_info.print_algorithm_highlights)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "   setting up %s data structures for horizon %d",
			   (genuine_flag ? "full-fledged" : "skeletal"),
			   hn);

	// decide what type of patch system this one should be
	const enum patch_system::patch_system_type ps_type
	   = STRING_EQUAL(patch_system_type[hn], "match Cactus grid symmetry")
	     ? // choose a patch system type based on grid:: parameters
	       // ... we inherit from grid, and we ask for some of its
	       //     parameters in our param.ccl file; these appear as
	       //     ordinary Cactus parameters here, so (eg)
	       //     grid::domain is just "domain" here
	       choose_patch_system_type(/* grid:: */ domain,
					/* grid:: */ bitant_plane,
					/* grid:: */ quadrant_direction,
					/* grid:: */ rotation_axis,
					origin_x[hn],
					origin_y[hn],
					origin_z[hn])
	     : patch_system::type_of_name(patch_system_type[hn]);

	// create the patch system
	AH_data.ps_ptr
	   = new patch_system(origin_x[hn], origin_y[hn], origin_z[hn],
			      ps_type,
			      ghost_zone_width, patch_overlap_width,
			      N_zones_per_right_angle[hn],
			      gfns::nominal_min_gfn,
			      (genuine_flag ? gfns::nominal_max_gfn
					    : gfns::skeletal_nominal_max_gfn),
			      gfns::ghosted_min_gfn, gfns::ghosted_max_gfn,
			      ip_interp_handle, ip_interp_param_table_handle,
			      surface_interp_handle,
			      surface_interp_param_table_handle,
			      true, verbose_info.print_algorithm_details);
	patch_system& ps = *AH_data.ps_ptr;
	if (genuine_flag)
	   then ps.set_gridfn_to_constant(0.0, gfns::gfn__zero);
	if (genuine_flag)
	   then ps.set_gridfn_to_constant(1.0, gfns::gfn__one);

	AH_data.Jac_ptr = genuine_flag
			  ? new_Jacobian(Jac_info.Jacobian_store_solve_method,
					 ps,
					 verbose_info.print_algorithm_details)
			  : NULL;

        AH_data.compute_info.surface_definition =
          STRING_EQUAL(surface_definition[hn], "expansion")
          ? definition_expansion
          : STRING_EQUAL(surface_definition[hn], "inner expansion")
          ? definition_inner_expansion
          : STRING_EQUAL(surface_definition[hn], "mean curvature")
          ? definition_mean_curvature
          : STRING_EQUAL(surface_definition[hn], "expansion product")
          ? definition_expansion_product
          : (CCTK_WARN (0, "internal error"), definition_error);
        AH_data.compute_info.surface_modification =
          STRING_EQUAL(surface_modification[hn], "none")
          ? modification_none
          : STRING_EQUAL(surface_modification[hn], "radius")
          ? modification_radius
          : STRING_EQUAL(surface_modification[hn], "radius^2")
          ? modification_radius2
#if 0
          : STRING_EQUAL(surface_modification[hn], "mean radius")
          ? modification_mean_radius
          : STRING_EQUAL(surface_modification[hn], "areal radius")
          ? modification_areal_radius
#endif
          : (CCTK_WARN (0, "internal error"), modification_error);
        AH_data.compute_info.surface_selection = 
          STRING_EQUAL(surface_selection[hn], "definition")
          ? selection_definition
          : STRING_EQUAL(surface_selection[hn], "mean coordinate radius")
          ? selection_mean_coordinate_radius
          : STRING_EQUAL(surface_selection[hn], "areal radius")
          ? selection_areal_radius
          : STRING_EQUAL(surface_selection[hn], "expansion times mean coordinate radius")
          ? selection_expansion_mean_coordinate_radius
          : STRING_EQUAL(surface_selection[hn], "expansion times areal radius")
          ? selection_expansion_areal_radius
          : (CCTK_WARN (0, "internal error"), selection_error);
	AH_data.compute_info.desired_value = desired_value[hn];

        AH_data.move_origins               = move_origins;

        AH_data.use_pretracking            = use_pretracking[hn];
        AH_data.pretracking_max_iterations = pretracking_max_iterations[hn];

        AH_data.pretracking_value          = pretracking_value[hn];
        AH_data.pretracking_minimum_value  = pretracking_minimum_value[hn];
        AH_data.pretracking_maximum_value  = pretracking_maximum_value[hn];
        AH_data.pretracking_delta          = pretracking_delta[hn];
        AH_data.pretracking_minimum_delta  = pretracking_minimum_delta[hn];
        AH_data.pretracking_maximum_delta  = pretracking_maximum_delta[hn];

        AH_data.depends_on           = depends_on[hn];
        AH_data.desired_value_factor = desired_value_factor[hn];
        AH_data.desired_value_offset = desired_value_offset[hn];

        AH_data.shiftout_factor = shiftout_factor[hn];
        AH_data.smoothing_factor = smoothing_factor[hn];

	// AH_data.initial_find_flag = genuine_flag;
	// AH_data.really_initial_find_flag = AH_data.initial_find_flag;
	AH_data.initial_find_flag = true;
	AH_data.really_initial_find_flag = AH_data.initial_find_flag;

	if (genuine_flag)
	   then {
		if (verbose_info.print_algorithm_details)
		   then CCTK_VInfo(CCTK_THORNSTRING,
			   "      setting initial guess parameters etc");
		set_initial_guess_parameters(AH_data, hn, /* irrelevant here; leave at zero */0, 0, 0);
		}

	AH_data.search_flag = false;
	AH_data.found_flag = false;
	AH_data.h_files_written = false;
	AH_data.BH_diagnostics_fileptr = NULL;
	}
	  }

  // Initialise the centroids.  These values may later be overwritten
  // when they are recovered from a checkpoint.  However, if new
  // horizons are enabled during recovery, they are then correctly
  // initialised.
  for (int n = 0; n < N_horizons; ++ n) {
    // Horizon centroids are not valid
    ah_centroid_x[n] = 0.0;
    ah_centroid_y[n] = 0.0;
    ah_centroid_z[n] = 0.0;
    ah_centroid_t[n] = 0.0;
    ah_centroid_valid[n] = 0;
    ah_centroid_iteration[n] = -1;
    
    ah_centroid_x_p[n] = 0.0;
    ah_centroid_y_p[n] = 0.0;
    ah_centroid_z_p[n] = 0.0;
    ah_centroid_t_p[n] = 0.0;
    ah_centroid_valid_p[n] = 0;
    ah_centroid_iteration_p[n] = -1;
    
    struct AH_data& AH_data = *state.AH_data_array[n+1];
    ah_initial_find_flag[n]        = AH_data.initial_find_flag;
    ah_really_initial_find_flag[n] = AH_data.really_initial_find_flag;
    ah_search_flag[n]              = AH_data.search_flag;
    ah_found_flag[n]               = AH_data.found_flag;
    if (verbose_info.print_algorithm_details) {
      printf ("AHF setup %d initial_find_flag=%d\n",        n+1, (int) AH_data.initial_find_flag);
      printf ("AHF setup %d really_initial_find_flag=%d\n", n+1, (int) AH_data.really_initial_find_flag);
      printf ("AHF setup %d search_flag=%d\n",              n+1, (int) AH_data.search_flag);
      printf ("AHF setup %d found_flag=%d\n",               n+1, (int) AH_data.found_flag);
    }
  }

  // Save in grid array in case a recovery only recovers some horizons
  for (int n = 0; n < N_horizons; ++ n) {
    struct AH_data& AH_data = *state.AH_data_array[n+1];
    const struct BH_diagnostics& BH_diagnostics = AH_data.BH_diagnostics;
    BH_diagnostics.save(cctkGH, n+1);
  }
}


//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function decodes the  method  parameter (string) into an
// internal enum for future use.
//
namespace {
enum method
  decode_method(const char method_string[])
{
if	(STRING_EQUAL(method_string, "evaluate expansions"))
   then return method__evaluate_expansions;
else if (STRING_EQUAL(method_string, "test expansion Jacobians"))
   then return method__test_expansion_Jacobians;
else if (STRING_EQUAL(method_string, "find horizons"))
   then return method__find_horizons;
else	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
		    "decode_method(): unknown method_string=\"%s\"!",
		    method_string);				/*NOTREACHED*/
}
	  }

//******************************************************************************

//
// This function decodes the  verbose_level  parameter (string) into an
// internal enum for future use.
//
namespace {
enum verbose_level
  decode_verbose_level(const char verbose_level_string[])
{
if	(STRING_EQUAL(verbose_level_string, "physics highlights"))
   then return verbose_level__physics_highlights;
else if (STRING_EQUAL(verbose_level_string, "physics details"))
   then return verbose_level__physics_details;
else if (STRING_EQUAL(verbose_level_string, "algorithm highlights"))
   then return verbose_level__algorithm_highlights;
else if (STRING_EQUAL(verbose_level_string, "algorithm details"))
   then return verbose_level__algorithm_details;
else if (STRING_EQUAL(verbose_level_string, "algorithm debug"))
   then return verbose_level__algorithm_debug;
else	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
"decode_verbose_level(): unknown verbose_level_string=\"%s\"!",
		   verbose_level_string);			/*NOTREACHED*/
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function allocates horizons to a processor.  That is, it chooses
// which horizons a given processor (in practice, the current processor)
// will attempt to find.  It records this choice in a  horizon_sequence
// object.
//
// Arguments:
// N_procs = The total number of processors.
// my_proc = My processor number (0-origin).
// N_horizons = The total number of (genuine) horizons to be found.
// multiproc_flag = true to spread the horizons over multiple processors,
//		    false to allocate all horizons to processor #0.
// my_hs = (out) This will be set to the sequence of (genuine) horizons
//		 allocated to the specified processor.
// verbose_info = Specifies what information to print about the allocation.
//
// Results:
// This function returns the total number of active processors, i.e. the
// total number of processors which are allocated one or more genuine
// horizons allocated.
//
namespace {
int allocate_horizons_to_processor(int N_procs, int my_proc,
				   int N_horizons, bool multiproc_flag,
                                   const CCTK_INT depends_on[],
				   horizon_sequence& my_hs,
				   const struct verbose_info& verbose_info)
{
const int N_active_procs = multiproc_flag ? jtutil::min(N_procs, N_horizons)
					  : 1;

//
// Implementation note:
// We allocate the horizons to active processors in round-robin order.
//
std::vector<int> proc_of_horizon (N_horizons+1);
	for (int hn = 1 ; hn <= N_horizons ; ++hn)
        {
        proc_of_horizon.at(hn) = -1;
        }

int proc = 0;
	for (int hn = 1 ; hn <= N_horizons ; ++hn)
	{
        if (depends_on[hn] < 0 || depends_on[hn] > N_horizons) {
          CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "horizon %d depends on a horizon with the illegal index %d",
                     hn, int(depends_on[hn]));
        } else if (depends_on[hn] == hn) {
          CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "horizon %d depends on itself",
                     hn);
        } else if (depends_on[hn] > hn) {
          CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "horizon %d depends on a horizon with a larger index %d",
                     hn, int(depends_on[hn]));
        }
        const int this_horizons_proc
          = depends_on[hn] == 0 ? proc : proc_of_horizon.at(depends_on[hn]);
        assert (this_horizons_proc >= 0 && this_horizons_proc < N_procs);
        proc_of_horizon.at(hn) = this_horizons_proc;
	if (verbose_info.print_algorithm_highlights)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "   allocating horizon %d to processor #%d",
			   hn, this_horizons_proc);
	if (this_horizons_proc == my_proc)
           then my_hs.append_hn(hn);
	if (++proc >= N_active_procs)
	   then proc = 0;
	}

return N_active_procs;
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function chooses a patch system type based on the Cactus grid
// symmetry and the patch system's origin position.
//
// Arguments:
// grid_* = The values of the Cactus parameters
//		grid::domain
//		grid::bitant_plane
//		grid::quadrant_direction
//		grid::rotation_axis
// origin_[xyz] = The origin position for this patch system.
//
namespace {
enum patch_system::patch_system_type
  choose_patch_system_type(const char grid_domain[],
			   const char grid_bitant_plane[],
			   const char grid_quadrant_direction[],
			   const char grid_rotation_axis[],
			   fp origin_x, fp origin_y, fp origin_z)
{
if	(STRING_EQUAL(grid_domain, "full"))
   then return patch_system::patch_system__full_sphere;

else if ( STRING_EQUAL(grid_domain, "bitant")
	  && STRING_EQUAL(grid_bitant_plane, "xy") )
   then {
	if (origin_z == 0.0)
	   then return patch_system::patch_system__plus_z_hemisphere;
	   else {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "      Cactus has bitant mode about xy plane");
		CCTK_VInfo(CCTK_THORNSTRING,
			   "      but patch system origin_z=%g != 0",
			   double(origin_z));
		CCTK_VInfo(CCTK_THORNSTRING,
			   "      ==> using \"full sphere\" patch system");
		return patch_system::patch_system__full_sphere;
		}
	}

else if ( STRING_EQUAL(grid_domain, "quadrant")
	  && STRING_EQUAL(grid_quadrant_direction, "z") )
   then {
	if ((origin_x == 0.0) && (origin_y == 0.0))
	   then return patch_system::patch_system__plus_xy_quadrant_mirrored;
	   else {
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      Cactus has quadrant mode about z axis");
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      but patch system origin_(x,y)=(%g,%g) != (0,0)",
		   double(origin_x), double(origin_y));
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      ==> using \"full sphere\" patch system");
		return patch_system::patch_system__full_sphere;
		}
	}

else if ( STRING_EQUAL(grid_domain, "quadrant_reflect_rotate")
	  && STRING_EQUAL(grid_quadrant_direction, "z")
	  && STRING_EQUAL(grid_rotation_axis, "z") )
   then {
	if ((origin_x == 0.0) && (origin_y == 0.0))
	   then return patch_system::patch_system__plus_xy_quadrant_rotating;
	   else {
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      Cactus has rotating quadrant mode about z axis");
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      but patch system origin_(x,y)=(%g,%g) != (0,0)",
		   double(origin_x), double(origin_y));
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      ==> using \"full sphere\" patch system");
		return patch_system::patch_system__full_sphere;
		}
	}

else if (STRING_EQUAL(grid_domain, "octant"))
   then {
	if	((origin_x == 0.0) && (origin_y == 0.0) && (origin_z == 0.0))
	   then return patch_system::patch_system__plus_xyz_octant_mirrored;
	else if ((origin_x == 0.0) && (origin_y == 0.0))
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      Cactus has mirrored octant mode");
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      but patch system origin_z=%g != 0",
		   double(origin_z));
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      ==> using \"+xy quadrant (mirrored)\" patch system");
		return patch_system::patch_system__plus_xy_quadrant_mirrored;
		}
	else if ((origin_x == 0.0) && (origin_z == 0.0))
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      Cactus has mirrored octant mode");
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      but patch system origin_y=%g != 0",
		   double(origin_z));
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      ==> using \"+xz quadrant (mirrored)\" patch system");
		return patch_system::patch_system__plus_xz_quadrant_mirrored;
		}
	else if ((origin_y == 0.0) && (origin_z == 0.0))
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      Cactus has mirrored octant mode");
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      but patch system origin_x=%g != 0",
		   double(origin_z));
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      ==> using \"+yz quadrant (mirrored)\" patch system");
		//return patch_system::patch_system__plus_yz_quadrant_mirrored;
                CCTK_WARN (0, "This patch system type is not implemented");
		return patch_system::patch_system__full_sphere;
		}
	else	{
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      Cactus has mirrored octant mode");
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      but patch system origin_(x,y,z)=(%g,%g,%g) != (0,0,0)",
		   double(origin_x), double(origin_y), double(origin_z));
		CCTK_VInfo(CCTK_THORNSTRING,
		   "      ==> using \"full sphere\" patch system");
		return patch_system::patch_system__full_sphere;
		}
	}

else	{
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   choose_patch_system_type()\n"
"        grid::domain = \"%s\" not supported!\n"
"        ==> will try using \"full sphere\" patch system\n"
"            but this may or may not work..."
,
		   grid_domain);
	return patch_system::patch_system__full_sphere;
	}
}
	  }

//******************************************************************************

	  }	// namespace AHFinderDirect
