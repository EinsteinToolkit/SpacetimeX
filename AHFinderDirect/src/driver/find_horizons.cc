// find_horizons.cc -- top level driver for finding apparent horizons
// $Header$
//
// <<<access to persistent data>>>
// <<<prototypes for functions local to this file>>>
// AHFinderDirect_find_horizons - top-level driver to find apparent horizons
///
/// find_horizon - find a horizon
/// do_evaluate_expansions
/// do_test_expansion_Jacobian
///

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

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
// ***** prototypes for functions local to this file
//
namespace {
void do_evaluate_expansions(int my_proc, int N_horizons,
			    horizon_sequence& hs,
			       struct AH_data* const AH_data_array[],
			    const struct cactus_grid_info& cgi,
			    const struct geometry_info& gi,
			    const struct IO_info& IO_info,
			    const struct error_info& error_info,
			    const struct verbose_info& verbose_info,
			    int timer_handle);
void do_test_expansion_Jacobians(int my_proc, int N_horizons,
				 struct AH_data* const AH_data_array[],
				 const struct cactus_grid_info& cgi,
				 const struct geometry_info& gi,
				       struct Jacobian_info& Jac_info,
				 bool test_all_Jacobian_compute_methods,
				 const struct IO_info& IO_info,
				 const struct error_info& error_info,
				 const struct verbose_info& verbose_info,
				 int timer_handle);
	  }

//******************************************************************************

//
// This function is called by the Cactus scheduler to import
// the excision mask.
//
extern "C"
  void AHFinderDirect_import_mask(CCTK_ARGUMENTS)
{
DECLARE_CCTK_ARGUMENTS_AHFinderDirect_import_mask
DECLARE_CCTK_PARAMETERS

assert(ahmask != 0);

for (int k=0; k<cctk_lsh[2]; ++k)
for (int j=0; j<cctk_lsh[1]; ++j)
for (int i=0; i<cctk_lsh[0]; ++i)
{
        const int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);
        // zero means: point can be used,
        // non-zero means: point must be avoided
        ahmask[ind] = 0;
        if (use_mask)
           // grid points with mask values of 1.0 can be used,
           // values of 0.0 and 0.5 must be avoided.
           // the excision boundary cannot be used because
           // (a) it is inaccurate
           // (b) it does not respect the symmetries e.g. in Kerr.
           then ahmask[ind] = fabs(emask[ind] - 1.0) > 0.01;
}
}

//******************************************************************************

//
// This function is called by the Cactus scheduler to find the apparent
// horizon or horizons in the current slice.
//
extern "C"
  void AHFinderDirect_find_horizons(CCTK_ARGUMENTS)
{
DECLARE_CCTK_ARGUMENTS_AHFinderDirect_find_horizons
DECLARE_CCTK_PARAMETERS

// determine whether a horizon should be found at this iteration
bool find_any = false;
for (int hn = 1; hn <= state.my_hs->N_horizons(); ++ hn)
{
  // only try to find horizons every  find_every  time steps
  const int my_find_after = find_after_individual[hn];
  const int my_dont_find_after = dont_find_after_individual[hn];
  const fp my_find_after_time = find_after_individual_time[hn];
  const fp my_dont_find_after_time = dont_find_after_individual_time[hn];
  const int my_find_every = (find_every_individual[hn] >= 0
                             ? find_every_individual[hn]
                             : find_every);
  const bool find_this =    cctk_iteration >= my_find_after
                         && (my_dont_find_after < 0
                             ? true
                             : cctk_iteration <= my_dont_find_after)
                         && cctk_time >= my_find_after_time
                         && (my_dont_find_after_time <= my_find_after_time
                             ? true
                             : cctk_time <= my_dont_find_after_time)
                         && my_find_every > 0
                         && cctk_iteration % my_find_every == 0
                         && ! disable_horizon[hn];
  struct AH_data& AH_data = *state.AH_data_array[hn];
  AH_data.search_flag = find_this;
  find_any = find_any || find_this;
}
if (! find_any) return;

if (state.timer_handle >= 0)
   then CCTK_TimerResetI(state.timer_handle);

const int my_proc = state.my_proc;
horizon_sequence& hs = *state.my_hs;
const bool active_flag = hs.has_genuine_horizons();
const bool broadcast_horizon_shape = true;

      struct cactus_grid_info&          cgi = state.cgi;
const struct    geometry_info&           gi = state.gi;
      struct    Jacobian_info&     Jac_info = state.Jac_info;
      struct          IO_info&      IO_info = state.IO_info;
const struct       error_info&   error_info = state.error_info;
const struct     verbose_info& verbose_info = state.verbose_info;

// what are the semantics of the Cactus gxx variables? (these may
// change from one call to another, so we have to re-check each time)
if      (CCTK_Equals(metric_type, "physical"))
   then cgi.use_Cactus_conformal_metric = false;
else if (CCTK_Equals(metric_type, "static conformal"))
   then cgi.use_Cactus_conformal_metric = (*conformal_state > 0);
else	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_find_horizons(): unknown metric_type=\"%s\"!",
		   metric_type);				/*NOTREACHED*/

// update parameters
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

// get the Cactus time step and decide if we want to output h and/or Theta now
IO_info.time_iteration = cctk_iteration;
IO_info.time           = cctk_time;
IO_info.output_h
   = (IO_info.output_h_every > 0)
     && ((IO_info.time_iteration % IO_info.output_h_every) == 0);
IO_info.output_Theta
   = (IO_info.output_Theta_every > 0)
     && ((IO_info.time_iteration % IO_info.output_Theta_every) == 0);
IO_info.output_mean_curvature
   = (IO_info.output_mean_curvature_every > 0)
     && ((IO_info.time_iteration % IO_info.output_mean_curvature_every) == 0);

// set initial guess for any (genuine) horizons that need it,
// i.e. for any (genuine) horizons where we didn't find the horizon previously
	for (int hn = hs.init_hn() ; hs.is_genuine() ; hn = hs.next_hn())
	{
	assert( state.AH_data_array[hn] != NULL );
	struct AH_data& AH_data = *state.AH_data_array[hn];
        if (verbose_info.print_algorithm_details) {
          printf ("AHF find_horizons[%d] initial_find_flag=%d\n", hn, (int) AH_data.initial_find_flag);
          printf ("AHF find_horizons[%d] really_initial_find_flag=%d\n", hn, (int) AH_data.really_initial_find_flag);
          printf ("AHF find_horizons[%d] search_flag=%d\n", hn, (int) AH_data.search_flag);
          printf ("AHF find_horizons[%d] found_flag=%d\n", hn, (int) AH_data.found_flag);
        }
	if (AH_data.found_flag)
           then {
                AH_data.initial_find_flag = false;
                AH_data.really_initial_find_flag = false;
                }
	   else {
                if (AH_data.really_initial_find_flag
                    || AH_data.initial_guess_info.reset_horizon_after_not_finding)
                   then {
		        patch_system& ps = *AH_data.ps_ptr;
                        if (verbose_info.print_algorithm_details) {
                          printf ("AHF find_horizons[%d] setup_initial_guess\n", hn);
                        }
                        if (track_origin_from_grid_scalar[hn] && state.method == method__find_horizons) {
                           track_origin(cctkGH, ps, &AH_data, hn, verbose_info.print_algorithm_details);
                           set_initial_guess_parameters(AH_data, hn, 
                                                        ps.origin_x(), ps.origin_y(), ps.origin_z());
                        }
        		setup_initial_guess(ps,
        		          	    AH_data.initial_guess_info,
        				    IO_info,
        				    hn, N_horizons, verbose_info);
        		if (active_flag && IO_info.output_initial_guess)
        		   then output_gridfn(ps, gfns::gfn__h,
                                              "h", cgi.GH,
        				      IO_info, IO_info.h_base_file_name,
                                              IO_info.h_min_digits,
        				      hn, verbose_info
        				      .print_algorithm_highlights);
        		AH_data.initial_find_flag = true;
        		}
                }
	}

//
// now the main horizon finding (or other computation)
//
switch	(state.method)
	{
case method__evaluate_expansions:
	do_evaluate_expansions(my_proc, N_horizons,
			       *state.my_hs, state.AH_data_array,
			       cgi, gi, IO_info,
			       error_info, verbose_info,
			       state.timer_handle);
	break;

case method__test_expansion_Jacobians:
	do_test_expansion_Jacobians(my_proc, N_horizons,
				    state.AH_data_array,
				    cgi, gi, Jac_info,
				    (test_all_Jacobian_compute_methods != 0),
				    IO_info, error_info, verbose_info,
				    state.timer_handle);
	break;

case method__find_horizons:
	  {
	if (state.timer_handle >= 0)
	   then CCTK_TimerStartI(state.timer_handle);
	Newton(cctkGH,
	       state.N_procs, state.N_active_procs, my_proc,
	       *state.my_hs, state.AH_data_array,
	       cgi, gi, Jac_info, state.solver_info,
	       IO_info, state.BH_diagnostics_info, broadcast_horizon_shape,
	       error_info, verbose_info,
	       state.isb);
	if (state.timer_handle >= 0)
	   then CCTK_TimerStopI(state.timer_handle);
	break;
	  }

default:
	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   find_horizons(): unknown method=(int)%d!\n"
"                    (this should never happen!)"
		   ,
		   int(state.method));				/*NOTREACHED*/
	}

if (state.timer_handle >= 0)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "timer stats for computation:");
	CCTK_TimerPrintDataI(state.timer_handle, -1);
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function implements  AHFinderDirect::method == "horizon function":
// On processor #0 it evaluates the Theta(h) function for each apparent
// horizon (and does any I/O desired); on other processors it does N_horizons
// dummy evaluations on horizon #0.
//
// Note that if we decide to output h, we output it *after* any Theta(h)
// evaluation or horizon finding has been done, to ensure that all the
// ghost zones are filled in in case we need to print them.
//
// Arguments:
// timer_handle = a valid Cactus timer handle if we want to time the
//		  apparent horizon process, or -ve to skip this
//		  (we only time the computation, not the file I/O)
//
namespace {
void do_evaluate_expansions(int my_proc, int N_horizons,
			    horizon_sequence& hs,
			       struct AH_data* const AH_data_array[],
			    const struct cactus_grid_info& cgi,
			    const struct geometry_info& gi,
			    const struct IO_info& IO_info,
			    const struct error_info& error_info,
			    const struct verbose_info& verbose_info,
			    int timer_handle)
{
const bool active_flag = (my_proc == 0);

if (active_flag)
   then {
	assert( hs.N_horizons() == N_horizons );
	assert( hs.my_N_horizons() == N_horizons );

		for (int hn = hs.init_hn() ;
		     hs.is_genuine() ;
		     hn = hs.next_hn())
		{
		assert( AH_data_array[hn] != NULL );
		struct AH_data& AH_data = *AH_data_array[hn];
		patch_system& ps = *AH_data.ps_ptr;

		if (timer_handle >= 0)
		   then CCTK_TimerStartI(timer_handle);
		jtutil::norm<fp> Theta_norms;
		const bool Theta_ok = expansion(&ps,
                                                AH_data.compute_info,
						cgi, gi,
						error_info, true,// initial eval
						false,	// no Jacobian coeffs
						true,	// yes, print msgs
						&Theta_norms);
		if (timer_handle >= 0)
		   then CCTK_TimerStopI(timer_handle);

		if (IO_info.output_h)
		   then output_gridfn(ps, gfns::gfn__h,
                                      "h", cgi.GH,
				      IO_info, IO_info.h_base_file_name,
                                      IO_info.h_min_digits,
				      hn, verbose_info.print_algorithm_details);

		if (Theta_ok)
		   then {
			CCTK_VInfo(CCTK_THORNSTRING,
			   "   Theta(h) rms-norm %.2e, infinity-norm %.2e",
			   Theta_norms.rms_norm(), Theta_norms.infinity_norm());
			if (IO_info.output_Theta)
			   then output_gridfn(ps, gfns::gfn__Theta,
                                              "Theta", cgi.GH,
					      IO_info, IO_info
						       .Theta_base_file_name,
                                              IO_info.h_min_digits,
					      hn, verbose_info
						  .print_algorithm_details);
			if (IO_info.output_mean_curvature)
			   then output_gridfn(ps, gfns::gfn__mean_curvature,
                                              "mean_curvature", cgi.GH,
					      IO_info, IO_info
						       .mean_curvature_base_file_name,
                                              IO_info.h_min_digits,
					      hn, verbose_info
						  .print_algorithm_details);
			}
		}
	}
   else {
                struct what_to_compute new_compute_info;
		for (int i = 0 ; i < N_horizons ; ++i)
		{
                expansion(NULL, new_compute_info,
			  cgi, gi,
			  error_info, true);	// initial evaluation
		}
	}
}
	  }

//******************************************************************************

//
// This function implements
//  AHFinderDirect::method == "test expansion Jacobians":
// On processor #0 it computes and prints the Jacobian matrix J[Theta(h)]
// function for horizon #1; on other processors it does dummy Jacobian
// computations.
//
// The Jacobian computation may optionally be done in several different
// ways, in which case all the resulting Jacobian matrices are printed,
// as are their differences.  Alternatively, only
// the numerical perturbation computation may be done/printed.
//
// Arguments:
// timer_handle = a valid Cactus timer handle if we want to time the
//		  apparent horizon process, or -ve to skip this
//		  (we only time the computation, not the file I/O)
// test_all_Jacobian_compute_methods
//	= true ==> Test all known methods of computing the Jacobian
//		   matrix, and print all the resulting Jacobian matrices
//		   and their differences.
//	  false ==> Just test/print the numerical perturbation calculation.
//		    (This may be useful if one or more of the other methods
//		    is broken.)
//
namespace {
void do_test_expansion_Jacobians(int my_proc, int N_horizons,
				 struct AH_data* const AH_data_array[],
				 const struct cactus_grid_info& cgi,
				 const struct geometry_info& gi,
				       struct Jacobian_info& Jac_info,
				 bool test_all_Jacobian_compute_methods,
				 const struct IO_info& IO_info,
				 const struct error_info& error_info,
				 const struct verbose_info& verbose_info,
				 int timer_handle)
{
const bool active_flag = (my_proc == 0);
assert(N_horizons >= 1);

const bool print_msg_flag = true;
const int hn = 1;

struct AH_data* const AH_data_ptr = active_flag ? AH_data_array[hn]   : NULL;
patch_system*   const      ps_ptr = active_flag ? AH_data_ptr->ps_ptr : NULL;
struct what_to_compute dummy_compute_info;
struct what_to_compute & compute_info =
	active_flag
	? AH_data_ptr->compute_info
	: dummy_compute_info;

//
// numerical-perturbation Jacobian
//
Jacobian* Jac_NP_ptr = active_flag ? AH_data_ptr->Jac_ptr : NULL;
expansion(ps_ptr, compute_info,
	  cgi, gi,
	  error_info, true);		// initial evaluation
Jac_info.Jacobian_compute_method = Jacobian__numerical_perturbation;
expansion_Jacobian(ps_ptr, Jac_NP_ptr,
                   compute_info,
		   cgi, gi, Jac_info,
		   error_info, true,	// initial evaluation
		   print_msg_flag);

Jacobian* Jac_SD_FDdr_ptr = NULL;
if (test_all_Jacobian_compute_methods)
   then {
	// symbolic differentiation with finite diff d/dr
	Jac_SD_FDdr_ptr = active_flag
			  ? new_Jacobian(Jac_info.Jacobian_store_solve_method,
					 *ps_ptr,
					 verbose_info.print_algorithm_details)
			  : NULL;
	expansion(ps_ptr, compute_info,
		  cgi, gi,
		  error_info, true,	// initial evaluation
		  true);		// compute SD Jacobian coeffs
	Jac_info.Jacobian_compute_method = Jacobian__symbolic_diff_with_FD_dr;
	expansion_Jacobian(ps_ptr, Jac_SD_FDdr_ptr,
                           compute_info,
			   cgi, gi, Jac_info,
			   error_info, true,	// initial evaluation
			   print_msg_flag);
	}

if (active_flag)
   then output_Jacobians(*ps_ptr,
			 Jac_NP_ptr, Jac_SD_FDdr_ptr,
			 IO_info, IO_info.Jacobian_base_file_name,
			 IO_info.h_min_digits,
			 hn, print_msg_flag);
}
	  }

//******************************************************************************

	  }	// namespace AHFinderDirect
