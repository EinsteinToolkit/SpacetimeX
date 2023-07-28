// Newton.cc -- solve Theta(h) = 0 via Newton's method
// $Header$
//
// Newton_solve - driver to solve Theta(h) = 0 for all our horizons
//
/// broadcast_status - broadcast status from active processors to all processors
/// broadcast_horizon_data - broadcast AH data to all processors
///
/// print_status - print Newton-iteration status on each active proc
/// Newton_step - take the Newton step, scaling it down if it's too large
//

#include <iostream>

#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
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
// prototypes for functions local to this file
//

namespace {
bool broadcast_status(const cGH* GH,
		      int N_procs, int N_active_procs,
		      int my_proc, bool my_active_flag,
		      int hn, int iteration,
		      enum expansion_status expansion_status,
		      fp mean_horizon_radius, fp infinity_norm,
		      bool found_this_horizon, bool I_need_more_iterations,
		      struct iteration_status_buffers& isb);
void broadcast_horizon_data(const cGH* GH,
			    bool broadcast_flag, bool broadcast_horizon_shape,
                            struct AH_data& AH_data,
			    struct BH_diagnostics& BH_diagnostics,
			    patch_system& ps,
			    struct horizon_buffers& horizon_buffers);

void print_status(int N_active_procs,
		  const struct iteration_status_buffers& isb);
void Newton_step(patch_system& ps,
		 fp mean_horizon_radius, fp max_allowable_Delta_h_over_h,
		 const struct verbose_info& verbose_info);
	  }

void track_origin(const cGH* const cctkGH, patch_system& ps, 
                  struct AH_data* const AH_data_ptr, 
                  const int hn, const bool print_algorithm_details);


//******************************************************************************
//
// This function reads a coordinate origin from a grid scalar, 
// and sets the patch system's origin to that new value
//
//
void track_origin(const cGH* const cctkGH, patch_system& ps, 
                  struct AH_data* const AH_data_ptr, 
                  const int hn, const bool print_algorithm_details)
{
DECLARE_CCTK_ARGUMENTS;
DECLARE_CCTK_PARAMETERS;

if (AH_data_ptr->depends_on == 0) {
        // move the origin as specified in the grid scalars
        fp const * const ox =
          static_cast<CCTK_REAL const *>
          (CCTK_VarDataPtr (cctkGH, 0, track_origin_source_x[hn]));
        assert (ox);
        fp const * const oy =
          static_cast<CCTK_REAL const *>
          (CCTK_VarDataPtr (cctkGH, 0, track_origin_source_y[hn]));
        assert (oy);
        fp const * const oz =
          static_cast<CCTK_REAL const *>
          (CCTK_VarDataPtr (cctkGH, 0, track_origin_source_z[hn]));
        assert (oz);
        if (print_algorithm_details) {
          std::cout << "AHF tracked position ox " << *ox << " oy " << *oy << " oz " << *oz << std::endl;
        }
        ps.origin_x (*ox);
        ps.origin_y (*oy);
        ps.origin_z (*oz);
}

}


//******************************************************************************

//
// This function tries to finds each horizon assigned to this processor,
// by solving Theta(h) = 0 via Newton's method.  For each horizon found,
// it computes the BH diagnostics, optionally prints them (via CCTK_VInfo()),
// and optionally writes them to a BH diagnostics output file.  It also
// optionally writes the horizon shape itself to a data file.
//
// This function must be called synchronously across all processors,
// and it operates synchronously, returning only when every processor
// is done with all its horizons.  See ./README.parallel for a discussion
// of the parallel/multiprocessor issues and algorithm.
//
void Newton(const cGH* GH,
	    int N_procs, int N_active_procs, int my_proc,
	    horizon_sequence& hs, struct AH_data* const AH_data_array[],
	    const struct cactus_grid_info& cgi,
	    const struct geometry_info& gi,
	    const struct Jacobian_info& Jacobian_info,
	    const struct solver_info& solver_info,
	    const struct IO_info& IO_info,
	    const struct BH_diagnostics_info& BH_diagnostics_info,
	    bool broadcast_horizon_shape,
	    const struct error_info& error_info,
	    const struct verbose_info& verbose_info,
	    struct iteration_status_buffers& isb)
{
cGH const * const cctkGH = GH;
DECLARE_CCTK_ARGUMENTS;
DECLARE_CCTK_PARAMETERS;

const bool my_active_flag = hs.has_genuine_horizons();

// print out which horizons we're finding on this processor
if (hs.has_genuine_horizons())
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "proc %d: searching for horizon%s %s/%d",
		   my_proc,
		   (hs.my_N_horizons() > 1 ? "s" : ""),
		   hs.sequence_string(","), int(N_horizons));
   else CCTK_VInfo(CCTK_THORNSTRING,
		   "proc %d: dummy horizon(s) only",
		   my_proc);

    //
    // each pass through this loop finds a single horizon,
    // or does corresponding dummy-horizon calls
    //
    // note there is no explicit exit test, rather we exit from the middle
    // of the loop (only) when all processors are done with all their genuine
    // horizons
    //
    for (int hn = hs.init_hn() ; ; hn = hs.next_hn())
    {
    if (verbose_info.print_algorithm_details)
       then CCTK_VInfo(CCTK_THORNSTRING,
		       "Newton_solve(): processor %d working on horizon %d",
		       my_proc, hn);

    // only try to find horizons every  find_every  time steps
    const bool horizon_is_genuine =
      hs.is_genuine() && AH_data_array[hn]->search_flag;
    // this is only a pessimistic approximation
    const bool there_is_another_genuine_horizon = hs.is_next_genuine();
    if (verbose_info.print_algorithm_details)
       then {
	    CCTK_VInfo(CCTK_THORNSTRING,
		       "                horizon_is_genuine=%d",
		       int(horizon_is_genuine));
	    CCTK_VInfo(CCTK_THORNSTRING,
		       "                there_is_another_genuine_horizon=%d",
		       int(there_is_another_genuine_horizon));
	    }

    struct AH_data* AH_data_ptr
	                  = horizon_is_genuine ? AH_data_array[hn]   : NULL;

    patch_system* const  ps_ptr = horizon_is_genuine
				  ? AH_data_ptr->ps_ptr : NULL;
    Jacobian*     const Jac_ptr = horizon_is_genuine
				  ? AH_data_ptr->Jac_ptr: NULL;

    if (horizon_is_genuine) {
      // deal with dependent horizons
      if (AH_data_ptr->depends_on > 0) {
        assert (AH_data_ptr->depends_on != hn);
        assert (AH_data_ptr->depends_on < hn);
        // check that the other horizon is on the same processor!
        AH_data *AH_other_ptr = AH_data_array[AH_data_ptr->depends_on];
        assert (AH_other_ptr);
        AH_data_ptr->compute_info.desired_value
          = AH_other_ptr->compute_info.desired_value
          * AH_data_ptr->desired_value_factor
          + AH_data_ptr->desired_value_offset;
        const int gnp = ps_ptr->ghosted_N_grid_points();
        assert (AH_other_ptr->ps_ptr->ghosted_N_grid_points() == gnp);
        memcpy (ps_ptr->ghosted_gridfn_data(gfns::gfn__h),
                AH_other_ptr->ps_ptr->ghosted_gridfn_data(gfns::gfn__h),
                gnp * sizeof(fp));
        ps_ptr->origin_x (AH_other_ptr->ps_ptr->origin_x());
        ps_ptr->origin_y (AH_other_ptr->ps_ptr->origin_y());
        ps_ptr->origin_z (AH_other_ptr->ps_ptr->origin_z());
      }
    }

    if (horizon_is_genuine) {
      if (AH_data_ptr->move_origins
          && AH_data_ptr->depends_on == 0
          && AH_data_ptr->found_flag)
      {
        // move the origin to the centre
        patch_system& ps = *ps_ptr;
        fp cx=0, cy=0, cz=0;    // change of origin
        fp sx=0, sy=0, sz=0;    // change of shape
        if (! predict_origin_movement) {
          size_t N=0;
          for (int pn=0; pn<ps.N_patches(); ++pn) {
            patch& p = ps.ith_patch(pn);
            for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
              for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
                const fp rho = p.rho_of_irho (irho);
                const fp sigma = p.sigma_of_isigma (isigma);
                const fp r = p.ghosted_gridfn (gfns::gfn__h, irho, isigma);
                fp x, y, z;
                p.xyz_of_r_rho_sigma (r, rho, sigma, x, y, z);
                cx += x;
                cy += y;
                cz += z;
                ++N;
              }
            }
          }
          assert (N > 0);
          cx /= N;
          cy /= N;
          cz /= N;
          sx=cx; sy=cy; sz=cz;
        } else {                // if predict_origin_movement
          if (ah_centroid_valid[hn-1]) {
            if (verbose_info.print_algorithm_details) {
              std::cout << "AHF centroid x " << (ah_centroid_x[hn-1]) << " y " << (ah_centroid_y[hn-1]) << " z " << (ah_centroid_z[hn-1]) << " t " << (ah_centroid_t[hn-1]) << std::endl;
            }
            if (ah_centroid_valid_p[hn-1]) {
              // have two previous origins: linear extrapolation
              if (verbose_info.print_algorithm_details) {
                std::cout << "AHF centroid_p x " << (ah_centroid_x_p[hn-1]) << " y " << (ah_centroid_y_p[hn-1]) << " z " << (ah_centroid_z_p[hn-1]) << " t " << (ah_centroid_t_p[hn-1]) << std::endl;
              }
              fp const dt   = ah_centroid_t[hn-1] - ah_centroid_t_p[hn-1];
              fp const dt_n = cctk_time - ah_centroid_t  [hn-1];
              fp const dt_p = cctk_time - ah_centroid_t_p[hn-1];
              fp const timescale =
                fabs (dt) + fabs (dt_p) + fabs (dt_n) +
                0.001 * fabs (cctk_delta_time);
              if (10 * fabs (dt  ) < timescale ||
                  10 * fabs (dt_p) < timescale ||
                  10 * fabs (dt_n) < timescale)
              {
                // the times are too similar
                if (verbose_info.print_algorithm_details) {
                  std::cout << "AHF toosim" << std::endl;
                }
                cx = ah_centroid_x[hn-1];
                cy = ah_centroid_y[hn-1];
                cz = ah_centroid_z[hn-1];
              } else {
                fp const fact_n = + dt_p / dt;
                fp const fact_p = - dt_n / dt;
                if (fabs (fact_n) > 5 || fabs (fact_p) > 5) {
                  // don't trust a large extrapolation
                  if (verbose_info.print_algorithm_details) {
                    std::cout << "AHF notrust" << std::endl;
                  }
                  cx = ah_centroid_x[hn-1];
                  cy = ah_centroid_y[hn-1];
                  cz = ah_centroid_z[hn-1];
                } else {
                  if (verbose_info.print_algorithm_details) {
                    std::cout << "AHF extrap fn " << fact_n << " fp " << fact_p << std::endl;
                  }
                  cx = fact_n * ah_centroid_x[hn-1] + fact_p * ah_centroid_x_p[hn-1];
                  cy = fact_n * ah_centroid_y[hn-1] + fact_p * ah_centroid_y_p[hn-1];
                  cz = fact_n * ah_centroid_z[hn-1] + fact_p * ah_centroid_z_p[hn-1];
                  if (verbose_info.print_algorithm_details) {
                    std::cout << "AHF xp " << (ah_centroid_x_p[hn-1]) << " x " << (ah_centroid_x[hn-1]) << " cx " << cx << std::endl;
                  }
                }
              }
            } else {
              // have one previous origin: constant extrapolation
              if (verbose_info.print_algorithm_details) {
                std::cout << "AHF const" << std::endl;
              }
              cx = ah_centroid_x[hn-1];
              cy = ah_centroid_y[hn-1];
              cz = ah_centroid_z[hn-1];
            }
            if (verbose_info.print_algorithm_details) {
              std::cout << "AHF predicted position cx " << cx << " cy " << cy << " cz " << cz << std::endl;
            }
            cx -= ps.origin_x();
            cy -= ps.origin_y();
            cz -= ps.origin_z();
          }
          sx = ah_centroid_x[hn-1] - ps.origin_x();
          sy = ah_centroid_y[hn-1] - ps.origin_y();
          sz = ah_centroid_z[hn-1] - ps.origin_z();
        } // if predict_origin_movement
        switch (ps.type()) {
        case patch_system::patch_system__full_sphere:
          break;                // do nothing
        case patch_system::patch_system__plus_z_hemisphere:
          cz = 0; sz = 0; break;
        case patch_system::patch_system__plus_xy_quadrant_mirrored:
          cx = cy = 0; sx = sy = 0; break;
        case patch_system::patch_system__plus_xy_quadrant_rotating:
          cx = cy = 0; sx = sy = 0; break;
        case patch_system::patch_system__plus_xz_quadrant_mirrored:
          cx = cz = 0; sx = sz = 0; break;
        case patch_system::patch_system__plus_xz_quadrant_rotating:
          cx = cz = 0; sx = sz = 0; break;
        case patch_system::patch_system__plus_xyz_octant_mirrored:
          cx = cy = cz = 0; sx = sy = sz = 0; break;
        case patch_system::patch_system__plus_xyz_octant_rotating:
          cx = cy = cz = 0; sx = sy = sz = 0; break;
        default: assert(0);
        }
        ps.origin_x (ps.origin_x() + cx);
        ps.origin_y (ps.origin_y() + cy);
        ps.origin_z (ps.origin_z() + cz);
        if (reshape_while_moving) {
          for (int pn=0; pn<ps.N_patches(); ++pn) {
            patch& p = ps.ith_patch(pn);
            for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
              for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
                const fp rho = p.rho_of_irho (irho);
                const fp sigma = p.sigma_of_isigma (isigma);
                fp & r = p.ghosted_gridfn (gfns::gfn__h, irho, isigma);
                fp x, y, z;
                p.xyz_of_r_rho_sigma (r, rho, sigma, x, y, z);
                fp const proj = (sx*x + sy*y + sz*z) / r;
                r -= proj;
              }
            }
          }
          ps.synchronize();
        }
      }
      if (track_origin_from_grid_scalar[hn]) {
         track_origin(cctkGH, *ps_ptr, AH_data_ptr, hn, verbose_info.print_algorithm_details);
      }
      
      // modify the initial guess
      jtutil::norm<fp> norms;
      ps_ptr->ghosted_gridfn_norms (gfns::gfn__h, norms);
      // smoothing:
      ps_ptr->scale_ghosted_gridfn
        (1.0 - AH_data_ptr->smoothing_factor, gfns::gfn__h);
      ps_ptr->add_to_ghosted_gridfn
        (AH_data_ptr->smoothing_factor * norms.mean(), gfns::gfn__h);
      // enlarging:
      ps_ptr->scale_ghosted_gridfn
        (AH_data_ptr->shiftout_factor, gfns::gfn__h);
    }

    bool I_am_pretracking = horizon_is_genuine && AH_data_ptr->use_pretracking;
    bool I_was_pretracking = false;
    bool pretracking_have_upper_bound = false;
    bool pretracking_have_lower_bound = false;
    bool pretracking_was_successful = false;
    fp const old_pretracking_value = I_am_pretracking ? AH_data_ptr->pretracking_value : 0.0;
    fp pretracking_upper_bound;
    fp pretracking_lower_bound;
    bool pretracking_have_horizon_info;
    fp pretracking_mean_expansion;
    for (int pretracking_iter = 0;
         I_am_pretracking
           ? pretracking_iter < AH_data_ptr->pretracking_max_iterations
           : true;
         ++pretracking_iter)
      {
        bool found_this_horizon;
        if (I_am_pretracking || I_was_pretracking) {
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Pretracking: iteration %d", pretracking_iter);
          if (pretracking_have_lower_bound) {
            if (pretracking_have_upper_bound) {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: looking for value %g in [%g,%g]",
                         double(AH_data_ptr->pretracking_value),
                           double(pretracking_lower_bound),
                         double(pretracking_upper_bound));
            } else {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: looking for value %g in [%g,*]",
                         double(AH_data_ptr->pretracking_value),
                         double(pretracking_lower_bound));
            }
          } else {
            if (pretracking_have_upper_bound) {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: looking for value %g in [*,%g]",
                         double(AH_data_ptr->pretracking_value),
                           double(pretracking_upper_bound));
            } else {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: looking for value %g in [*,*]",
                         double(AH_data_ptr->pretracking_value));
            }
          }
          AH_data_ptr->compute_info.desired_value = AH_data_ptr->pretracking_value;
          ps_ptr->ghosted_gridfn_copy (gfns::gfn__h, gfns::gfn__save_h);
          pretracking_have_horizon_info = false;
        }

        if (! I_am_pretracking) {
          if (horizon_is_genuine) {
            if (! AH_data_ptr->initial_guess_info.reset_horizon_after_not_finding) {
              // save the surface for possible backtracking later
              ps_ptr->ghosted_gridfn_copy (gfns::gfn__h, gfns::gfn__save_h);
            }
          }
        }

	struct what_to_compute dummy_compute_info;
	struct what_to_compute & compute_info =
		horizon_is_genuine
		? AH_data_ptr->compute_info
		: dummy_compute_info;

    if (horizon_is_genuine) {
      if (ps_ptr->N_additional_points()) {
        const int gnp = ps_ptr->ghosted_N_grid_points();
        ps_ptr->ghosted_gridfn_data(gfns::gfn__h)[gnp] = 0;
      }
    }

    const int max_iterations
       = horizon_is_genuine
	 ? (AH_data_ptr->initial_find_flag
	    ? solver_info.max_Newton_iterations__initial
	    : solver_info.max_Newton_iterations__subsequent)
	 : INT_MAX;

    int num_Theta_growth_iterations = 0;
    fp previous_Theta_norm = 1.0e30;
    int num_Theta_nonshrink_iterations = 0;
    fp best_Theta_norm = 1.0e30;

	//
	// each pass through this loop does a single Newton iteration
	// on the current horizon (which might be either genuine or dummy)
	//
        bool do_return = false;
	for (int iteration = 1 ; ; ++iteration)
	{
	if (verbose_info.print_algorithm_debug)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "beginning iteration %d (horizon_is_genuine=%d)",
			   iteration, int(horizon_is_genuine));


	//
	// evaluate the expansion Theta(h) and its norms
	// *** this is a synchronous operation across all processors ***
	//
	if (horizon_is_genuine
	    && solver_info.debugging_output_at_each_Newton_iteration)
	   then output_gridfn(*ps_ptr, gfns::gfn__h,
                              "h", GH,
			      IO_info, IO_info.h_base_file_name,
                              IO_info.h_min_digits,
			      hn, verbose_info.print_algorithm_highlights,
			      iteration);

        // calculate the norms also for a surface a bit further out,
        // so that we know how they vary in space.
        // perform this calculation first, so that the real values
        // do not have to be saved.
        const fp epsilon = Jacobian_info.perturbation_amplitude;
	jtutil::norm<fp> shifted_Theta_norms;
        jtutil::norm<fp> shifted_expansion_Theta_norms;
        jtutil::norm<fp> shifted_inner_expansion_Theta_norms;
        jtutil::norm<fp> shifted_product_expansion_Theta_norms;
        jtutil::norm<fp> shifted_mean_curvature_Theta_norms;
        fp shifted_area;
	jtutil::norm<fp> shifted2_Theta_norms;
        jtutil::norm<fp> shifted2_expansion_Theta_norms;
        jtutil::norm<fp> shifted2_inner_expansion_Theta_norms;
        jtutil::norm<fp> shifted2_product_expansion_Theta_norms;
        jtutil::norm<fp> shifted2_mean_curvature_Theta_norms;
        fp shifted2_area;
        if (solver_info.want_expansion_gradients) {
          if (horizon_is_genuine) {
            ps_ptr->add_to_ghosted_gridfn(epsilon, gfns::gfn__h);
            // ps_ptr->scale_ghosted_gridfn(1.0+epsilon, gfns::gfn__h);
          }
          const enum expansion_status raw_shifted_expansion_status
            = expansion(ps_ptr, compute_info,
                        cgi, gi,
                        error_info, (iteration == 1),
                        false,	// no, we don't want Jacobian coeffs
                        false,
                        &shifted_Theta_norms,
                        &shifted_expansion_Theta_norms,
                        &shifted_inner_expansion_Theta_norms,
                        &shifted_product_expansion_Theta_norms,
                        &shifted_mean_curvature_Theta_norms);
          if (horizon_is_genuine) {
            shifted_area = ps_ptr->integrate_gridfn
              (gfns::gfn__one, true, true, true,
               gfns::gfn__h,
               gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                   gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                       gfns::gfn__g_dd_33,
               BH_diagnostics_info.integral_method);
            ps_ptr->add_to_ghosted_gridfn(-epsilon, gfns::gfn__h);
            // ps_ptr->scale_ghosted_gridfn(1.0/(1.0+epsilon), gfns::gfn__h);
            
            ps_ptr->add_to_ghosted_gridfn(-epsilon, gfns::gfn__h);
            // ps_ptr->scale_ghosted_gridfn(1.0/(1.0+epsilon), gfns::gfn__h);
          }
          const enum expansion_status raw_shifted2_expansion_status
            = expansion(ps_ptr, compute_info,
                        cgi, gi,
                        error_info, (iteration == 1),
                        false,	// no, we don't want Jacobian coeffs
                        false,
                        &shifted2_Theta_norms,
                        &shifted2_expansion_Theta_norms,
                        &shifted2_inner_expansion_Theta_norms,
                        &shifted2_product_expansion_Theta_norms,
                        &shifted2_mean_curvature_Theta_norms);
          if (horizon_is_genuine) {
            shifted2_area = ps_ptr->integrate_gridfn
              (gfns::gfn__one, true, true, true,
               gfns::gfn__h,
               gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                   gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                       gfns::gfn__g_dd_33,
               BH_diagnostics_info.integral_method);
            ps_ptr->add_to_ghosted_gridfn(epsilon, gfns::gfn__h);
            // ps_ptr->scale_ghosted_gridfn(1.0+epsilon, gfns::gfn__h);
          }
        } // if want_expansion_gradients

        // now calculate the real values.
	jtutil::norm<fp> Theta_norms;
        jtutil::norm<fp> expansion_Theta_norms;
        jtutil::norm<fp> inner_expansion_Theta_norms;
        jtutil::norm<fp> product_expansion_Theta_norms;
        jtutil::norm<fp> mean_curvature_Theta_norms;
	const enum expansion_status raw_expansion_status
               = expansion(ps_ptr, compute_info,
			   cgi, gi,
			   error_info, (iteration == 1),
			   true,	// yes, we want Jacobian coeffs
			   verbose_info.print_algorithm_details,
                           &Theta_norms,
                           &expansion_Theta_norms,
                           &inner_expansion_Theta_norms,
                           &product_expansion_Theta_norms,
                           &mean_curvature_Theta_norms);
        fp area;
        if (horizon_is_genuine)
	        {
	        area = ps_ptr->integrate_gridfn
			(gfns::gfn__one, true, true, true,
			 gfns::gfn__h,
		         gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
					     gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
							         gfns::gfn__g_dd_33,
	                 BH_diagnostics_info.integral_method);
	        }
	const bool Theta_is_ok = (raw_expansion_status == expansion_success);
        const bool norms_are_ok = horizon_is_genuine && Theta_is_ok;

        if (norms_are_ok) {
          const fp this_norm = Theta_norms.infinity_norm();
          if (this_norm > previous_Theta_norm) {
            ++ num_Theta_growth_iterations;
          } else {
            num_Theta_growth_iterations = 0;
          }
          previous_Theta_norm = this_norm;
          
          if (this_norm >= best_Theta_norm) {
            ++ num_Theta_nonshrink_iterations;
          } else {
            num_Theta_nonshrink_iterations = 0;
            best_Theta_norm = this_norm;
          }
        }

        if (I_am_pretracking && norms_are_ok) {
          pretracking_mean_expansion = expansion_Theta_norms.mean();
          pretracking_have_horizon_info = true;
        }

	if (verbose_info.print_algorithm_debug)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "   Newton_solve(): Theta_is_ok=%d",
			   Theta_is_ok);
		if (norms_are_ok)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "   Theta rms-norm %.1e, infinity-norm %.1e",
				   double(Theta_norms.rms_norm()),
				   double(Theta_norms.infinity_norm()));
		}
	if (horizon_is_genuine && Theta_is_ok
	    && solver_info.debugging_output_at_each_Newton_iteration)
	   then {
		output_gridfn(*ps_ptr, gfns::gfn__Theta,
                              "Theta", GH,
			      IO_info, IO_info.Theta_base_file_name,
                              IO_info.h_min_digits,
			      hn, verbose_info.print_algorithm_highlights,
			      iteration);
		output_gridfn(*ps_ptr, gfns::gfn__mean_curvature,
                              "mean_curvature", GH,
			      IO_info, IO_info.mean_curvature_base_file_name,
                              IO_info.h_min_digits,
			      hn, verbose_info.print_algorithm_highlights,
			      iteration);
		}


	//
	// have we found this horizon?
	// if so, compute and output BH diagnostics
	//
        found_this_horizon
           = (norms_are_ok
              && (I_was_pretracking
                  ? pretracking_was_successful
                  : Theta_norms.infinity_norm() <= solver_info.Theta_norm_for_convergence));
	if (horizon_is_genuine)
	   then AH_data_ptr->found_flag = found_this_horizon;

	if (found_this_horizon && ! I_am_pretracking)
	   then {
		struct BH_diagnostics& BH_diagnostics
			= AH_data_ptr->BH_diagnostics;
                const fp mean_expansion       = expansion_Theta_norms.mean();
                const fp mean_inner_expansion = inner_expansion_Theta_norms.mean();
                const fp mean_product_expansion = product_expansion_Theta_norms.mean();
                const fp mean_mean_curvature  = mean_curvature_Theta_norms.mean();
                // const fp area_gradient                 = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_area                               - area                              ) / epsilon;
                // const fp mean_expansion_gradient       = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_expansion_Theta_norms.mean()       - expansion_Theta_norms.mean()      ) / epsilon;
                // const fp mean_inner_expansion_gradient = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_inner_expansion_Theta_norms.mean() - inner_expansion_Theta_norms.mean()) / epsilon;
                // const fp mean_product_expansion_gradient = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_product_expansion_Theta_norms.mean() - product_expansion_Theta_norms.mean()) / epsilon;
                // const fp mean_mean_curvature_gradient  = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_mean_curvature_Theta_norms.mean()  - mean_curvature_Theta_norms.mean() ) / epsilon;
                const fp area_gradient                 = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_area                               - shifted2_area                              ) / (2*epsilon);
                const fp mean_expansion_gradient       = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_expansion_Theta_norms.mean()       - shifted2_expansion_Theta_norms.mean()      ) / (2*epsilon);
                const fp mean_inner_expansion_gradient = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_inner_expansion_Theta_norms.mean() - shifted2_inner_expansion_Theta_norms.mean()) / (2*epsilon);
                const fp mean_product_expansion_gradient = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_product_expansion_Theta_norms.mean() - shifted2_product_expansion_Theta_norms.mean()) / (2*epsilon);
                const fp mean_mean_curvature_gradient  = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_mean_curvature_Theta_norms.mean()  - shifted2_mean_curvature_Theta_norms.mean() ) / (2*epsilon);
		BH_diagnostics.compute(*ps_ptr,
                                       area,
                                       mean_expansion,
                                       mean_inner_expansion,
                                       mean_product_expansion,
                                       mean_mean_curvature,
                                       area_gradient,
                                       mean_expansion_gradient,
                                       mean_inner_expansion_gradient,
                                       mean_product_expansion_gradient,
                                       mean_mean_curvature_gradient,
                                       BH_diagnostics_info);
		if (IO_info.output_BH_diagnostics)
		   then {
			if (AH_data_ptr->BH_diagnostics_fileptr == NULL)
			   then AH_data_ptr->BH_diagnostics_fileptr
				  = BH_diagnostics.setup_output_file
					(IO_info, N_horizons, hn);
			BH_diagnostics.output(AH_data_ptr
					      ->BH_diagnostics_fileptr,
					      IO_info);
			}
		}


	//
	// see if the expansion is too big
	// (if so, we'll give up on this horizon)
	//
	const bool expansion_is_too_large
		= norms_are_ok
                  && (   Theta_norms.infinity_norm() > solver_info.max_allowable_Theta
                      || (   solver_info.max_allowable_Theta_growth_iterations > 0
                          && num_Theta_growth_iterations > solver_info.max_allowable_Theta_growth_iterations)
                      || (   solver_info.max_allowable_Theta_nonshrink_iterations > 0
                          && num_Theta_nonshrink_iterations > solver_info.max_allowable_Theta_nonshrink_iterations)
                         );


	//
	// compute the mean horizon radius, and if it's too large,
	// then pretend expansion() returned a "surface too large" error status
	//
	jtutil::norm<fp> h_norms;
	if (horizon_is_genuine) {
          ps_ptr->ghosted_gridfn_norms(gfns::gfn__h, h_norms);
        }
	const fp mean_horizon_radius
		= horizon_is_genuine ? h_norms.mean() : 0.0;
	const bool horizon_is_too_large
		= (mean_horizon_radius > solver_info
					 .max_allowable_horizon_radius[hn]);

	const enum expansion_status effective_expansion_status
		= horizon_is_too_large ? expansion_failure__surface_too_large
				       : raw_expansion_status;


	//
	// see if we need more iterations (either on this or another horizon)
	//

	// does *this* horizon need more iterations?
	// i.e. has this horizon's Newton iteration not yet converged?
        const bool this_horizon_needs_more_iterations
	   = horizon_is_genuine && Theta_is_ok
	     && !found_this_horizon
	     && !expansion_is_too_large
	     && !horizon_is_too_large
	     && (iteration < max_iterations);

	// do I (this processor) need to do more iterations
	// on this or a following horizon?
	const bool I_need_more_iterations
	   = this_horizon_needs_more_iterations
	     || there_is_another_genuine_horizon
             || I_am_pretracking;

	if (verbose_info.print_algorithm_details)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "   flags: found_this_horizon=%d",
			   int(found_this_horizon));
		CCTK_VInfo(CCTK_THORNSTRING,
			   "          this_horizon_needs_more_iterations=%d",
			   int(this_horizon_needs_more_iterations));
		CCTK_VInfo(CCTK_THORNSTRING,
			   "          I_need_more_iterations=%d",
			   int(I_need_more_iterations));
		}

	//
	// broadcast iteration status from each active processor
	// to all processors, and inclusive-or the "we need more iterations"
	// flags to see if *any* (active) processor needs more iterations
	//
	const bool any_proc_needs_more_iterations
	  = broadcast_status(GH,
			     N_procs, N_active_procs,
			     my_proc, my_active_flag,
			     hn, iteration, effective_expansion_status,
			     mean_horizon_radius,
			     (norms_are_ok ? Theta_norms.infinity_norm() : 0.0),
			     found_this_horizon, I_need_more_iterations,
			     isb);
	// set found-this-horizon flags
	// for all active processors' non-dummy horizons
		  {
		for (int found_proc = 0 ;
		     found_proc < N_active_procs ;
		     ++found_proc)
		{
		const int found_hn = isb.hn_buffer[found_proc];
		if (found_hn > 0)
		   then AH_data_array[found_hn]->found_flag
				= isb.found_horizon_buffer[found_proc];
		}
		  }
	if (verbose_info.print_algorithm_details)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "          ==> any_proc_needs_more_iterations=%d",
			   int(any_proc_needs_more_iterations));
		}


	//
	// print the iteration status info
	//
	if ((my_proc == 0) && verbose_info.print_algorithm_highlights)
	   then print_status(N_active_procs, isb);


	//
	// for each processor which found a horizon,
	// broadcast its horizon info to all processors
	// and print the BH diagnostics on processor 0
	//
		  {
		for (int found_proc = 0 ;
		     found_proc < N_active_procs ;
		     ++found_proc)
		{
		if (! isb.found_horizon_buffer[found_proc])
		   then continue;

		const int found_hn = isb.hn_buffer[found_proc];
		struct AH_data& found_AH_data = *AH_data_array[found_hn];

		if (verbose_info.print_algorithm_details)
		   then CCTK_VInfo(CCTK_THORNSTRING,
			   "   broadcasting proc %d/horizon %d diagnostics%s",
				   found_proc, found_hn,
				   (broadcast_horizon_shape ? "+shape" : ""));

		broadcast_horizon_data(GH,
				       my_proc == found_proc,
				       broadcast_horizon_shape,
                                       found_AH_data,
				       found_AH_data.BH_diagnostics,
				       *found_AH_data.ps_ptr,
				       found_AH_data.horizon_buffers);

		if ((my_proc == 0) && verbose_info.print_physics_details)
		   then found_AH_data.BH_diagnostics
				     .print(N_horizons, found_hn);
		}
		  }


	//
	// if we found our horizon, maybe output the horizon shape?
	//
	if (found_this_horizon && ! I_am_pretracking)
	   then {
		// printf("will output h/Th/mc: %d/%d/%d\n", IO_info.output_h, IO_info.output_Theta, IO_info.output_mean_curvature); //xxxxxxxxxxxx
		if (IO_info.output_h)
		   then {
			// if this is the first time we've output h for this
			// horizon, maybe output an OpenDX control file?
			if (!AH_data_ptr->h_files_written)
			   then setup_h_files(*ps_ptr, IO_info, hn);
			output_gridfn(*ps_ptr, gfns::gfn__h,
                                      "h", GH,
				      IO_info, IO_info.h_base_file_name,
                                      IO_info.h_min_digits,
				      hn, verbose_info
					  .print_algorithm_highlights);
			}
		if (IO_info.output_Theta)
		   then output_gridfn(*ps_ptr, gfns::gfn__Theta,
                                      "Theta", GH,
				      IO_info, IO_info.Theta_base_file_name,
                                      IO_info.h_min_digits,
				      hn, verbose_info
					  .print_algorithm_highlights);
		if (IO_info.output_mean_curvature)
		   then output_gridfn(*ps_ptr, gfns::gfn__mean_curvature,
                                      "mean_curvature", GH,
				      IO_info, IO_info.mean_curvature_base_file_name,
                                      IO_info.h_min_digits,
				      hn, verbose_info
					  .print_algorithm_highlights);
		}


	//
	// are all processors done with all their genuine horizons?
	// or if this is a single-processor run, are we done with this horizon?
	//
	if (! any_proc_needs_more_iterations)
	   then {
		if (verbose_info.print_algorithm_details)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "==> all processors are finished!");
                do_return = true;
		break;					// *** LOOP EXIT ***
		}
	if ((N_procs == 1) && !this_horizon_needs_more_iterations)
	   then {
		if (verbose_info.print_algorithm_debug)
		   then CCTK_VInfo(CCTK_THORNSTRING,
			   "==> [single processor] Skipping to next horizon");
		break;					// *** LOOP EXIT ***
		}


	//
	// compute the Jacobian matrix
	// *** this is a synchronous operation across all processors ***
	//
	if (verbose_info.print_algorithm_debug)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "   computing Jacobian: genuine/dummy flag %d",
			   int(this_horizon_needs_more_iterations));
	const enum expansion_status
	  Jacobian_status
	     = expansion_Jacobian
		 (this_horizon_needs_more_iterations ? ps_ptr : NULL,
		  this_horizon_needs_more_iterations ? Jac_ptr : NULL,
		  compute_info,
		  cgi, gi, Jacobian_info,
		  error_info, (iteration == 1),
		  verbose_info.print_algorithm_details);
	const bool Jacobian_is_ok = (Jacobian_status == expansion_success);


	//
	// skip to the next horizon unless
	// this is a genuine Jacobian computation, and it went ok
	//
	if (! (this_horizon_needs_more_iterations && Jacobian_is_ok))
	   then {
		if (verbose_info.print_algorithm_debug)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "   skipping to next horizon");
		break;					// *** LOOP EXIT ***
		}


	//
	// compute the Newton step
	//
	if (verbose_info.print_algorithm_details)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "solving linear system for Delta_h (%d unknowns)",
			   Jac_ptr->N_rows());
	const fp rcond
	   = Jac_ptr->solve_linear_system(gfns::gfn__Theta, gfns::gfn__Delta_h,
					  solver_info.linear_solver_pars,
					  verbose_info.print_algorithm_details);
	if	((rcond >= 0.0) && (rcond < 100.0*FP_EPSILON))
	   then {
		CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "Newton_solve: Jacobian matrix is numerically singular!");
		// give up on this horizon
		break;				// *** LOOP CONTROL ***
		}
	if (verbose_info.print_algorithm_details)
	   then {
		if (rcond > 0.0)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "done solving linear system (rcond=%.1e)",
				   double(rcond));
		   else CCTK_VInfo(CCTK_THORNSTRING,
				   "done solving linear system");
		}

	if (solver_info.debugging_output_at_each_Newton_iteration)
	   then output_gridfn(*ps_ptr, gfns::gfn__Delta_h,
                              "Delta_h", GH,
			      IO_info, IO_info.Delta_h_base_file_name,
                              IO_info.h_min_digits,
			      hn, verbose_info.print_algorithm_details,
			      iteration);


	//
	// take the Newton step (scaled if need be)
	//
	Newton_step(*ps_ptr,
		    mean_horizon_radius, solver_info
					 .max_allowable_Delta_h_over_h,
		    verbose_info);

	// end of this Newton iteration
	}

        if (! I_am_pretracking) {
          if (horizon_is_genuine) {
            if (! AH_data_ptr->initial_guess_info.reset_horizon_after_not_finding) {
              if (! found_this_horizon) {
                // the surface failed; backtrack and continue
                AH_data_ptr->ps_ptr->ghosted_gridfn_copy
                  (gfns::gfn__save_h, gfns::gfn__h);
              }
            }
          }
        }

        // exit
        if (do_return) return;				// *** NORMAL RETURN ***

        // break out of the pretracking loop if we are not pretracking
        if (! I_am_pretracking) break;
        
        if (! found_this_horizon) {
          // the surface failed; backtrack and continue
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Pretracking: solving failed; backtracking");
          AH_data_ptr->ps_ptr->ghosted_gridfn_copy (gfns::gfn__save_h, gfns::gfn__h);
          if (pretracking_have_lower_bound) {
            assert (AH_data_ptr->pretracking_value >= pretracking_lower_bound - 1.0e-10 * fabs(pretracking_lower_bound));
          }
          pretracking_have_lower_bound = true;
          pretracking_lower_bound = AH_data_ptr->pretracking_value;
          if (pretracking_have_upper_bound) {
            AH_data_ptr->pretracking_delta = 0.5 * (pretracking_upper_bound - pretracking_lower_bound);
            AH_data_ptr->pretracking_value = pretracking_lower_bound + 0.5 * (pretracking_upper_bound - pretracking_lower_bound);
          } else {
            AH_data_ptr->pretracking_delta = fabs(AH_data_ptr->pretracking_delta);
            AH_data_ptr->pretracking_delta *= 2.0;
            if (AH_data_ptr->pretracking_delta > AH_data_ptr->pretracking_maximum_delta) {
              AH_data_ptr->pretracking_delta = AH_data_ptr->pretracking_maximum_delta;
            }
            AH_data_ptr->pretracking_value += AH_data_ptr->pretracking_delta;
            if (AH_data_ptr->pretracking_value > AH_data_ptr->pretracking_maximum_value) {
              AH_data_ptr->pretracking_value = AH_data_ptr->pretracking_maximum_value;
            }
          }
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Pretracking: new value %g, delta %g",
                     double(AH_data_ptr->pretracking_value),
                     double(AH_data_ptr->pretracking_delta));
          if (pretracking_lower_bound >= AH_data_ptr->pretracking_maximum_value * 0.9999999999) {
            // give up
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Pretracking: upper bound reached, giving up.");
            I_am_pretracking = false;
            I_was_pretracking = true;
            pretracking_was_successful = false;
            // restore old pretracking goal
            AH_data_ptr->pretracking_value = old_pretracking_value;
          } else if (AH_data_ptr->pretracking_delta < AH_data_ptr->pretracking_minimum_delta) {
            // give up
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Pretracking: step size too small, giving up.");
            I_am_pretracking = false;
            I_was_pretracking = true;
            pretracking_was_successful = true;
          }
        } else {
          // the surface was okay
          // get mean expansion
          assert (pretracking_have_horizon_info);
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Pretracking: solving succeeded; expansion is now %g",
                     double(pretracking_mean_expansion));
//           if (fabs(AH_data_ptr->pretracking_value) > 2.0 * AH_data_ptr->pretracking_minimum_delta) {
          if (AH_data_ptr->pretracking_value > AH_data_ptr->pretracking_minimum_value + 1.0e-10 * AH_data_ptr->pretracking_minimum_delta) {
            // it is not yet a horizon
            if (pretracking_have_upper_bound) {
              assert (AH_data_ptr->pretracking_value <= pretracking_upper_bound + 1.0e-10 * fabs(pretracking_upper_bound));
            }
            pretracking_have_upper_bound = true;
            pretracking_upper_bound = AH_data_ptr->pretracking_value;
            if (pretracking_have_lower_bound) {
#if 1
              // TODO
              // move lower bound further down
              pretracking_lower_bound -= pretracking_upper_bound - pretracking_lower_bound;
              if (pretracking_lower_bound < AH_data_ptr->pretracking_minimum_value) pretracking_lower_bound = AH_data_ptr->pretracking_minimum_value;
#endif
              AH_data_ptr->pretracking_delta = 0.5 * (pretracking_lower_bound - pretracking_upper_bound);
              AH_data_ptr->pretracking_value = pretracking_lower_bound + 0.5 * (pretracking_upper_bound - pretracking_lower_bound);
            } else {
              AH_data_ptr->pretracking_delta = - fabs(AH_data_ptr->pretracking_delta);
              AH_data_ptr->pretracking_delta *= 2.0;
              if (- AH_data_ptr->pretracking_delta > AH_data_ptr->pretracking_maximum_delta) {
                AH_data_ptr->pretracking_delta = - AH_data_ptr->pretracking_maximum_delta;
              }
              AH_data_ptr->pretracking_value += AH_data_ptr->pretracking_delta;
              if (AH_data_ptr->pretracking_value < AH_data_ptr->pretracking_minimum_value) {
                AH_data_ptr->pretracking_value = AH_data_ptr->pretracking_minimum_value;
              }
            }
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Pretracking: new value radius %g, delta %g",
                       double(AH_data_ptr->pretracking_value),
                       double(AH_data_ptr->pretracking_delta));
            if (pretracking_upper_bound <= AH_data_ptr->pretracking_minimum_value * 1.00000000001) {
              // give up
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: lower bound reached, giving up.");
              I_am_pretracking = false;
              I_was_pretracking = true;
              pretracking_was_successful = false;
              // restore old pretracking goal
              AH_data_ptr->pretracking_value = old_pretracking_value;
            } else if (- AH_data_ptr->pretracking_delta < AH_data_ptr->pretracking_minimum_delta) {
              // give up
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: step size too small, giving up.");
              I_am_pretracking = false;
              I_was_pretracking = true;
              pretracking_was_successful = true;
            }
          } else {
            // a true horizon was found; we are done
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Pretracking: done.");
            I_am_pretracking = false;
            I_was_pretracking = true;
            pretracking_was_successful = true;
          }
        } // if surface found

    // end of pretracking loop
    }
    I_am_pretracking = false;

    // end of this horizon
    }

// we should never get to here
assert( false );
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function (which must be called on *every* processor) broadcasts
// the per-iteration status information from each active processor to
// all processors.
//
// The present implementation of this function uses the Cactus reduction
// API.  If AHFinderDirect is ported to some other software environment,
// it's probably best to re-implement this function on top of whatever
// interprocessor-broadcast facility that environment provides.
//
// Arguments:
// GH --> The Cactus grid hierarchy.
// N_procs = The total number of processors.
// N_active_procs = The number of active processors.
// my_active_flag = Is this processor an active processor?
// hn,iteration,effective_expansion_status,
// mean_horizon_radius,infinity_norm,
// found_this_horizon,I_need_more_iterations
//	= On an active processors, these are the values to be broadcast.
//	  On a dummy processors, these are ignored.
// isb = (out) Holds both user buffers (set to the broadcast results)
//	       and low-level buffers (used within this function).  If the
//	       buffer pointers are NULL then the buffers are allocated.
//
// Results:
// This function returns the inclusive-or over all active processors,
// of the broadcast I_need_more_iterations flags.
//
namespace {
bool broadcast_status(const cGH* GH,
		      int N_procs, int N_active_procs,
		      int my_proc, bool my_active_flag,
		      int hn, int iteration,
		      enum expansion_status effective_expansion_status,
		      fp mean_horizon_radius, fp infinity_norm,
		      bool found_this_horizon, bool I_need_more_iterations,
		      struct iteration_status_buffers& isb)
{
assert( my_proc >= 0 );
assert( my_proc < N_procs );

//
// We do the broadcast via a Cactus reduction operation (this is a KLUDGE,
// but until Cactus gets a generic interprocessor communications mechanism
// it's the best we can do without mucking with MPI ourself).  To do the
// gather via a reduction, we set up a 2-D buffer whose entries are all
// 0s, except that on each processor the [my_proc] row has the values we
// want to gather.  Then we do a sum-reduction of the buffers across
// processors.  For the actual reduction we treat the buffers as 1-D
// arrays on each processor (this slightly simplifies the code).
//
// To reduce overheads, we do the entire operation with a single (CCTK_REAL)
// Cactus reduce, casting values to CCTK_REAL as necessary (and encoding
// Boolean flags into signs of known-to-be-positive values.
//
// Alas MPI (and thus Cactus) requires the input and output reduce buffers
// to be distinct, so we need two copies of the buffers on each processor.
//

// the buffers are actually 2-D arrays; these are the column numbers
// ... if we wanted to, we could pack hn, iteration, and
//     effective_expansion_status all into a single (64-bit)
//     floating-point value, but it's not worth the trouble...
enum	{
	// CCTK_INT buffer
	buffer_var__hn = 0,	// also encodes found_this_horizon flag
				// in sign: +=true, -=false
	buffer_var__iteration,	// also encodes I_need_more_iterations flag
				// in sign: +=true, -=false
	buffer_var__expansion_status,
	buffer_var__mean_horizon_radius,
	buffer_var__Theta_infinity_norm,
	N_buffer_vars // no comma
	};

//
// allocate buffers if this is the first use
//
if (isb.hn_buffer == NULL)
   then {
	isb.hn_buffer                  = new int [N_active_procs];
	isb.iteration_buffer           = new int [N_active_procs];
	isb.expansion_status_buffer = new enum expansion_status[N_active_procs];
	isb.mean_horizon_radius_buffer = new fp  [N_active_procs];
	isb.Theta_infinity_norm_buffer = new fp  [N_active_procs];
	isb.found_horizon_buffer       = new bool[N_active_procs];

	isb.send_buffer_ptr    = new jtutil::array2d<CCTK_REAL>
						    (0, N_active_procs-1,
						     0, N_buffer_vars-1);
	isb.receive_buffer_ptr = new jtutil::array2d<CCTK_REAL>
						    (0, N_active_procs-1,
						     0, N_buffer_vars-1);
	}
jtutil::array2d<CCTK_REAL>& send_buffer    = *isb.send_buffer_ptr;
jtutil::array2d<CCTK_REAL>& receive_buffer = *isb.receive_buffer_ptr;

//
// pack this processor's values into the reduction buffer
//
jtutil::zero_C_array(send_buffer.N_array(), send_buffer.data_array());
if (my_active_flag)
   then {
	assert( send_buffer.is_valid_i(my_proc) );
	assert( hn >= 0 );		// encoding scheme assumes this
	assert( iteration > 0 );	// encoding scheme assumes this
	send_buffer(my_proc, buffer_var__hn)
		= found_this_horizon ? +hn : -hn;
	send_buffer(my_proc, buffer_var__iteration)
		= I_need_more_iterations ? +iteration : -iteration;
	send_buffer(my_proc, buffer_var__expansion_status)
		= int(effective_expansion_status);
	send_buffer(my_proc, buffer_var__mean_horizon_radius)
		= mean_horizon_radius;
	send_buffer(my_proc, buffer_var__Theta_infinity_norm)
		= infinity_norm;
	}

//
// do the reduction
//

// this name is appropriate for PUGHReduce, caveat user for other drivers :)
const int reduction_handle = CCTK_ReductionArrayHandle("sum");
if (reduction_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "broadcast_status(): can't get sum-reduction handle!");
								/*NOTREACHED*/

const int reduction_status
   = CCTK_ReduceLocArrayToArray1D(GH,
				  -1,	// broadcast results to all processors
				  reduction_handle,
				  static_cast<const void*>
					     (send_buffer   .data_array()),
				  static_cast<      void*>
					     (receive_buffer.data_array()),
				  send_buffer.N_array(),
				  CCTK_VARIABLE_REAL);
if (reduction_status < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "broadcast_status(): error status %d from reduction!",
		   reduction_status);				/*NOTREACHED*/

//
// unpack the reduction buffer back to the high-level result buffers and
// compute the inclusive-or of the broadcast I_need_more_iterations flags
//
bool any_proc_needs_more_iterations = false;
	for (int proc = 0 ; proc < N_active_procs ; ++proc)
	{
	const int hn_temp = static_cast<int>(
			      receive_buffer(proc, buffer_var__hn)
					    );
	isb.hn_buffer[proc] = jtutil::abs(hn_temp);
	isb.found_horizon_buffer[proc] = (hn_temp > 0);

	const int iteration_temp = static_cast<int>(
				     receive_buffer(proc, buffer_var__iteration)
						   );
	isb.iteration_buffer[proc] = jtutil::abs(iteration_temp);
	const bool proc_needs_more_iterations = (iteration_temp > 0);
	any_proc_needs_more_iterations |= proc_needs_more_iterations;

	isb.expansion_status_buffer[proc]
		= static_cast<enum expansion_status>(
		    static_cast<int>(
		      receive_buffer(proc, buffer_var__expansion_status)
				    )
						    );

	isb.mean_horizon_radius_buffer[proc]
		= receive_buffer(proc, buffer_var__mean_horizon_radius);
	isb.Theta_infinity_norm_buffer[proc]
		= receive_buffer(proc, buffer_var__Theta_infinity_norm);
	}

return any_proc_needs_more_iterations;
}
	  }

//******************************************************************************

//
// This function (which must be called on *every* processor) broadcasts
// the BH diagnostics and (ghosted) horizon shape from a specified processor
// to all processors.
//
// The present implementation of this function uses the Cactus reduction
// API.  If AHFinderDirect is ported to some other software environment,
// it's probably best to re-implement this function on top of whatever
// interprocessor-broadcast facility that environment provides.
//
// Arguments:
// GH --> The Cactus grid hierarchy.
// broadcast_flag = true on the broadcasting processor,
//		    false on all other processors
// BH_diagnostics = On the broadcasting processor, this is the BH diagnostics
//		    to broadcast; on all other processors, it's set to the
//		    broadcast BH diagnostics.
// ps = On the broadcasting processor,  gfn__h  is broadcast;
//      on all other processors,  gfn__h  is set to the broadcast values.
// horizon_buffers = Internal buffers for use in the broadcast;
//		     if  N_buffer == 0  then we set N_buffer and allocate
//		     the buffers.
//
namespace {
void broadcast_horizon_data(const cGH* GH,
			    bool broadcast_flag, bool broadcast_horizon_shape,
                            struct AH_data& AH_data,
			    struct BH_diagnostics& BH_diagnostics,
			    patch_system& ps,
			    struct horizon_buffers& horizon_buffers)
{
//
// Implementation notes:
//
// We do the send via a Cactus sum-reduce where the data are 0 on
// all processors except the sending one.
//
// To reduce the interprocessor-communications cost, we actually only
// broadcast the nominal-grid horizon shape; we then do a synchronize
// operation on the patch system to recover the full ghosted-grid shape.
//

if (horizon_buffers.N_buffer == 0)
   then {
	// allocate the buffers
	horizon_buffers.N_buffer
		= BH_diagnostics::N_buffer
		  + (broadcast_horizon_shape ? ps.N_grid_points() : 0)
                  + 4;
	horizon_buffers.send_buffer
		= new CCTK_REAL[horizon_buffers.N_buffer];
	horizon_buffers.receive_buffer
		= new CCTK_REAL[horizon_buffers.N_buffer];
	}

if (broadcast_flag)
   then {
	// pack the data to be broadcast into the send buffer
	BH_diagnostics.copy_to_buffer(horizon_buffers.send_buffer);
        int posn = BH_diagnostics::N_buffer;
	if (broadcast_horizon_shape)
	   then {
			for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
			{
			const patch& p = ps.ith_patch(pn);
				for (int irho = p.min_irho() ;
				     irho <= p.max_irho() ;
				     ++irho)
				{
				for (int isigma = p.min_isigma() ;
				     isigma <= p.max_isigma() ;
				     ++isigma)
				{
				assert( posn < horizon_buffers.N_buffer );
				horizon_buffers.send_buffer[posn++]
				  = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
				}
				}
			}
		}
        horizon_buffers.send_buffer[posn++] = AH_data.initial_find_flag;
        horizon_buffers.send_buffer[posn++] = AH_data.really_initial_find_flag;
        horizon_buffers.send_buffer[posn++] = AH_data.search_flag;
        horizon_buffers.send_buffer[posn++] = AH_data.found_flag;
        assert( posn == horizon_buffers.N_buffer );
	}
   else jtutil::zero_C_array(horizon_buffers.N_buffer,
			     horizon_buffers.send_buffer);

// this name is appropriate for PUGHReduce, caveat user for other drivers :)
const int reduction_handle = CCTK_ReductionArrayHandle("sum");
if (reduction_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "broadcast_horizon_data(): can't get sum-reduction handle!");
								/*NOTREACHED*/

const int status
   = CCTK_ReduceLocArrayToArray1D(GH,
				  -1,	// result broadcast to all processors
				  reduction_handle,
				  static_cast<const void*>
					     (horizon_buffers.send_buffer),
				  static_cast<      void*>
					     (horizon_buffers.receive_buffer),
				  horizon_buffers.N_buffer,
				  CCTK_VARIABLE_REAL);
if (status < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "broadcast_horizon_data(): error status %d from reduction!",
		   status);					/*NOTREACHED*/

if (!broadcast_flag)
   then {
	// unpack the data from the receive buffer
	BH_diagnostics.copy_from_buffer(horizon_buffers.receive_buffer);
        ps.origin_x(BH_diagnostics.origin_x);
        ps.origin_y(BH_diagnostics.origin_y);
        ps.origin_z(BH_diagnostics.origin_z);
	int posn = BH_diagnostics::N_buffer;
	if (broadcast_horizon_shape)
	   then {
			for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
			{
			patch& p = ps.ith_patch(pn);
				for (int irho = p.min_irho() ;
				     irho <= p.max_irho() ;
				     ++irho)
				{
				for (int isigma = p.min_isigma() ;
				     isigma <= p.max_isigma() ;
				     ++isigma)
				{
				assert( posn < horizon_buffers.N_buffer );
				p.ghosted_gridfn(gfns::gfn__h, irho,isigma)
				   = horizon_buffers.receive_buffer[posn++];
				}
				}
			}

		// recover the full ghosted-grid horizon shape
		// (we only broadcast the nominal-grid shape)
		ps.synchronize();
		}
        AH_data.initial_find_flag        = horizon_buffers.receive_buffer[posn++];
        AH_data.really_initial_find_flag = horizon_buffers.receive_buffer[posn++];
        AH_data.search_flag              = horizon_buffers.receive_buffer[posn++];
        AH_data.found_flag               = horizon_buffers.receive_buffer[posn++];
	assert( posn == horizon_buffers.N_buffer );
	}
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function is called on processor #0 to print the status of the
// Newton iteration on each active processor.
//
// Arguments:
// N_active_procs = The number of active processors.
// isb = The high-level buffers here give the information to be printed.
//
namespace {
void print_status(int N_active_procs,
		  const struct iteration_status_buffers& isb)
{
	for (int proc = 0 ; proc < N_active_procs ; ++proc)
	{
	// don't print anything for processors doing dummy evaluations
	if (isb.hn_buffer[proc] == 0)
	   then continue;

	if (isb.expansion_status_buffer[proc] == expansion_success)
	   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   proc %d/horizon %d:it %d r_grid=%#.3g ||Theta||=%.1e",
			   proc, isb.hn_buffer[proc],
			   isb.iteration_buffer[proc],
			   double(isb.mean_horizon_radius_buffer[proc]),
			   double(isb.Theta_infinity_norm_buffer[proc]));
	   else CCTK_VInfo(CCTK_THORNSTRING,
		   "   proc %d/horizon %d: %s",
			   proc, isb.hn_buffer[proc],
			   expansion_status_string(
			     isb.expansion_status_buffer[proc]
						  ));
	}
}
	  }

//******************************************************************************

//
// This function takes the Newton step, scaling it down if it's too large.
//
// Arguments:
// ps = The patch system containing the gridfns h and Delta_h.
// mean_horizon_radius = ||h||_mean
// max_allowable_Delta_h_over_h = The maximum allowable
//				     ||Delta_h||_infinity / ||h||_mean
//				  Any step over this is internally clamped
//				  (scaled down) to this size.
//
namespace {
void Newton_step(patch_system& ps,
		 fp mean_horizon_radius, fp max_allowable_Delta_h_over_h,
		 const struct verbose_info& verbose_info)
{
//
// compute scale factor (1 for small steps, <1 for large steps)
//

const fp max_allowable_Delta_h
	= max_allowable_Delta_h_over_h * mean_horizon_radius;

jtutil::norm<fp> Delta_h_norms;
ps.gridfn_norms(gfns::gfn__Delta_h, Delta_h_norms);
const fp max_Delta_h = Delta_h_norms.infinity_norm();

const fp scale = (max_Delta_h <= max_allowable_Delta_h)
		 ? 1.0
		 : max_allowable_Delta_h / max_Delta_h;

if (verbose_info.print_algorithm_details)
   then {

	if (scale == 1.0)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "h += Delta_h (rms-norm=%.1e, infinity-norm=%.1e)",
			   Delta_h_norms.rms_norm(),
			   Delta_h_norms.infinity_norm());
	   else CCTK_VInfo(CCTK_THORNSTRING,
			   "h += %g * Delta_h (infinity-norm clamped to %.2g)",
			   scale,
			   scale * Delta_h_norms.infinity_norm());
	}


//
// take the Newton step (scaled if necessary)
//
	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		p.ghosted_gridfn(gfns::gfn__h, irho,isigma)
			-= scale * p.gridfn(gfns::gfn__Delta_h, irho,isigma);
		}
		}
	}

if (ps.N_additional_points())
	{
	const int np = ps.N_grid_points();
	const int gnp = ps.ghosted_N_grid_points();
	ps.ghosted_gridfn_data(gfns::gfn__h)[gnp]
		-= scale * ps.gridfn_data(gfns::gfn__Delta_h)[np];
	}
}
	  }

//******************************************************************************

	  }	// namespace AHFinderDirect
