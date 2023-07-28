// expansion.cc -- evaluate expansion Theta of trial horizon surface
// $Header$
//
// <<<prototypes for functions local to this file>>>
// expansion - top-level driver
//
/// setup_xyz_posns - setup global xyz posns of grid points
/// interpolate_geometry - interpolate g_ij and K_ij from Cactus 3-D grid
/// convert_conformal_to_physical - convert conformal gij to physical
///
/// h_is_finite - does function h contain any NaNs/infinities?
/// geometry_is_finite - do geometry vars contain NaN/infty?
///
/// compute_Theta - compute Theta(h) given earlier setup
///

//
// debugging flags
//
#undef GEOMETRY_INTERP_DEBUG	// define this for verbose debugging
				// of geometry interpolator calls/results
#undef GEOMETRY_INTERP_DEBUG2	// define this for even more verbose debugging
				// of geometry interpolator calls/results
#undef GEOMETRY_INTERP_SYNC_SLEEP	// use sleep() to try to ensure we're
					// synchronized across all processors
					// before calling geometry interpolator
#undef GEOMETRY_INTERP_SYNC_BARRIER	// use CCTK_Barrier() to ensure we're
					// synchronized across all processors
					// before calling geometry interpolator
#undef COMPUTE_THETA_DEBUG	// define this for verbose debugging
				// in compute_Theta()

#ifdef GEOMETRY_INTERP_SYNC_SLEEP
#include <unistd.h>		// sleep()
#endif

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "util_Table.h"
#include "cctk.h"
#include "cctk_Arguments.h"

#include "config.h"
#include "stdc.h"
#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"
using jtutil::error_exit;
using jtutil::pow2;
using jtutil::pow3;
using jtutil::pow4;

#include "../patch/coords.hh"
#include "../patch/grid.hh"
#include "../patch/fd_grid.hh"
#include "../patch/patch.hh"
#include "../patch/patch_edge.hh"
#include "../patch/patch_interp.hh"
#include "../patch/ghost_zone.hh"
#include "../patch/patch_system.hh"

#include "../elliptic/Jacobian.hh"

#include "gfns.hh"
#include "gr.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** prototypes for functions local to this file *****
//

namespace {
void setup_xyz_posns(patch_system& ps, bool print_msg_flag);
enum expansion_status
  interpolate_geometry(patch_system* ps_ptr,
		       const struct cactus_grid_info& cgi,
		       const struct geometry_info& gi,
		       const struct error_info& error_info, bool initial_flag,
		       bool print_msg_flag);
void convert_conformal_to_physical(patch_system& ps,
				   bool print_msg_flag);
void Schwarzschild_EF_geometry(patch_system& ps,
			       const struct cactus_grid_info& cgi,
			       const struct geometry_info& gi,
			       bool print_msg_flag);

bool h_is_finite(patch_system& ps,
		 const struct error_info& error_info, bool initial_flag,
		 bool print_msg_flag);
bool geometry_is_finite(patch_system& ps,
			const struct error_info& error_info, bool initial_flag,
			bool print_msg_flag);

bool compute_Theta(patch_system& ps, const struct what_to_compute& comput_info,
		   bool Jacobian_flag,
                   jtutil::norm<fp>* Theta_norms_ptr,
                   jtutil::norm<fp>* expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* inner_expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* product_expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* mean_curvature_Theta_norms_ptr,
		   const struct error_info& error_info, bool initial_flag,
		   bool print_msg_flag);
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// If ps_ptr != NULL, this function computes the LHS function Theta(h),
// and optionally also its Jacobian coefficients (from which the Jacobian
// matrix may be computed later).
//
// If ps_ptr == NULL, this function does a dummy computation, described
// below.
//
// Inputs (angular gridfns, on ghosted grid):
// ... defined on ghosted grid
// ... only values on nominal grid are actually used as input
//	h				# shape of trial surface
//
// Inputs (Cactus 3-D gridfns):
//	gxx,gxy,gxz,gyy,gyz,gzz		# 3-metric $g_{ij}$
//	kxx,kxy,kxz,kyy,kyz,kzz		# extrinsic curvature $K_{ij}$
//
// Outputs (temporaries computed at each grid point)
//	## computed by hand-written code
//	global_[xyz]			# xyz positions of grid points
//	X_ud_*, X_udd_*			# xyz derivative coefficients
//	## computed by Maple-generated code
//	g_uu_{11,12,13,22,23,33}	# $g^{ij}$
//	K				# $K$
//	K_dd_{11,12,13,22,23,33}	# $K^{ij}$
//	partial_d_ln_sqrt_g_d		# $\partial_i \ln \sqrt{g}$
//	partial_d_g_uu_{1,2,3}{11,12,13,22,23,33}	# $\partial_k g^{ij}$
//
// Outputs (angular gridfns, all on the nominal grid):
//	## interpolated from 3-D Cactus grid
//	g_dd_{11,12,13,22,23,33}			# $g_{ij}$
//	K_dd_{11,12,13,22,23,33}			# $K_{ij}$
//	partial_d_g_dd_{1,2,3}{11,12,13,22,23,33}	# $\partial_k g_{ij}$
//	## computed at the nominal grid points
//	Theta						# $\Theta = \Theta(h)$
//
// Arguments:
// ps_ptr --> The patch system, or == NULL to do (only) a dummy computation
//	      in which only the parameter-table setup and a dummy geometry
//	      interpolator call are done, the latter with the number of
//	      interpolation points is set to 0 and all the output array
//	      pointers set to NULL.
// add_to_expansion = A real number which is added to the computed expansion
//		      at each grid point.
// initial_flag = true if this is the first evaluation of  expansion()
//		       for this horizon,
//		  false otherwise;
//		  this is used (only) to select which elements of  error_info
//		  are relevant
// Jacobian_flag = true to compute the Jacobian coefficients,
//		   false to skip this.
// print_msg_flag = true to print status messages,
//		    false to skip this.
// Theta_norms_ptr = (out) If this pointer is non-NULL, the norm object it
//			   points to is updated with all the Theta values
//			   in the grid.  This norm object can then be used
//			   to compute various (gridwise) norms of Theta.
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.
//
enum expansion_status
  expansion(patch_system* ps_ptr,
            const struct what_to_compute& compute_info,
	    const struct cactus_grid_info& cgi,
	    const struct geometry_info& gi,
	    const struct error_info& error_info, bool initial_flag,
	    bool Jacobian_flag /* = false */,
	    bool print_msg_flag /* = false */,
	    jtutil::norm<fp>* Theta_norms_ptr /* = NULL */,
	    jtutil::norm<fp>* expansion_Theta_norms_ptr /* = NULL */,
	    jtutil::norm<fp>* inner_expansion_Theta_norms_ptr /* = NULL */,
	    jtutil::norm<fp>* product_expansion_Theta_norms_ptr /* = NULL */,
	    jtutil::norm<fp>* mean_curvature_Theta_norms_ptr /* = NULL */)
{
const bool active_flag = (ps_ptr != NULL);
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   %sexpansion",
		   active_flag ? "" : "dummy ");

if (active_flag)
   then {
	//
	// normal computation
	//

	// fill in values of all ghosted gridfns in ghost zones
	ps_ptr->synchronize();

	if (gi.check_that_h_is_finite && !h_is_finite(*ps_ptr,
						      error_info, initial_flag,
						      print_msg_flag))
	   then return expansion_failure__surface_nonfinite;
							// *** ERROR RETURN ***

	// set up xyz positions of grid points
	setup_xyz_posns(*ps_ptr, print_msg_flag);
	}

// compute the "geometry" g_ij, K_ij, and partial_k g_ij
if (gi.hardwire_Schwarzschild_EF_geometry)
   then {
	if (active_flag)
	   then Schwarzschild_EF_geometry(*ps_ptr,
					  gi,
					  print_msg_flag);
	}
   else {
	// this is the only function we call unconditionally; it looks at
	// ps_ptr (non-NULL vs NULL) to choose a normal vs dummy computation
	const enum expansion_status status
		= interpolate_geometry(ps_ptr,
				       cgi, gi,
				       error_info, initial_flag,
				       print_msg_flag);
	if (status != expansion_success)
	   then return status;				// *** ERROR RETURN ***
	if (active_flag && cgi.use_Cactus_conformal_metric)
	   then convert_conformal_to_physical(*ps_ptr,
					      print_msg_flag);
	}

if (active_flag)
   then {


	if (gi.check_that_geometry_is_finite
	    && !geometry_is_finite(*ps_ptr,
				   error_info, initial_flag,
				   print_msg_flag))
	   then return expansion_failure__geometry_nonfinite;
							// *** ERROR RETURN ***

        // Ensure that there is a norm object
        const bool want_norms = Theta_norms_ptr;
        jtutil::norm<fp> norms;
        if (compute_info.surface_selection != selection_definition)
           then if (! Theta_norms_ptr) Theta_norms_ptr = &norms;


	// compute remaining gridfns --> $\Theta$
	// and optionally also the Jacobian coefficients
	// by algebraic ops and angular finite differencing
        what_to_compute this_compute_info (compute_info);
        this_compute_info.surface_selection = selection_definition;
	if (!compute_Theta(*ps_ptr, this_compute_info,
			   Jacobian_flag, Theta_norms_ptr,
                           expansion_Theta_norms_ptr,
                           inner_expansion_Theta_norms_ptr,
                           product_expansion_Theta_norms_ptr,
                           mean_curvature_Theta_norms_ptr,
			   error_info, initial_flag,
			   print_msg_flag))
	   then return expansion_failure__gij_not_positive_definite;
							// *** ERROR RETURN ***

        if (compute_info.surface_selection != selection_definition) {
          //
          // Apply correction to find a surface by its areal radius
          //
          // get mean expansion
          fp mean_expansion;
          fp areal_radius;
          switch (compute_info.surface_selection) {
          case selection_mean_coordinate_radius: {
            const int np = ps_ptr->N_grid_points();
            fp sum_expansion = 0;
            fp sum_radius = 0;
            for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
              patch& p = ps_ptr->ith_patch(pn);
              for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
                for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
                  sum_expansion += p.gridfn(gfns::gfn__Theta, irho,isigma);
                  sum_radius += p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
                }
              }
            }
            mean_expansion = sum_expansion / np;
            areal_radius = sum_radius / np;
            break;
          }
          case selection_areal_radius: {
            // get surface area
            const fp area = ps_ptr->integrate_gridfn
              (gfns::gfn__one, true, true, true,
               gfns::gfn__h,
               gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                   gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                       gfns::gfn__g_dd_33,
               patch::integration_method__automatic_choice);
            mean_expansion = Theta_norms_ptr->mean();
            areal_radius = sqrt(area / (4.0*PI));
            break;
          }
          case selection_expansion_mean_coordinate_radius: {
            const int np = ps_ptr->N_grid_points();
            fp sum_expansion = 0;
            fp sum_radius = 0;
            for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
              patch& p = ps_ptr->ith_patch(pn);
              for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
                for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
                  sum_expansion += p.gridfn(gfns::gfn__Theta, irho,isigma);
                  sum_radius += p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
                }
              }
            }
            mean_expansion = sum_expansion / np;
            areal_radius = mean_expansion * sum_radius / np;
            break;
          }
          case selection_expansion_areal_radius: {
            // get surface area
            const fp area = ps_ptr->integrate_gridfn
              (gfns::gfn__one, true, true, true,
               gfns::gfn__h,
               gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                   gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                       gfns::gfn__g_dd_33,
               patch::integration_method__automatic_choice);
            mean_expansion = Theta_norms_ptr->mean();
            areal_radius = mean_expansion * sqrt(area / (4.0*PI));
            break;
          }
          default:
            assert (0);
          } // switch areal_radius_definition
          
          if (! ps_ptr->N_additional_points()) {
            // calculate correction
            const fp correction
              = (- mean_expansion
                 + areal_radius - compute_info.desired_value);
            // apply correction
            ps_ptr->add_to_gridfn(-correction, gfns::gfn__Theta);
          } else {
            const int np = ps_ptr->N_grid_points();
            const int gnp = ps_ptr->ghosted_N_grid_points();
            // apply correction
            const fp correction
              = ps_ptr->ghosted_gridfn_data(gfns::gfn__h)[gnp];
            ps_ptr->add_to_gridfn(-correction, gfns::gfn__Theta);
            ps_ptr->gridfn_data(gfns::gfn__Theta)[np]
              = (mean_expansion - correction
                 - areal_radius + compute_info.desired_value);
          }
          if (want_norms) {
            // recalculate norms
            ps_ptr->gridfn_norms (gfns::gfn__Theta, *Theta_norms_ptr);
          } // if want_norms
        } else {
          //
          // do not apply correction
          //
          if (ps_ptr->N_additional_points()) {
            const int np = ps_ptr->N_grid_points();
            const int gnp = ps_ptr->ghosted_N_grid_points();
            ps_ptr->gridfn_data(gfns::gfn__Theta)[np]
              = ps_ptr->ghosted_gridfn_data(gfns::gfn__h)[gnp];
          }
        } // if use-areal-radius
	}

return expansion_success;				// *** NORMAL RETURN ***
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function sets up the global xyz positions of the grid points
// in the gridfns global_[xyz].  These will be used by interplate_geometry().
//
namespace {
void setup_xyz_posns(patch_system& ps, const bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      xyz positions and derivative coefficients");

#ifdef GEOMETRY_INTERP_DEBUG2
        printf ("AH exp\n");
#endif
	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
                const fp r = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);
		fp local_x, local_y, local_z;
		p.xyz_of_r_rho_sigma(r,rho,sigma, local_x,local_y,local_z);
#ifdef GEOMETRY_INTERP_DEBUG2
                 printf ("   pn=%d irho=%d isigma=%d   x=%g y=%g z=%g\n",
                         pn, irho, isigma, local_x, local_y, local_z);
#endif

		const fp global_x = ps.origin_x() + local_x;
		const fp global_y = ps.origin_y() + local_y;
		const fp global_z = ps.origin_z() + local_z;

		p.gridfn(gfns::gfn__global_x, irho,isigma) = global_x;
		p.gridfn(gfns::gfn__global_y, irho,isigma) = global_y;
		p.gridfn(gfns::gfn__global_z, irho,isigma) = global_z;

		const fp global_xx = global_x * global_x;
		const fp global_xy = global_x * global_y;
		const fp global_xz = global_x * global_z;
		const fp global_yy = global_y * global_y;
		const fp global_yz = global_y * global_z;
		const fp global_zz = global_z * global_z;

		p.gridfn(gfns::gfn__global_xx, irho,isigma) = global_xx;
		p.gridfn(gfns::gfn__global_xy, irho,isigma) = global_xy;
		p.gridfn(gfns::gfn__global_xz, irho,isigma) = global_xz;
		p.gridfn(gfns::gfn__global_yy, irho,isigma) = global_yy;
		p.gridfn(gfns::gfn__global_yz, irho,isigma) = global_yz;
		p.gridfn(gfns::gfn__global_zz, irho,isigma) = global_zz;
		}
		}
	}
}
	  }

//******************************************************************************

//
// If ps_ptr != NULL, this function interpolates the Cactus gridfns
//	gxx...gzz
//	kxx...kzz
//	psi			# optional
// to determine the nominal-grid angular gridfns
//	g_dd_ij
//	partial_d_g_dd_kij
//	K_dd_ij
//	psi			# optional
//	partial_d_psi_k		# optional
// at the nominal-grid trial horizon surface positions given by the
// global_(x,y,z) angular gridfns in the patch system *ps_ptr.  The psi
// interpolation is only done if the cgi.use_Cactus_conformal_metric flag
// is set.  Note that this function ignores the physical-vs-conformal
// semantics of the gridfns; it just interpolates and takes derivatives
// of the stored gridfn values.
//
// If ps_ptr == NULL, this function does (only) the parameter-table
// setup and a a dummy interpolator call, as described in the comments
// to  expansion()  above.
//
// The interpolation is done via  CCTK_InterpGridArrays() .  This has the
// option to return both an overall interpolation status, and a "local"
// status which gives the results of interpolating only the points requested
// on *this* processor; if the local status is available we use it, otherwise
// we fall back to the overall status.
//
// Inputs (angular gridfns, all on the nominal grid):
//	global_[xyz]			# xyz positions of grid points
//
// Inputs (Cactus 3-D gridfns):
//      ahmask                          # excision mask
//	gxx,gxy,gxz,gyy,gyz,gzz		# 3-metric $g_{ij}$
//					# (may be either physical or conformal)
//	kxx,kxy,kxz,kyy,kyz,kzz		# extrinsic curvature $K_{ij}$
//	psi				# optional conformal factor $\psi$
//
// Outputs (angular gridfns, all on the nominal grid):
//      mask                                    # excision mask
//      partial_d_mask_[123]                    # derivatives of the mask
//	g_dd_{11,12,13,22,23,33}		# $\stored{g}_{ij}$
//	K_dd_{11,12,13,22,23,33}		# $K_{ij}$
//	partial_d_g_dd_[123]{11,12,13,22,23,33}	# $\partial_k \stored{g}_{ij}$
//	psi					# (optional) $\psi$
//	partial_d_psi_[123]			# (optional) $\partial_k \psi$
//
// This function may also modify the interpolator parameter table.
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.  Possible
// failure codes are
// * expansion_failure__surface_outside_grid
// * expansion_failure__surface_in_excised_region	// not implemented yet
//
namespace {
enum expansion_status
  interpolate_geometry(patch_system* ps_ptr,
			  const struct cactus_grid_info& cgi,
			  const struct geometry_info& gi,
			  const struct error_info& error_info, bool initial_flag,
			  bool print_msg_flag)
{
const bool active_flag = (ps_ptr != NULL);
const bool psi_flag = cgi.use_Cactus_conformal_metric;

//
// Implementation Notes:
//
// To handle the optional interpolation of psi, we set up all the data
// type and pointer arrays to include psi, but with the psi entries at
// the end.  We then choose the array sizes passed to the interpolator
// to either include or exclude the psi entries as appropriate.
//
// We remember whether or not psi was interpolated on the previous call,
// and only modify the interpolator parameter table if this changes (or
// if this is our first call).
//

if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      interpolating %s from Cactus grid",
		   (psi_flag ? "{g_ij, K_ij, psi}" : "{g_ij, K_ij}"));

int status;

#define CAST_PTR_OR_NULL(type_,ptr_)	\
	(ps_ptr == NULL) ? NULL : static_cast<type_>(ptr_)


//
// ***** interpolation points *****
//
const int N_interp_points = (ps_ptr == NULL) ? 0 : ps_ptr->N_grid_points();
const int interp_coords_type_code = CCTK_VARIABLE_REAL;
const void* const interp_coords[N_GRID_DIMS]
  = {
    CAST_PTR_OR_NULL(const void*, ps_ptr->gridfn_data(gfns::gfn__global_x)),
    CAST_PTR_OR_NULL(const void*, ps_ptr->gridfn_data(gfns::gfn__global_y)),
    CAST_PTR_OR_NULL(const void*, ps_ptr->gridfn_data(gfns::gfn__global_z)),
    };


//
// ***** input arrays *****
//

const CCTK_INT input_array_variable_indices[]
	= {
          cgi.mask_varindex,
	  cgi.g_dd_11_varindex, cgi.g_dd_12_varindex, cgi.g_dd_13_varindex,
				cgi.g_dd_22_varindex, cgi.g_dd_23_varindex,
						      cgi.g_dd_33_varindex,
	  cgi.K_dd_11_varindex, cgi.K_dd_12_varindex, cgi.K_dd_13_varindex,
				cgi.K_dd_22_varindex, cgi.K_dd_23_varindex,
						      cgi.K_dd_33_varindex,
	  cgi.psi_varindex,
	  };
const int N_input_arrays_for_psi = 1;
const int N_input_arrays_dim =   sizeof(input_array_variable_indices)
			       / sizeof(input_array_variable_indices[0]);
const int N_input_arrays_use
	= psi_flag ? N_input_arrays_dim
		   : N_input_arrays_dim - N_input_arrays_for_psi;


//
// ***** output arrays *****
//

const CCTK_INT output_array_type_codes[]
	= {
 // mask             $\partial_x$ mask   $\partial_y$ mask   $\partial_z$ mask
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 // $g_{ij}$         $\partial_x g_{ij}$ $\partial_y g_{ij}$ $\partial_z g_{ij}$
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 // $K_{ij}$
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
		     CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
					 CCTK_VARIABLE_REAL,
 // $\psi$           $\partial_x \psi$   $\partial_y \psi$   $\partial_z \psi$
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
	  };

const CCTK_INT operand_indices[]
	= {
	  0, 0, 0, 0,		// mask, partial_[xyz] mask
	  1, 1, 1, 1,		// g_dd_11, partial_[xyz] g_dd_11
	  2, 2, 2, 2,		// g_dd_12, partial_[xyz] g_dd_12
	  3, 3, 3, 3,		// g_dd_13, partial_[xyz] g_dd_13
	  4, 4, 4, 4,		// g_dd_22, partial_[xyz] g_dd_22
	  5, 5, 5, 5,		// g_dd_23, partial_[xyz] g_dd_23
	  6, 6, 6, 6,		// g_dd_33, partial_[xyz] g_dd_33
	  7, 8, 9, 10, 11, 12,	// K_dd_{11,12,13,22,23,33}
	  13, 13, 13, 13,	// psi, partial_[xyz] psi
	  };
#define DERIV(x)	x
const CCTK_INT operation_codes[]
  = {
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // mask, partial_[xyz] mask
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_11, partial_[xyz] g_dd_11
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_12, partial_[xyz] g_dd_12
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_13, partial_[xyz] g_dd_13
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_22, partial_[xyz] g_dd_22
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_23, partial_[xyz] g_dd_23
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_33, partial_[xyz] g_dd_33
    DERIV(0), DERIV(0), DERIV(0), DERIV(0), DERIV(0), DERIV(0),
					    // K_dd_{11,12,13,22,23,33}
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // psi, partial_[xyz] psi
    };

void* const output_arrays[]
  = {
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__mask)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_mask_1)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_mask_2)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_mask_3)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_11)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_111)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_211)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_311)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_12)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_112)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_212)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_312)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_13)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_113)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_213)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_313)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_22)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_122)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_222)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_322)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_23)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_123)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_223)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_323)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_33)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_133)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_233)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_333)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_11)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_12)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_13)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_22)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_23)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_33)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__psi)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_psi_1)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_psi_2)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_psi_3)),
    };

const int N_output_arrays_for_psi = 4;
const int N_output_arrays_dim
	= sizeof(output_arrays) / sizeof(output_arrays[0]);
const int N_output_arrays_use
	= psi_flag ? N_output_arrays_dim
		   : N_output_arrays_dim - N_output_arrays_for_psi;


//
// ***** parameter table *****
//

// this flag is true if and only if the parameter table already has the
//	suppress_warnings
//	operand_indices
//	operand_codes
// entries for  psi_flag == par_table_psi_flag .
static bool par_table_setup = false;

// if  par_table_setup,
//    this flag is the value of  psi_flag  for the parameter table entries
// otherwise this flag is ignored
static bool par_table_psi_flag = false;

if (par_table_setup && (psi_flag == par_table_psi_flag))
   then {
	// parameter table is already set to just what we need
	// ==> no-op here
	}
   else {
	// store derivative info in interpolator parameter table
	if (print_msg_flag)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "         setting up interpolator derivative info");

	// set a dummy value for the key "suppress_warnings"
	// to tell CCTK_InterpGridArrays() not to print warning messages
	// for points outside the grid
	status = Util_TableSetInt(gi.param_table_handle,
				  0,
				  "suppress_warnings");
	if (status < 0)
	   then error_exit(ERROR_EXIT,
"***** interpolate_geometry():\n"
"        unable to set \"suppress_warnings\" key in interpolator parameter table!\n"
"        Util_TableSetInt() status=%d\n"
			   ,
			   status);				/*NOTREACHED*/

	status = Util_TableSetIntArray(gi.param_table_handle,
				       N_output_arrays_use, operand_indices,
				       "operand_indices");
	if (status < 0)
	   then error_exit(ERROR_EXIT,
"***** interpolate_geometry():\n"
"        unable to set operand_indices in interpolator parameter table!\n"
"        Util_TableSetIntArray() status=%d\n"
			   ,
			   status);				/*NOTREACHED*/

	status = Util_TableSetIntArray(gi.param_table_handle,
				       N_output_arrays_use, operation_codes,
				       "operation_codes");
	if (status < 0)
	   then error_exit(ERROR_EXIT,
"***** interpolate_geometry():\n"
"        unable to set operation_codes in interpolator parameter table!\n"
"        Util_TableSetIntArray() status=%d\n"
			   ,
			   status);				/*NOTREACHED*/

	par_table_setup = true;
	par_table_psi_flag = psi_flag;
	}


//
// ***** the actual interpolation *****
//
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         calling geometry interpolator (%s%d points)",
		   (active_flag ? "" : "dummy: "), N_interp_points);

#ifdef GEOMETRY_INTERP_DEBUG2
	  {
	printf("AHFinderDirect:: proc %d: CCTK_InterpGridArrays() coordinates are:\n",
	       int(CCTK_MyProc(cgi.GH)));
	for (int pt = 0 ; pt < N_interp_points ; ++ pt)
	{
	printf("   pt=%d   [x,y,z]=[%g,%g,%g]\n",
               pt,
               double(((const CCTK_REAL*)interp_coords[0])[pt]),
               double(((const CCTK_REAL*)interp_coords[1])[pt]),
               double(((const CCTK_REAL*)interp_coords[2])[pt]));
	}
	  }
#endif	/* GEOMETRY_INTERP_DEBUG2 */

#ifdef GEOMETRY_INTERP_DEBUG
printf("AHFinderDirect:: proc %d: initializing interpolator outputs to 999.999\n",
       int(CCTK_MyProc(cgi.GH)));
	  {
	for (int pt = 0 ; pt < N_interp_points ; ++pt)
	{
		for (int out = 0 ; out < N_output_arrays_use ; ++out)
		{
		CCTK_REAL* const out_ptr
			= static_cast<CCTK_REAL*>(output_arrays[out]);
		out_ptr[pt] = 999.999;
		}
	}
	  }
#endif

#ifdef GEOMETRY_INTERP_DEBUG
printf("AHFinderDirect:: proc %d: calling CCTK_InterpGridArrays(N_interp_points=%d)\n",
       int(CCTK_MyProc(cgi.GH)), N_interp_points);
fflush(stdout);
#endif

#ifdef GEOMETRY_INTERP_SYNC_SLEEP
sleep(1);
#endif
#ifdef GEOMETRY_INTERP_SYNC_BARRIER
CCTK_Barrier(cgi.GH);
#endif

status = CCTK_InterpGridArrays(cgi.GH, N_GRID_DIMS,
			       gi.operator_handle, gi.param_table_handle,
			       cgi.coord_system_handle,
			       N_interp_points,
				  interp_coords_type_code,
				  interp_coords,
			       N_input_arrays_use,
				  input_array_variable_indices,
			       N_output_arrays_use,
				  output_array_type_codes,
				  output_arrays);

#ifdef GEOMETRY_INTERP_DEBUG
printf("AHFinderDirect:: proc %d: CCTK_InterpGridArrays() returned status=%d\n",
       int(CCTK_MyProc(cgi.GH)), status);
fflush(stdout);
#endif

#ifdef GEOMETRY_INTERP_DEBUG2
	  {
	for (int pt = 0 ; pt < N_interp_points ; pt = 2*pt + (pt == 0))
	{
	printf("AHFinderDirect:: proc %d: CCTK_InterpGridArrays() results for pt=%d at [x,y,z]=[%g,%g,%g]:\n",
	       int(CCTK_MyProc(cgi.GH)), pt,
               double(((const CCTK_REAL*)interp_coords[0])[pt]),
               double(((const CCTK_REAL*)interp_coords[1])[pt]),
               double(((const CCTK_REAL*)interp_coords[2])[pt]));
		for (int out = 0 ; out < N_output_arrays_use ; ++out)
		{
		const CCTK_REAL* const out_ptr
			= static_cast<const CCTK_REAL*>(output_arrays[out]);
		printf("   out=%d   result=%g\n", out, double(out_ptr[pt]));
		}
	}
	  }
#endif	/* GEOMETRY_INTERP_DEBUG2 */


//
// ***** get the local interpolator status (if available)
//
// CCTK_InterpGridArrays() returns a status that reflects the interpolation
// on *all* processors, but what we really want here is just the status
// for the interpolation of the points that *we* (this processor) asked for
// ==> get the local status if it's available
//
CCTK_INT local_interpolator_status;
const int table_status = Util_TableGetInt(gi.param_table_handle,
					  &local_interpolator_status,
					  "local_interpolator_status");
if	(table_status == 1)
   then {
	// we got the local interpolator status successfully
	// ==> use it for all further checks
	#ifdef GEOMETRY_INTERP_DEBUG
	printf("AHFinderDirect:: proc %d: replacing status=%d with local_interpolator_status=%d\n",
	       int(CCTK_MyProc(cgi.GH)), status, int(local_interpolator_status));
	fflush(stdout);
	#endif
	status = local_interpolator_status;
	}
else if (table_status == UTIL_ERROR_TABLE_NO_SUCH_KEY)
   then {
	// evidently the interpolators don't provide the
	// local interpolator status
	// ==> stick with the status have
	// ==> no-op here
	}
else	error_exit(ERROR_EXIT,
"***** interpolate_geometry():\n"
"        error return %d trying to get local interpolator status\n"
"        from parameter table!  (CCTK_InterpGridArrays() status=%d)\n"
		   ,
		   table_status, status);		/*NOTREACHED*/


//
// ***** handle any interpolation errors *****
//
if (status == CCTK_ERROR_INTERP_POINT_OUTSIDE)
   then {
	if (print_msg_flag)
	   then {
		// see if we can get further info
		const int warn_level
		   = initial_flag
		     ? error_info.warn_level__point_outside__initial
		     : error_info.warn_level__point_outside__subsequent;


		CCTK_VWarn(warn_level, __LINE__, __FILE__, CCTK_THORNSTRING,
"interpolate_geometry():\n"
"        one or more points on the trial horizon surface point\n"
"        is/are outside the grid (or too close to the grid boundary)\n"
"        (in a single-processor run, this may also mean that\n"
"         driver::ghost_size is too small for this geometry interpolator)\n");
		}

	return expansion_failure__surface_outside_grid;	// *** ERROR RETURN ***
	}

else if (status == CCTK_ERROR_INTERP_GHOST_SIZE_TOO_SMALL)
   then error_exit(ERROR_EXIT,
"***** interpolate_geometry(): driver::ghost_size is too small\n"
"                              for this geometry interpolator!\n");
								/*NOTREACHED*/

else if (status < 0)
   then error_exit(ERROR_EXIT,
"***** interpolate_geometry(): error return %d from interpolator!\n",
		   status);					/*NOTREACHED*/

//
// ***** check the interpolated excision mask *****
//
{
bool did_use_excised_gridpoint = false;
if (active_flag)
   then {
        for (int pn = 0 ; pn < ps_ptr->N_patches() ; ++pn)
            {
            patch& p = ps_ptr->ith_patch(pn);
        
            for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
            for (int isigma = p.min_isigma() ; isigma <= p.max_isigma() ; ++isigma)
        	{
                const fp m = p.gridfn(gfns::gfn__mask, irho,isigma);
                const fp m1 = p.gridfn(gfns::gfn__partial_d_mask_1, irho,isigma);
                const fp m2 = p.gridfn(gfns::gfn__partial_d_mask_2, irho,isigma);
                const fp m3 = p.gridfn(gfns::gfn__partial_d_mask_3, irho,isigma);
                if (fabs(m) > 1.0e-12
                    || fabs(m1) > 1.0e-12 || fabs(m2) > 1.0e-12 || fabs(m3) > 1.0e-12)
                   then did_use_excised_gridpoint = true;
                }
            }
        }
if (gi.mask_is_noshrink && did_use_excised_gridpoint)
   then {
	if (print_msg_flag)
	   then {
		// see if we can get further info
		const int warn_level
		   = initial_flag
		     ? error_info.warn_level__point_outside__initial
		     : error_info.warn_level__point_outside__subsequent;


		CCTK_VWarn(warn_level, __LINE__, __FILE__, CCTK_THORNSTRING,
"interpolate_geometry():\n"
"        one or more points on the trial horizon surface point\n"
"        is/are in an excised region (or too close to the excision boundary)\n");
		}

	return expansion_failure__surface_in_excised_region;	// *** ERROR RETURN ***
	}
}

return expansion_success;				// *** NORMAL RETURN ***
}
	  }

//******************************************************************************

//
// This function converts the g_dd_ij and partial_d_g_dd_kij gridfns
// from the Cactus conformal semantics to the physical semantics.  As
// documented in the CactusEinstein/ConformalState thorn guide, the
// Cactus conformal semantics are
//	physical_gij = psi^4 * stored_gij
// From this it's trivial to derive the transformation for the derivatives,
//	partial_k physical_gij = 4 psi^3 (partial_k psi) stored_gij
//				 + psi^4 partial_k stored_gij
//
// Inputs (angular gridfns, all on the nominal grid):
//   psi					# $\psi$
//   partial_d_psi_{1,2,3}			# $\partial_k \psi$
//   g_dd_{11,12,13,22,23,33}			# $\stored{g}_{ij}$
//   partial_d_g_dd_{1,2,3}{11,12,13,22,23,33}	# $\partial_k \stored{g}_{ij}$
//
//
// Outputs (angular gridfns, all on the nominal grid):
//   g_dd_{11,12,13,22,23,33}			# $g_{ij}$
//   partial_d_g_dd_{1,2,3}{11,12,13,22,23,33}	# $\partial_k g_{ij}$
//
namespace {
void convert_conformal_to_physical(patch_system& ps, bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      converting Cactus conformal gij --> physical gij");

    for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
    {
    patch& p = ps.ith_patch(pn);

	for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
	{
	for (int isigma = p.min_isigma() ; isigma <= p.max_isigma() ; ++isigma)
	{
	const fp psi = p.gridfn(gfns::gfn__psi, irho,isigma);
	const fp psi3 = jtutil::pow3(psi);
	const fp psi4 = jtutil::pow4(psi);

	const fp partial_d_psi_1
		= p.gridfn(gfns::gfn__partial_d_psi_1, irho,isigma);
	const fp partial_d_psi_2
		= p.gridfn(gfns::gfn__partial_d_psi_2, irho,isigma);
	const fp partial_d_psi_3
		= p.gridfn(gfns::gfn__partial_d_psi_3, irho,isigma);

	const fp stored_g_dd_11 = p.gridfn(gfns::gfn__g_dd_11, irho, isigma);
	const fp stored_g_dd_12 = p.gridfn(gfns::gfn__g_dd_12, irho, isigma);
	const fp stored_g_dd_13 = p.gridfn(gfns::gfn__g_dd_13, irho, isigma);
	const fp stored_g_dd_22 = p.gridfn(gfns::gfn__g_dd_22, irho, isigma);
	const fp stored_g_dd_23 = p.gridfn(gfns::gfn__g_dd_23, irho, isigma);
	const fp stored_g_dd_33 = p.gridfn(gfns::gfn__g_dd_33, irho, isigma);

	p.gridfn(gfns::gfn__g_dd_11, irho, isigma) *= psi4;
	p.gridfn(gfns::gfn__g_dd_12, irho, isigma) *= psi4;
	p.gridfn(gfns::gfn__g_dd_13, irho, isigma) *= psi4;
	p.gridfn(gfns::gfn__g_dd_22, irho, isigma) *= psi4;
	p.gridfn(gfns::gfn__g_dd_23, irho, isigma) *= psi4;
	p.gridfn(gfns::gfn__g_dd_33, irho, isigma) *= psi4;

	p.gridfn(gfns::gfn__partial_d_g_dd_111, irho,isigma)
		= 4.0*psi3*partial_d_psi_1*stored_g_dd_11
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_111, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_112, irho,isigma)
		= 4.0*psi3*partial_d_psi_1*stored_g_dd_12
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_112, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_113, irho,isigma)
		= 4.0*psi3*partial_d_psi_1*stored_g_dd_13
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_113, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_122, irho,isigma)
		= 4.0*psi3*partial_d_psi_1*stored_g_dd_22
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_122, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_123, irho,isigma)
		= 4.0*psi3*partial_d_psi_1*stored_g_dd_23
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_123, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_133, irho,isigma)
		= 4.0*psi3*partial_d_psi_1*stored_g_dd_33
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_133, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_211, irho,isigma)
		= 4.0*psi3*partial_d_psi_2*stored_g_dd_11
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_211, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_212, irho,isigma)
		= 4.0*psi3*partial_d_psi_2*stored_g_dd_12
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_212, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_213, irho,isigma)
		= 4.0*psi3*partial_d_psi_2*stored_g_dd_13
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_213, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_222, irho,isigma)
		= 4.0*psi3*partial_d_psi_2*stored_g_dd_22
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_222, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_223, irho,isigma)
		= 4.0*psi3*partial_d_psi_2*stored_g_dd_23
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_223, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_233, irho,isigma)
		= 4.0*psi3*partial_d_psi_2*stored_g_dd_33
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_233, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_311, irho,isigma)
		= 4.0*psi3*partial_d_psi_3*stored_g_dd_11
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_311, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_312, irho,isigma)
		= 4.0*psi3*partial_d_psi_3*stored_g_dd_12
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_312, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_313, irho,isigma)
		= 4.0*psi3*partial_d_psi_3*stored_g_dd_13
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_313, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_322, irho,isigma)
		= 4.0*psi3*partial_d_psi_3*stored_g_dd_22
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_322, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_323, irho,isigma)
		= 4.0*psi3*partial_d_psi_3*stored_g_dd_23
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_323, irho,isigma);
	p.gridfn(gfns::gfn__partial_d_g_dd_333, irho,isigma)
		= 4.0*psi3*partial_d_psi_3*stored_g_dd_33
		  + psi4*p.gridfn(gfns::gfn__partial_d_g_dd_333, irho,isigma);
	}
	}
    }
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function whether the horizon shape function  h  (nominal grid only)
// contains all finite floating-point numbers, or whether it contain one
// or more NaNs or infinities.
//
// Results:
#ifdef HAVE_FINITE
// This function returns true if all the h values are finite, false
// otherwise (i.e. if it contains any NaNs or infinities).
#else
// This function is a no-op, and just prints a warning message and
// returns true.
#endif
//
namespace {
bool h_is_finite(patch_system& ps,
		 const struct error_info& error_info, bool initial_flag,
		 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING, "      checking that h is finite");

#ifdef HAVE_FINITE
	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		const fp h = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
		if (!std::isfinite(h))
		   then {
			const fp rho = p.rho_of_irho(irho);
			const fp sigma = p.sigma_of_isigma(isigma);
			const fp drho   = jtutil::degrees_of_radians(rho);
			const fp dsigma = jtutil::degrees_of_radians(sigma);
			CCTK_VWarn(error_info.warn_level__nonfinite_geometry,
				   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   h=%g isn't finite!\n"
"   %s patch (rho,sigma)=(%g,%g) (drho,dsigma)=(%g,%g)\n"
				   ,
				   double(h),
				   p.name(), double(rho), double(sigma),
					     double(drho), double(dsigma));
			return false;			// *** found a NaN ***
			}
		}
		}
	}
return true;					// *** all values finite ***
#else
CCTK_VWarn(error_info.warn_level__skipping_finite_check,
	   __LINE__, __FILE__, CCTK_THORNSTRING,
           "      no finite() fn ==>  skipping is-h-finite check");
return true;					// *** no check possible ***
#endif
}
	  }

//******************************************************************************

//
// This function whether the geometry variables
//	g_dd_ij
//	partial_d_g_dd_kij
//	K_dd_ij
// are all finite floating-point numbers, or whether they contain one
// or more NaNs or infinities.
//
// Results:
#ifdef HAVE_FINITE
// This function returns true if all the geometry variables are finite,
// false otherwise (i.e. if they contain any NaNs or infinities).
#else
// This function is a no-op, and just prints a warning message and
// returns true.
#endif
//
namespace {
bool geometry_is_finite(patch_system& ps,
			const struct error_info& error_info, bool initial_flag,
			bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING, "      checking that geometry is finite");

#ifdef HAVE_FINITE
    for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
    {
    patch& p = ps.ith_patch(pn);

	for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
	{
	for (int isigma = p.min_isigma() ;
	     isigma <= p.max_isigma() ;
	     ++isigma)
	{
	const fp g_dd_11 = p.gridfn(gfns::gfn__g_dd_11, irho,isigma);
	const fp g_dd_12 = p.gridfn(gfns::gfn__g_dd_12, irho,isigma);
	const fp g_dd_13 = p.gridfn(gfns::gfn__g_dd_13, irho,isigma);
	const fp g_dd_22 = p.gridfn(gfns::gfn__g_dd_22, irho,isigma);
	const fp g_dd_23 = p.gridfn(gfns::gfn__g_dd_23, irho,isigma);
	const fp g_dd_33 = p.gridfn(gfns::gfn__g_dd_33, irho,isigma);

	const fp K_dd_11 = p.gridfn(gfns::gfn__K_dd_11, irho,isigma);
	const fp K_dd_12 = p.gridfn(gfns::gfn__K_dd_12, irho,isigma);
	const fp K_dd_13 = p.gridfn(gfns::gfn__K_dd_13, irho,isigma);
	const fp K_dd_22 = p.gridfn(gfns::gfn__K_dd_22, irho,isigma);
	const fp K_dd_23 = p.gridfn(gfns::gfn__K_dd_23, irho,isigma);
	const fp K_dd_33 = p.gridfn(gfns::gfn__K_dd_33, irho,isigma);

	const fp partial_d_g_dd_111
		= p.gridfn(gfns::gfn__partial_d_g_dd_111, irho,isigma);
	const fp partial_d_g_dd_112
		= p.gridfn(gfns::gfn__partial_d_g_dd_112, irho,isigma);
	const fp partial_d_g_dd_113
		= p.gridfn(gfns::gfn__partial_d_g_dd_113, irho,isigma);
	const fp partial_d_g_dd_122
		= p.gridfn(gfns::gfn__partial_d_g_dd_122, irho,isigma);
	const fp partial_d_g_dd_123
		= p.gridfn(gfns::gfn__partial_d_g_dd_123, irho,isigma);
	const fp partial_d_g_dd_133
		= p.gridfn(gfns::gfn__partial_d_g_dd_133, irho,isigma);
	const fp partial_d_g_dd_211
		= p.gridfn(gfns::gfn__partial_d_g_dd_211, irho,isigma);
	const fp partial_d_g_dd_212
		= p.gridfn(gfns::gfn__partial_d_g_dd_212, irho,isigma);
	const fp partial_d_g_dd_213
		= p.gridfn(gfns::gfn__partial_d_g_dd_213, irho,isigma);
	const fp partial_d_g_dd_222
		= p.gridfn(gfns::gfn__partial_d_g_dd_222, irho,isigma);
	const fp partial_d_g_dd_223
		= p.gridfn(gfns::gfn__partial_d_g_dd_223, irho,isigma);
	const fp partial_d_g_dd_233
		= p.gridfn(gfns::gfn__partial_d_g_dd_233, irho,isigma);
	const fp partial_d_g_dd_311
		= p.gridfn(gfns::gfn__partial_d_g_dd_311, irho,isigma);
	const fp partial_d_g_dd_312
		= p.gridfn(gfns::gfn__partial_d_g_dd_312, irho,isigma);
	const fp partial_d_g_dd_313
		= p.gridfn(gfns::gfn__partial_d_g_dd_313, irho,isigma);
	const fp partial_d_g_dd_322
		= p.gridfn(gfns::gfn__partial_d_g_dd_322, irho,isigma);
	const fp partial_d_g_dd_323
		= p.gridfn(gfns::gfn__partial_d_g_dd_323, irho,isigma);
	const fp partial_d_g_dd_333
		= p.gridfn(gfns::gfn__partial_d_g_dd_333, irho,isigma);

	if (    !std::isfinite(g_dd_11) || !std::isfinite(g_dd_12)
	     || !std::isfinite(g_dd_13) || !std::isfinite(g_dd_22)
	     || !std::isfinite(g_dd_23) || !std::isfinite(g_dd_33)
	     || !std::isfinite(K_dd_11) || !std::isfinite(K_dd_12)
	     || !std::isfinite(K_dd_13) || !std::isfinite(K_dd_22)
	     || !std::isfinite(K_dd_23) || !std::isfinite(K_dd_33)
	     || !std::isfinite(partial_d_g_dd_111)
	     || !std::isfinite(partial_d_g_dd_112)
	     || !std::isfinite(partial_d_g_dd_113)
	     || !std::isfinite(partial_d_g_dd_122)
	     || !std::isfinite(partial_d_g_dd_123)
	     || !std::isfinite(partial_d_g_dd_133)
	     || !std::isfinite(partial_d_g_dd_211)
	     || !std::isfinite(partial_d_g_dd_212)
	     || !std::isfinite(partial_d_g_dd_213)
	     || !std::isfinite(partial_d_g_dd_222)
	     || !std::isfinite(partial_d_g_dd_223)
	     || !std::isfinite(partial_d_g_dd_233)
	     || !std::isfinite(partial_d_g_dd_311)
	     || !std::isfinite(partial_d_g_dd_312)
	     || !std::isfinite(partial_d_g_dd_313)
	     || !std::isfinite(partial_d_g_dd_322)
	     || !std::isfinite(partial_d_g_dd_323)
	     || !std::isfinite(partial_d_g_dd_333)    )
	   then {
		const fp h = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);
		const fp drho   = jtutil::degrees_of_radians(rho);
		const fp dsigma = jtutil::degrees_of_radians(sigma);
		fp local_x, local_y, local_z;
		p.xyz_of_r_rho_sigma(h,rho,sigma, local_x,local_y,local_z);
		const fp global_x = ps.origin_x() + local_x;
		const fp global_y = ps.origin_y() + local_y;
		const fp global_z = ps.origin_z() + local_z;
		CCTK_VWarn(error_info.warn_level__nonfinite_geometry,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   geometry isn't finite at %s patch\n"
"   h=%g (rho,sigma)=(%g,%g) (drho,dsigma)=(%g,%g)\n"
"   local_(x,y,z)=(%g,%g,%g)\n"
"   global_(x,y,z)=(%g,%g,%g)\n"
"   g_dd_11=%g   _12=%g   _13=%g\n"
"       _22=%g   _23=%g   _33=%g\n"
"   K_dd_11=%g   _12=%g   _13=%g\n"
"       _22=%g   _23=%g   _33=%g\n"
"   partial_d_g_dd_111=%g   _112=%g   _113=%g\n"
"                 _122=%g   _123=%g   _133=%g\n"
"   partial_d_g_dd_211=%g   _212=%g   _213=%g\n"
"                 _222=%g   _223=%g   _233=%g\n"
"   partial_d_g_dd_311=%g   _312=%g   _313=%g\n"
"                 _322=%g   _323=%g   _333=%g\n"
			   ,
			   p.name(),
			   double(h), double(rho), double(sigma),
				      double(drho), double(dsigma),
			   double(local_x), double(local_y), double(local_z),
			   double(global_x), double(global_y), double(global_z),
			   double(g_dd_11), double(g_dd_12), double(g_dd_13),
			   double(g_dd_22), double(g_dd_23), double(g_dd_33),
			   double(K_dd_11), double(K_dd_12), double(K_dd_13),
			   double(K_dd_22), double(K_dd_23), double(K_dd_33),
			   double(partial_d_g_dd_111),
			   double(partial_d_g_dd_112),
			   double(partial_d_g_dd_113),
			   double(partial_d_g_dd_122),
			   double(partial_d_g_dd_123),
			   double(partial_d_g_dd_133),
			   double(partial_d_g_dd_211),
			   double(partial_d_g_dd_212),
			   double(partial_d_g_dd_213),
			   double(partial_d_g_dd_222),
			   double(partial_d_g_dd_223),
			   double(partial_d_g_dd_233),
			   double(partial_d_g_dd_311),
			   double(partial_d_g_dd_312),
			   double(partial_d_g_dd_313),
			   double(partial_d_g_dd_322),
			   double(partial_d_g_dd_323),
			   double(partial_d_g_dd_333));
		return false;			// *** found a NaN ***
		}
	}
	}
    }
return true;						// *** no NaNs found ***
#else
CCTK_VWarn(error_info.warn_level__skipping_finite_check,
	   __LINE__, __FILE__, CCTK_THORNSTRING,
           "      no finite() ==>  skipping is-geometry-finite check");
return true;					// *** no check possible ***
#endif
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function computes the expansion Theta(h), and optionally also
// its Jacobian coefficients, (from which the Jacobian matrix may be
// computed later).  This function uses a mixture of algebraic operations
// and (rho,sigma) finite differencing.  The computation is done entirely
// on the nominal angular grid.
//
// N.b. This function #includes "cg.hh", which defines "dangerous" macros
//      which will stay in effect for the rest of this compilation unit!
//
// Arguments:
// Jacobian_flag = true to compute the Jacobian coefficients,
//		   false to skip this.
//
// Results:
// This function returns true for a successful computation, or false
// if the computation failed because Theta_D <= 0 (this means the interpolated
// g_ij isn't positive definite).
//
namespace {
bool compute_Theta(patch_system& ps, const struct what_to_compute& compute_info,
		   bool Jacobian_flag,
                   jtutil::norm<fp>* Theta_norms_ptr,
                   jtutil::norm<fp>* expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* inner_expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* product_expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* mean_curvature_Theta_norms_ptr,
		   const struct error_info& error_info, bool initial_flag,
		   bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING, "      computing Theta(h)");

#if 0
 fp mean_radius, areal_radius;
 switch (compute_info.surface_modification) {
 case modification_none:
 case modification_radius:
 case modification_radius2:
   // do nothing
   break;
 case modification_mean_radius: {
   // get average coordinate radius
   const int np = ps.N_grid_points();
   fp sum_radius = 0;
   for (int pn = 0; pn < ps.N_patches(); ++pn) {
     patch& p = ps.ith_patch(pn);
     for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
       for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
         sum_radius += p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
       }
     }
   }
   mean_radius = sum_radius / np;
   break;
 }
 case modification_areal_radius: {
   // get surface area
   const fp area = ps.integrate_gridfn
     (gfns::gfn__one, true, true, true,
      gfns::gfn__h,
      gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                          gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                              gfns::gfn__g_dd_33,
      patch::integration_method__automatic_choice);
   areal_radius = sqrt(area / (4.0*PI));
   break;
 }
 default:
   assert (0);
 }
#endif

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		//
		// compute the X_ud and X_udd derivative coefficients
		// ... n.b. this uses the *local* (x,y,z) coordinates
		//
		const fp r = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);
		fp xx, yy, zz;
		p.xyz_of_r_rho_sigma(r,rho,sigma, xx, yy, zz);
		#ifdef COMPUTE_THETA_DEBUG
		if (    ((irho == 0) && (isigma == 0))
		     || ((irho == 5) && (isigma == 5))    )
		   then {
			printf("AHFinderDirect:: computing at %s patch (%d,%d) r=%g ==> local xyz=(%g,%g,%g)\n",
			       p.name(), irho,isigma, double(r),
			       double(xx), double(yy), double(zz));
			printf("AHFinderDirect:: got g_dd_11=%g g_dd_12=%g g_dd_13=%g\n",
			       p.gridfn(gfns::gfn__g_dd_11, irho,isigma),
			       p.gridfn(gfns::gfn__g_dd_12, irho,isigma),
			       p.gridfn(gfns::gfn__g_dd_13, irho,isigma));
			printf("                                g_dd_22=%g g_dd_23=%g\n",
			       p.gridfn(gfns::gfn__g_dd_22, irho,isigma),
			       p.gridfn(gfns::gfn__g_dd_23, irho,isigma));
			printf("                                           g_dd_33=%g\n",
			       p.gridfn(gfns::gfn__g_dd_33, irho,isigma));
			fflush(stdout);
			}
		#endif

		// 1st derivative coefficients X_ud
		const fp X_ud_11 = p.partial_rho_wrt_x(xx, yy, zz);
		const fp X_ud_12 = p.partial_rho_wrt_y(xx, yy, zz);
		const fp X_ud_13 = p.partial_rho_wrt_z(xx, yy, zz);
		const fp X_ud_21 = p.partial_sigma_wrt_x(xx, yy, zz);
		const fp X_ud_22 = p.partial_sigma_wrt_y(xx, yy, zz);
		const fp X_ud_23 = p.partial_sigma_wrt_z(xx, yy, zz);

		// 2nd derivative coefficient gridfns X_udd
		const fp X_udd_111 = p.partial2_rho_wrt_xx(xx, yy, zz);
		const fp X_udd_112 = p.partial2_rho_wrt_xy(xx, yy, zz);
		const fp X_udd_113 = p.partial2_rho_wrt_xz(xx, yy, zz);
		const fp X_udd_122 = p.partial2_rho_wrt_yy(xx, yy, zz);
		const fp X_udd_123 = p.partial2_rho_wrt_yz(xx, yy, zz);
		const fp X_udd_133 = p.partial2_rho_wrt_zz(xx, yy, zz);
		const fp X_udd_211 = p.partial2_sigma_wrt_xx(xx, yy, zz);
		const fp X_udd_212 = p.partial2_sigma_wrt_xy(xx, yy, zz);
		const fp X_udd_213 = p.partial2_sigma_wrt_xz(xx, yy, zz);
		const fp X_udd_222 = p.partial2_sigma_wrt_yy(xx, yy, zz);
		const fp X_udd_223 = p.partial2_sigma_wrt_yz(xx, yy, zz);
		const fp X_udd_233 = p.partial2_sigma_wrt_zz(xx, yy, zz);

		//
		// "call" the Maple-generated code
		// ... each cg/*.c file has a separate set of temp variables,
		//     and so must be inside its own set of { } braces
		//

		// gridfn #defines
		#include "cg.hh"

		  {
		// g_uu
		#include "../gr.cg/inverse_metric.c"
		  }

		  {
		// K, K_uu
		#include "../gr.cg/extrinsic_curvature_trace_raise.c"
		  }

		  {
		// partial_d_g_uu
		#include "../gr.cg/inverse_metric_gradient.c"
		  }

		  {
		// partial_d_ln_sqrt_g
		#include "../gr.cg/metric_det_gradient.c"
		  }

		  {
		// Theta_A, Theta_B, Theta_C, Theta_D
		#include "../gr.cg/expansion.c"
		  }

		if (Theta_D <= 0)
		   then {
const int warn_level
  = initial_flag ? error_info.warn_level__gij_not_positive_definite__initial
		 : error_info.warn_level__gij_not_positive_definite__subsequent;
CCTK_VWarn(warn_level, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   compute_Theta(): Theta_D = $g^{ij} s_i s_j$ = %g <= 0\n"
"                    at %s patch rho=%g sigma=%g!\n"
"                    (i.e. the interpolated g_ij isn't positive definite)",
	   double(Theta_D),
	   p.name(), double(rho), double(sigma));
			return false;			// *** ERROR RETURN ***
			}

                assert (compute_info.surface_selection == selection_definition);

		// compute H via equation (14) of my 1996 horizon finding paper
		const fp sqrt_Theta_D = sqrt(Theta_D);

		const fp Theta_X = + Theta_A/(Theta_D*sqrt_Theta_D)
                                   + Theta_B/sqrt_Theta_D;
                const fp Theta_Y = + Theta_C/Theta_D
                                   - K;

#define mean_curvature	p.gridfn(gfns::gfn__mean_curvature, irho,isigma)
		mean_curvature = Theta_X;
#undef mean_curvature

                switch (compute_info.surface_definition) {
                case definition_expansion:
                  Theta = + Theta_X + Theta_Y;
                  break;
                case definition_inner_expansion:
                  Theta = - Theta_X + Theta_Y;
                  break;
                case definition_mean_curvature:
                  Theta = + Theta_X;
                  break;
                case definition_expansion_product:
                  Theta = (+ Theta_X + Theta_Y) * (- Theta_X + Theta_Y);
                  break;
                default:
                  assert (0);
                }

                switch (compute_info.surface_modification) {
                case modification_none:
                  // do nothing
                  break;
                case modification_radius:
                  // multiply by radius
                  Theta *= r;
                  break;
                case modification_radius2:
                  // multiply by radius^2
                  Theta *= pow2(r);
                  break;
#if 0
                case modification_mean_radius:
                  // multiply by average coordinate radius
                  Theta *= mean_radius;
                  break;
                case modification_areal_radius:
                  // multiply by areal radius
                  Theta *= areal_radius;
                  break;
#endif
                default:
                  assert (0);
                }

                Theta -= compute_info.desired_value;
                
		// update running norms of Theta(h) function
		if (Theta_norms_ptr != NULL)
		   then Theta_norms_ptr->data(Theta);

		if (expansion_Theta_norms_ptr != NULL)
                   then expansion_Theta_norms_ptr->data(+ Theta_X + Theta_Y);

		if (inner_expansion_Theta_norms_ptr != NULL)
                   then inner_expansion_Theta_norms_ptr->data(- Theta_X + Theta_Y);

		if (product_expansion_Theta_norms_ptr != NULL)
                   then product_expansion_Theta_norms_ptr->data((+ Theta_X + Theta_Y) * (- Theta_X + Theta_Y));

		if (mean_curvature_Theta_norms_ptr != NULL)
                   then mean_curvature_Theta_norms_ptr->data(+ Theta_X);

                fp partial_Theta_X_wrt_partial_d_h_1;
                fp partial_Theta_X_wrt_partial_d_h_2;
                fp partial_Theta_X_wrt_partial_dd_h_11;
                fp partial_Theta_X_wrt_partial_dd_h_12;
                fp partial_Theta_X_wrt_partial_dd_h_22;
                fp partial_Theta_Y_wrt_partial_d_h_1;
                fp partial_Theta_Y_wrt_partial_d_h_2;
                fp partial_Theta_Y_wrt_partial_dd_h_11;
                fp partial_Theta_Y_wrt_partial_dd_h_12;
                fp partial_Theta_Y_wrt_partial_dd_h_22;

		if (Jacobian_flag)
		   then {
			// partial_Theta_wrt_partial_d_h,
			// partial_Theta_wrt_partial_dd_h
			#include "../gr.cg/expansion_Jacobian.c"
			}

		if (Jacobian_flag) {
                  switch (compute_info.surface_definition) {
                    
                  case definition_expansion:
                    partial_Theta_wrt_partial_d_h_1
                      = (+ partial_Theta_X_wrt_partial_d_h_1
                         + partial_Theta_Y_wrt_partial_d_h_1);
                    partial_Theta_wrt_partial_d_h_2
                      = (+ partial_Theta_X_wrt_partial_d_h_2
                         + partial_Theta_Y_wrt_partial_d_h_2);
                    partial_Theta_wrt_partial_dd_h_11
                      = (+ partial_Theta_X_wrt_partial_dd_h_11
                         + partial_Theta_Y_wrt_partial_dd_h_11);
                    partial_Theta_wrt_partial_dd_h_12
                      = (+ partial_Theta_X_wrt_partial_dd_h_12
                         + partial_Theta_Y_wrt_partial_dd_h_12);
                    partial_Theta_wrt_partial_dd_h_22
                      = (+ partial_Theta_X_wrt_partial_dd_h_22
                         + partial_Theta_Y_wrt_partial_dd_h_22);
                    break;
                    
                  case definition_inner_expansion:
                    partial_Theta_wrt_partial_d_h_1
                      = (- partial_Theta_X_wrt_partial_d_h_1
                         + partial_Theta_Y_wrt_partial_d_h_1);
                    partial_Theta_wrt_partial_d_h_2
                      = (- partial_Theta_X_wrt_partial_d_h_2
                         + partial_Theta_Y_wrt_partial_d_h_2);
                    partial_Theta_wrt_partial_dd_h_11
                      = (- partial_Theta_X_wrt_partial_dd_h_11
                         + partial_Theta_Y_wrt_partial_dd_h_11);
                    partial_Theta_wrt_partial_dd_h_12
                      = (- partial_Theta_X_wrt_partial_dd_h_12
                         + partial_Theta_Y_wrt_partial_dd_h_12);
                    partial_Theta_wrt_partial_dd_h_22
                      = (- partial_Theta_X_wrt_partial_dd_h_22
                         + partial_Theta_Y_wrt_partial_dd_h_22);
                    break;
                    
                  case definition_mean_curvature:
                    partial_Theta_wrt_partial_d_h_1
                      = + partial_Theta_X_wrt_partial_d_h_1;
                    partial_Theta_wrt_partial_d_h_2
                      = + partial_Theta_X_wrt_partial_d_h_2;
                    partial_Theta_wrt_partial_dd_h_11
                      = + partial_Theta_X_wrt_partial_dd_h_11;
                    partial_Theta_wrt_partial_dd_h_12
                      = + partial_Theta_X_wrt_partial_dd_h_12;
                    partial_Theta_wrt_partial_dd_h_22
                      = + partial_Theta_X_wrt_partial_dd_h_22;
                    break;
                    
                  case definition_expansion_product: {
#define f(x,y,dx,dy)  (- x*x + y*y)
#define df(x,y,dx,dy) (- 2*x*dx + 2*y*dy)
                    partial_Theta_wrt_partial_d_h_1   = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_d_h_1  , partial_Theta_Y_wrt_partial_d_h_1  );
                    partial_Theta_wrt_partial_d_h_2   = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_d_h_2  , partial_Theta_Y_wrt_partial_d_h_2  );
                    partial_Theta_wrt_partial_dd_h_11 = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_dd_h_11, partial_Theta_Y_wrt_partial_dd_h_11);
                    partial_Theta_wrt_partial_dd_h_12 = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_dd_h_12, partial_Theta_Y_wrt_partial_dd_h_12);
                    partial_Theta_wrt_partial_dd_h_22 = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_dd_h_22, partial_Theta_Y_wrt_partial_dd_h_22);
#undef f
#undef df
                    break;
                  }
                    
                  default:
                    assert (0);
                  }
                }

		}
		}
	}

return true;						// *** NORMAL RETURN ***
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
