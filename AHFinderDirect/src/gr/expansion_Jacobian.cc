// expansion_Jacobian.cc -- evaluate Jacobian matrix of LHS function Theta(h)
// $Header$
//
// <<<prototypes for functions local to this file>>>
//
// expansion_Jacobian - top-level driver to compute the Jacobian
///
/// expansion_Jacobian_NP - compute the Jacobian by numerical perturbation
/// expansion_Jacobian_partial_SD - compute partial-deriv terms: symbolic diff
/// add_ghost_zone_Jacobian - add ghost zone dependencies to Jacobian
/// expansion_Jacobian_dr_FD - sum d/dr terms (compute via FD) into Jacobian
///

#include <stdio.h>
#include <assert.h>
#include <math.h>

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

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// ***** prototypes for functions local to this file *****
//

namespace {
enum expansion_status
  expansion_Jacobian_NP
	(patch_system& ps, Jacobian& Jac,
         const struct what_to_compute& comput_info,
	 const struct cactus_grid_info& cgi,
	 const struct geometry_info& gi,
	 const struct Jacobian_info& Jacobian_info,
	 const struct error_info& error_info, bool initial_flag,
	 bool print_msg_flag);

void expansion_Jacobian_partial_SD(patch_system& ps, Jacobian& Jac,
				   const struct cactus_grid_info& cgi,
				   const struct geometry_info& gi,
				   const struct Jacobian_info& Jacobian_info,
				   bool print_msg_flag);
void add_ghost_zone_Jacobian(const patch_system& ps,
			     Jacobian& Jac,
			     fp mol,
			     const patch& xp, const ghost_zone& xmgz,
			     int x_II,
			     int xm_irho, int xm_isigma);
enum expansion_status
  expansion_Jacobian_dr_FD
	(patch_system* ps_ptr, Jacobian* Jac_ptr,
         const struct what_to_compute& compute_info,
	 const struct cactus_grid_info& cgi,
	 const struct geometry_info& gi,
	 const struct Jacobian_info& Jacobian_info,
	 const struct error_info& error_info, bool initial_flag,
	 bool print_msg_flag);
	  }

//******************************************************************************

//
// If ps_ptr != NULL and Jac_ptr != NULL, this function computes the
// Jacobian matrix J[Theta(h)] of the expansion Theta(h).  We assume
// that Theta(h) has already been computed.
//
// If ps_ptr == NULL and Jac_ptr == NULL, this function does a dummy
// computation, in which only any expansion() (and hence geometry
// interpolator) calls are done, these with the number of interpolation
// points set to 0 and all the output array pointers set to NULL.
//
// It's illegal for one but not both of ps_ptr and Jac_ptr to be NULL.
//
// Only some values of  Jacobian_info.Jacobian_compute_method  support
// the dummy computation.
//
// Arguments:
// ps_ptr --> The patch system, or == NULL to do (only) a dummy computation.
// Jac_ptr --> The Jacobian, or == NULL to do (only) a dummy computation.
// add_to_expansion = A real number to add to the expansion.
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.
//
enum expansion_status
  expansion_Jacobian(patch_system* ps_ptr, Jacobian* Jac_ptr,
                     const struct what_to_compute& compute_info,
		     const struct cactus_grid_info& cgi,
		     const struct geometry_info& gi,
		     const struct Jacobian_info& Jacobian_info,
		     const struct error_info& error_info, bool initial_flag,
		     bool print_msg_flag /* = false */)
{
const bool active_flag = (ps_ptr != NULL) && (Jac_ptr != NULL);
enum expansion_status status;

switch	(Jacobian_info.Jacobian_compute_method)
	{
case Jacobian__numerical_perturbation:
	if (active_flag)
	   then {
		status = expansion_Jacobian_NP(*ps_ptr, *Jac_ptr,
					       compute_info,
					       cgi, gi, Jacobian_info,
					       error_info, initial_flag,
					       print_msg_flag);
		if (status != expansion_success)
		   then return status;			// *** ERROR RETURN ***
		break;
		}
	   else error_exit(ERROR_EXIT,
"***** expansion_Jacobian():\n"
"        dummy computation isn't supported for\n"
"        Jacobian_compute_method = \"numerical perturbation\"!\n");
								/*NOTREACHED*/

case Jacobian__symbolic_diff:
	error_exit(ERROR_EXIT,
"***** expansion_Jacobian():\n"
"        Jacobian_compute_method == \"symbolic differentiation\"\n"
"        isn't implemented (yet)!\n");				/*NOTREACHED*/

case Jacobian__symbolic_diff_with_FD_dr:
	if (active_flag)
	   then expansion_Jacobian_partial_SD(*ps_ptr, *Jac_ptr,
					      cgi, gi, Jacobian_info,
					      print_msg_flag);
	// this function looks at ps_ptr and Jac_ptr (non-NULL vs NULL)
	// to choose a normal vs dummy computation
	  {
	status = expansion_Jacobian_dr_FD(ps_ptr, Jac_ptr,
                                          compute_info,
					  cgi, gi, Jacobian_info,
					  error_info, initial_flag,
					  print_msg_flag);
	if (status != expansion_success)
	   then return status;				// *** ERROR RETURN ***
	  }
	break;

default:
	error_exit(PANIC_EXIT,
"***** expansion_Jacobian():\n"
"        unknown Jacobian_info.Jacobian_compute_method=(int)%d!\n"
"        (this should never happen!)\n"
,
		   int(Jacobian_info.Jacobian_compute_method));	/*NOTREACHED*/
	}

 if (active_flag) {

   switch (compute_info.surface_modification) {

   case modification_none:
     // do nothing
     break;

   case modification_radius: {
     // multiply with the coordinate radius
     // H_{(r)i} = H_i h_i
     // J_{(r)ij} = J_{ij} h_i + H_i \delta_{ij}
     const int np = ps_ptr->N_grid_points();
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           const int i = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
           const fp radius = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
           for (int j=0; j<np; ++j) {
             if (Jac_ptr->is_explicitly_stored (i, j)) {
               const fp val = Jac_ptr->element (i, j);
               Jac_ptr->set_element (i, j, val * radius);
             }
           }
           const fp Theta = (p.gridfn(gfns::gfn__Theta, irho, isigma)
                             + compute_info.desired_value) / radius;
           Jac_ptr->sum_into_element (i, i, Theta);
         }
       }
     }
     break;
   }

   case modification_radius2: {
     // multiply with the square of the coordinate radius
     // H_{(r2)i} = H_i h_i^2
     // J_{(r2)ij} = J_{ij} h_i^2 + 2 H_i h_i \delta_{ij}
     const int np = ps_ptr->N_grid_points();
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           const int i = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
           const fp radius = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
           const fp radius2 = radius * radius;
           for (int j=0; j<np; ++j) {
             if (Jac_ptr->is_explicitly_stored (i, j)) {
               const fp val = Jac_ptr->element (i, j);
               Jac_ptr->set_element (i, j, val * radius2);
             }
           }
           const fp Theta = (p.gridfn(gfns::gfn__Theta, irho, isigma)
                             + compute_info.desired_value) / radius2;
           Jac_ptr->sum_into_element (i, i, 2 * Theta * radius);
         }
       }
     }
     break;
   }

#if 0
   case modification_mean_radius: {
     // multiply with the average coordinate radius
     // H_{(\bar r)i} = H_i \bar r
     // J_{(\bar r)ij} = J_{ij} \bar r + H_i / N
     // calculate average coordinate radius
     const int np = ps_ptr->N_grid_points();
     fp sum_radius = 0;
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           sum_radius += p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
         }
       }
     }
     mean_radius = sum_radius / np;
     // correct Jacobian
     const int np = ps_ptr->N_grid_points();
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           const int i = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
           for (int j=0; j<np; ++j) {
             if (Jac_ptr->is_explicitly_stored (i, j)) {
               const fp val = Jac_ptr->element (i, j);
               Jac_ptr->set_element (i, j, val * mean_radius);
             }
           }
#error "unfinished"
           const fp Theta = (p.gridfn(gfns::gfn__Theta, irho, isigma)
                             + compute_info.desired_value) / areal_radius;
           const fp dRdh = 0.5 * areal_radius;
           Jac_ptr->sum_into_element (i, i, Theta * dRdh);
         }
       }
     }
     break;
   }

   case modification_areal_radius: {
     // multiply with the areal radius
     // H_{(R)i} = H_i R
     // J_{(R)ij} = J_{ij} R + H_i dR/dh_j
     // get surface area
     const fp area = ps_ptr->integrate_gridfn
       (gfns::gfn__one, true, true, true,
        gfns::gfn__h,
        gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                            gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                gfns::gfn__g_dd_33,
        patch::integration_method__automatic_choice);
     const fp areal_radius = sqrt(area / (4.0*PI));
     // correct Jacobian
     const int np = ps_ptr->N_grid_points();
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           const int i = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
           for (int j=0; j<np; ++j) {
             if (Jac_ptr->is_explicitly_stored (i, j)) {
               const fp val = Jac_ptr->element (i, j);
               Jac_ptr->set_element (i, j, val * areal_radius);
             }
           }
           const fp Theta = (p.gridfn(gfns::gfn__Theta, irho, isigma)
                             + compute_info.desired_value) / areal_radius;
           const fp dRdh = 0.5 * areal_radius;
           Jac_ptr->sum_into_element (i, i, Theta * dRdh);
         }
       }
     }
     break;
   }
#endif

   default:
     assert (0);
   } // switch surface_modification
   
   if (ps_ptr->N_additional_points()) {
     switch (compute_info.surface_selection) {

     case selection_definition: {
       // we want nothing special
       const int np = ps_ptr->N_grid_points();
       for (int i=0; i<np; ++i) {
         Jac_ptr->set_element (i, np, 0.0);
       }
       for (int j=0; j<np; ++j) {
         Jac_ptr->set_element (np, j, 0.0);
       }
       Jac_ptr->set_element (np, np, 1.0);
       break;
     }

     case selection_mean_coordinate_radius: {
       // Jac_ptr->set_element (II, JJ, x) == dTheta(II)/dh(JJ)
       // \frac{\partial R}{\partial h_j} = 1 / N
       const int np = ps_ptr->N_grid_points();
       for (int i=0; i<np; ++i) {
         Jac_ptr->set_element (i, np, -1.0);
       }
       for (int j=0; j<np; ++j) {
         fp val = 0;
         for (int k=0; k<np; ++k) {
           val += Jac_ptr->element (k, j) / np;
         }
         val -= 1.0 / np;
         Jac_ptr->set_element (np, j, val);
       }
       Jac_ptr->set_element (np, np, -1.0);
       break;
     }

     case selection_areal_radius: {
       // \frac{\partial R_a}{\partial h_j}
       //    = \sqrt{1 / 16 \pi A} \sum_k \sqrt{q_k} dS_k
       // The "trapezoid" method is faster
//        const enum patch::integration_method method
//          = patch::integration_method__automatic_choice;
       const enum patch::integration_method method
         = patch::integration_method__trapezoid;
       const fp area = ps_ptr->integrate_gridfn
         (gfns::gfn__one, true, true, true,
          gfns::gfn__h,
          gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
			      gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
						  gfns::gfn__g_dd_33,
          method);
       const int np = ps_ptr->N_grid_points();
       for (int i=0; i<np; ++i) {
         Jac_ptr->set_element (i, np, -1.0);
       }
       for (int j=0; j<np; ++j) {
         fp val = 0;
         for (int k=0; k<np; ++k) {
           val += Jac_ptr->element (k, j) / np;
         }
         Jac_ptr->set_element (np, j, val);
       }
       for (int jpn = 0; jpn < ps_ptr->N_patches(); ++jpn) {
         patch& jp = ps_ptr->ith_patch(jpn);
         for (int jrho = jp.min_irho(); jrho <= jp.max_irho(); ++jrho) {
           for (int jsigma = jp.min_isigma(); jsigma <= jp.max_isigma(); ++jsigma) {
             const int j = ps_ptr->gpn_of_patch_irho_isigma(jp, jrho,jsigma);
             // const fp radius = jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma);
             const fp epsilon = Jacobian_info.perturbation_amplitude;
             fp val1, val2;
#if 0
             // Re-calculate all points
             // (this is slow, but it works)
             val1 = area;
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) += epsilon;
             val2 = ps_ptr->integrate_gridfn
               (gfns::gfn__one, true, true, true,
                gfns::gfn__h,
                gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                    gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                        gfns::gfn__g_dd_33,
                method);
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) -= epsilon;
#else
             // Re-calculate all points with non-zero Jacobian entries
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) -= epsilon/2;
             val1 = 0;
             for (int ipn = 0; ipn < ps_ptr->N_patches(); ++ipn) {
               patch& ip = ps_ptr->ith_patch(ipn);
               for (int irho = ip.min_irho(); irho <= ip.max_irho(); ++irho) {
                 for (int isigma = ip.min_isigma(); isigma <= ip.max_isigma(); ++isigma) {
                   const int i = ps_ptr->gpn_of_patch_irho_isigma(ip, irho,isigma);
                   if (Jac_ptr->is_explicitly_stored (i, j)) {
                     val1 += ps_ptr->integrate_gridpoint
                       (gfns::gfn__one,
                        gfns::gfn__h,
                        gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                            gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                                gfns::gfn__g_dd_33,
                        method,
                        ipn, irho, isigma);
                   }
                 }
               }
             }
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) += epsilon;
             val2 = 0;
             for (int ipn = 0; ipn < ps_ptr->N_patches(); ++ipn) {
               patch& ip = ps_ptr->ith_patch(ipn);
               for (int irho = ip.min_irho(); irho <= ip.max_irho(); ++irho) {
                 for (int isigma = ip.min_isigma(); isigma <= ip.max_isigma(); ++isigma) {
                   const int i = ps_ptr->gpn_of_patch_irho_isigma(ip, irho,isigma);
                   if (Jac_ptr->is_explicitly_stored (i, j)) {
                     val2 += ps_ptr->integrate_gridpoint
                       (gfns::gfn__one,
                        gfns::gfn__h,
                        gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                            gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                                gfns::gfn__g_dd_33,
                        method,
                        ipn, irho, isigma);
                   }
                 }
               }
             }
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) -= epsilon/2;
#endif
             const fp val = 1 / sqrt(16*PI*area) * ps_ptr->integrate_correction(true, true, true) * (val2 - val1) / epsilon;
             Jac_ptr->sum_into_element (np, j, -val);
           }
         }
       }
       Jac_ptr->set_element (np, np, -1.0);
       break;
     }

     case selection_expansion_mean_coordinate_radius: {
       // Jac_ptr->set_element (II, JJ, x) == dTheta(II)/dh(JJ)
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
       for (int i=0; i<np; ++i) {
         Jac_ptr->set_element (i, np, -1.0);
       }
       for (int j=0; j<np; ++j) {
         fp val = 0;
         for (int k=0; k<np; ++k) {
           val += (Jac_ptr->element (k, j) / np) * (1.0 - sum_radius / np);
         }
         val -= (sum_expansion / np) / np;
         Jac_ptr->set_element (np, j, val);
       }
       Jac_ptr->set_element (np, np, -1.0);
       break;
     }

     case selection_expansion_areal_radius: {
       CCTK_WARN (0, "selection_expansion_areal_radius not implemented");
       break;
     }

     default:
       assert (0);
     } // switch surface_selection
   } else {
     assert (compute_info.surface_selection == selection_definition);
   }
   
 } // if active

return expansion_success;				// *** NORMAL RETURN ***
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function computes the Jacobian matrix of the expansion Theta(h)
// by numerical perturbation of the Theta(h) function.  The algorithm is
// as follows:
//
// we assume that Theta = Theta(h) has already been evaluated
// save_Theta = Theta
//	for each point (y,JJ)
//	{
//	const fp save_h_y = h at y;
//	h at y += perturbation_amplitude;
//	evaluate Theta(h) (silently)
//		for each point (x,II)
//		{
//		Jac(II,JJ) = (Theta(II) - save_Theta(II))
//			     / perturbation_amplitude;
//		}
//	h at y = save_h_y;
//	}
// Theta = save_Theta
//
// Inputs (angular gridfns, on ghosted grid):
//	h			# shape of trial surface
//	Theta			# Theta(h) assumed to already be computed
//
// Outputs:
//	The Jacobian matrix is stored in the Jacobian object Jac.
//	As implied by the above algorithm, it's traversed by columns.
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.
//
namespace {
enum expansion_status
  expansion_Jacobian_NP
        (patch_system& ps, Jacobian& Jac,
         const struct what_to_compute& compute_info,
	 const struct cactus_grid_info& cgi,
	 const struct geometry_info& gi,
	 const struct Jacobian_info& Jacobian_info,
	 const struct error_info& error_info, bool initial_flag,
	 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   horizon Jacobian (numerical perturbation)");
const fp epsilon = Jacobian_info.perturbation_amplitude;

ps.gridfn_copy(gfns::gfn__Theta, gfns::gfn__save_Theta);
ps.gridfn_copy(gfns::gfn__mean_curvature, gfns::gfn__save_mean_curvature);

	for (int ypn = 0 ; ypn < ps.N_patches() ; ++ypn)
	{
	patch& yp = ps.ith_patch(ypn);
	if (print_msg_flag)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "      perturbing in %s patch",
			   yp.name());

	for (int y_irho = yp.min_irho() ; y_irho <= yp.max_irho() ; ++y_irho)
	{
	for (int y_isigma = yp.min_isigma() ;
	     y_isigma <= yp.max_isigma() ;
	     ++y_isigma)
	{
	const int JJ = ps.gpn_of_patch_irho_isigma(yp, y_irho,y_isigma);

	const fp save_h_y = yp.ghosted_gridfn(gfns::gfn__h, y_irho,y_isigma);
	yp.ghosted_gridfn(gfns::gfn__h, y_irho,y_isigma) += epsilon;
	const
	  enum expansion_status status = expansion(&ps,
                                                   compute_info,
						   cgi, gi,
						   error_info, initial_flag);
	if (status != expansion_success)
	   then return status;				// *** ERROR RETURN ***

		for (int xpn = 0 ; xpn < ps.N_patches() ; ++xpn)
		{
		patch& xp = ps.ith_patch(xpn);

		for (int x_irho = xp.min_irho() ;
		     x_irho <= xp.max_irho() ;
		     ++x_irho)
		{
		for (int x_isigma = xp.min_isigma() ;
		     x_isigma <= xp.max_isigma() ;
		     ++x_isigma)
		{
		const int II = ps.gpn_of_patch_irho_isigma(xp, x_irho,x_isigma);
		const fp old_Theta = xp.gridfn(gfns::gfn__save_Theta,
					       x_irho,x_isigma);
		const fp new_Theta = xp.gridfn(gfns::gfn__Theta,
					       x_irho,x_isigma);
		Jac.set_element(II,JJ, (new_Theta - old_Theta) / epsilon);
		}
		}
		}

	if (ps.N_additional_points())
		{
		const int np = ps.N_grid_points();
		Jac.set_element(np,JJ, 0.0);	// insert dummy value
		}

	yp.ghosted_gridfn(gfns::gfn__h, y_irho,y_isigma) = save_h_y;
	}
	}
   	} 

ps.gridfn_copy(gfns::gfn__save_Theta, gfns::gfn__Theta);
ps.gridfn_copy(gfns::gfn__save_mean_curvature, gfns::gfn__mean_curvature);
return expansion_success;				// *** NORMAL RETURN ***
}
	  }

//******************************************************************************

//
// This function computes the partial derivative terms in the Jacobian
// matrix of the expansion Theta(h), by symbolic differentiation from
// the Jacobian coefficient (angular) gridfns.  The Jacobian is traversed
// by rows, using equation (25) of my 1996 apparent horizon finding paper.
//
// Inputs (angular gridfns, on ghosted grid):
//	h			# shape of trial surface
//	Theta			# Theta(h) assumed to already be computed
//	partial_Theta_wrt_partial_d_h	# Jacobian coefficients
//	partial_Theta_wrt_partial_dd_h	# (also assumed to already be computed)
//
// Outputs:
//	The Jacobian matrix is stored in the Jacobian object Jac.
//
namespace {
void expansion_Jacobian_partial_SD(patch_system& ps, Jacobian& Jac,
				   const struct cactus_grid_info& cgi,
				   const struct geometry_info& gi,
				   const struct Jacobian_info& Jacobian_info,
				   bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   horizon Jacobian: partial-deriv terms (symbolic diff)");

Jac.zero_matrix();
ps.compute_synchronize_Jacobian();

    for (int xpn = 0 ; xpn < ps.N_patches() ; ++xpn)
    {
    patch& xp = ps.ith_patch(xpn);

	for (int x_irho = xp.min_irho() ; x_irho <= xp.max_irho() ; ++x_irho)
	{
	for (int x_isigma = xp.min_isigma() ;
	x_isigma <= xp.max_isigma() ;
	++x_isigma)
	{
	//
	// compute the main Jacobian terms for this grid point, i.e.
	//	partial Theta(this point x, Jacobian row II)
	//	---------------------------------------------
	//	partial h(other points y, Jacobian column JJ)
	//

	// Jacobian row index
	const int II = ps.gpn_of_patch_irho_isigma(xp, x_irho, x_isigma);

	// Jacobian coefficients for this point
	const fp Jacobian_coeff_rho
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_1,
		       x_irho, x_isigma);
	const fp Jacobian_coeff_sigma
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_2,
		       x_irho, x_isigma);
	const fp Jacobian_coeff_rho_rho
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_11,
		       x_irho, x_isigma);
	const fp Jacobian_coeff_rho_sigma
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_12,
		       x_irho, x_isigma);
	const fp Jacobian_coeff_sigma_sigma
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_22,
		       x_irho, x_isigma);

	// partial_rho, partial_rho_rho
	      {
	    for (int m_irho = xp.molecule_min_m() ;
		 m_irho <= xp.molecule_max_m() ;
		 ++m_irho)
	    {
	    const int xm_irho = x_irho + m_irho;
	    const fp Jac_rho     = Jacobian_coeff_rho
				   * xp.partial_rho_coeff(m_irho);
	    const fp Jac_rho_rho = Jacobian_coeff_rho_rho
				   * xp.partial_rho_rho_coeff(m_irho);
	    const fp Jac_sum = Jac_rho + Jac_rho_rho;
	    if (xp.is_in_nominal_grid(xm_irho, x_isigma))
	       then {
		    const int xm_JJ
		       = Jac.II_of_patch_irho_isigma(xp,xm_irho,x_isigma);
		    Jac.sum_into_element(II, xm_JJ, Jac_sum);
		    }
	       else add_ghost_zone_Jacobian
			(ps, Jac,
			 Jac_sum,
			 xp, xp.minmax_rho_ghost_zone(m_irho < 0),
			 II, xm_irho, x_isigma);
	    }
	      }

	// partial_sigma, partial_sigma_sigma
	      {
	    for (int m_isigma = xp.molecule_min_m() ;
		 m_isigma <= xp.molecule_max_m() ;
		 ++m_isigma)
	    {
	    const int xm_isigma = x_isigma + m_isigma;
	    const fp Jac_sigma       = Jacobian_coeff_sigma
				       * xp.partial_sigma_coeff(m_isigma);
	    const fp Jac_sigma_sigma = Jacobian_coeff_sigma_sigma
				       * xp.partial_sigma_sigma_coeff(m_isigma);
	    const fp Jac_sum = Jac_sigma + Jac_sigma_sigma;
	    if (xp.is_in_nominal_grid(x_irho, xm_isigma))
	       then {
		    const int xm_JJ
		       = Jac.II_of_patch_irho_isigma(xp, x_irho, xm_isigma);
		    Jac.sum_into_element(II, xm_JJ, Jac_sum);
		    }
	       else add_ghost_zone_Jacobian
			(ps, Jac,
			 Jac_sum,
			 xp, xp.minmax_sigma_ghost_zone(m_isigma < 0),
			 II, x_irho, xm_isigma);
	    }
	      }

	// partial_rho_sigma
	      {
	    for (int m_irho = xp.molecule_min_m() ;
		 m_irho <= xp.molecule_max_m() ;
		 ++m_irho)
	    {
	    for (int m_isigma = xp.molecule_min_m() ;
		 m_isigma <= xp.molecule_max_m() ;
		 ++m_isigma)
	    {
	    const int xm_irho   = x_irho   + m_irho;
	    const int xm_isigma = x_isigma + m_isigma;
	    const fp Jac_rho_sigma
	       = Jacobian_coeff_rho_sigma
		 * xp.partial_rho_sigma_coeff(m_irho, m_isigma);
	    if (xp.is_in_nominal_grid(xm_irho, xm_isigma))
	       then {
		    const int xm_JJ
		       = Jac.II_of_patch_irho_isigma(xp, xm_irho, xm_isigma);
		    Jac.sum_into_element(II, xm_JJ, Jac_rho_sigma);
		    }
	       else {
		    const ghost_zone& xmgz
		       = xp.corner_ghost_zone_containing_point
				(m_irho < 0, m_isigma < 0,
				 xm_irho, xm_isigma);
		    add_ghost_zone_Jacobian(ps, Jac,
					    Jac_rho_sigma,
					    xp, xmgz,
					    II, xm_irho, xm_isigma);
		    }
	    }
	    }
	      }

	if (ps.N_additional_points())
		{
		const int np = ps.N_grid_points();
		Jac.set_element(II,np, 0.0);	// insert dummy value
		}

	}
	}
    }
}
	  }

//******************************************************************************

//
// This function adds the ghost-zone Jacobian dependency contributions
// for a single ghost-zone point, to a Jacobian matrix.
//
// Arguments:
// ps = The patch system.
// Jac = (out) The Jacobian matrix.
// mol = The molecule coefficient.
// xp = The patch containing the center point of the molecule.
// xmgz = If the x+m point is in a ghost zone, this must be that ghost zone.
//	  If the x+m point is not in a ghost zone, this argument is ignored.
// x_II = The Jacobian row of the x point.
// xm_(irho,isigma) = The coordinates (in xp) of the x+m point of the molecule.
//
namespace {
void add_ghost_zone_Jacobian(const patch_system& ps,
			     Jacobian& Jac,
			     fp mol,
			     const patch& xp, const ghost_zone& xmgz,
			     int x_II,
			     int xm_irho, int xm_isigma)
{
const patch_edge& xme = xmgz.my_edge();
const int xm_iperp = xme.iperp_of_irho_isigma(xm_irho, xm_isigma);
const int xm_ipar  = xme. ipar_of_irho_isigma(xm_irho, xm_isigma);

// FIXME: this won't change from one call to another
//        ==> it would be more efficient to reuse the same buffer
//            across multiple calls on this function
int global_min_ym, global_max_ym;
ps.synchronize_Jacobian_global_minmax_ym(global_min_ym, global_max_ym);
jtutil::array1d<fp> Jacobian_buffer(global_min_ym, global_max_ym);

// on what other points y does this molecule point xm depend
// via the patch_system::synchronize() operation?
int y_iperp;
int y_posn, min_ym, max_ym;
const patch_edge& ye = ps.synchronize_Jacobian(xmgz,
					       xm_iperp, xm_ipar,
					       y_iperp,
					       y_posn, min_ym, max_ym,
					       Jacobian_buffer);
patch& yp = ye.my_patch();

// add the Jacobian contributions from the ym points
	for (int ym = min_ym ; ym <= max_ym ; ++ym)
	{
	const int y_ipar = y_posn + ym;
	const int y_irho   = ye.  irho_of_iperp_ipar(y_iperp,y_ipar);
	const int y_isigma = ye.isigma_of_iperp_ipar(y_iperp,y_ipar);
	const int y_JJ = Jac.II_of_patch_irho_isigma(yp, y_irho, y_isigma);
	Jac.sum_into_element(x_II, y_JJ, mol*Jacobian_buffer(ym));
	}
}
	  }

//******************************************************************************

//
// If ps_ptr != NULL and Jac_ptr != NULL, this function sums the d/dr
// terms into the Jacobian matrix of the expansion Theta(h), computing
// those terms by finite differencing.
//
// If ps_ptr == NULL and Jac_ptr == NULL, this function does a dummy
// computation, in which only any expansion() (and hence geometry
// interpolator) calls are done, these with the number of interpolation
// points set to 0 and all the output array pointers set to NULL.
//
// It's illegal for one but not both of ps_ptr and Jac_ptr to be NULL.
//
// The basic algorithm is that
//	Jac += diag[ (Theta(h+epsilon/2) - Theta(h-epsilon/2)) / epsilon ]
//
// Inputs (angular gridfns, on ghosted grid):
//	h			# shape of trial surface
//	Theta			# Theta(h) assumed to already be computed
//				# (saved and restored, but not used)
//
// Outputs:
//	Jac += d/dr terms
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.
//
namespace {
enum expansion_status
  expansion_Jacobian_dr_FD
	(patch_system* ps_ptr, Jacobian* Jac_ptr,
         const struct what_to_compute& compute_info,
	 const struct cactus_grid_info& cgi,
	 const struct geometry_info& gi,
	 const struct Jacobian_info& Jacobian_info,
	 const struct error_info& error_info, bool initial_flag,
	 bool print_msg_flag)
{
const bool active_flag = (ps_ptr != NULL) && (Jac_ptr != NULL);
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   horizon Jacobian: %sd/dr terms (finite diff)",
		   active_flag ? "" : "dummy ");

const fp epsilon = Jacobian_info.perturbation_amplitude;

what_to_compute this_compute_info (compute_info);
this_compute_info.surface_modification = modification_none;
this_compute_info.surface_selection = selection_definition;
this_compute_info.desired_value = 0.0;

fp additional_save_Theta;

// compute Theta(h-epsilon/2)
if (active_flag)
   then {
	ps_ptr->gridfn_copy(gfns::gfn__Theta, gfns::gfn__save_Theta);
	ps_ptr->gridfn_copy(gfns::gfn__mean_curvature, gfns::gfn__save_mean_curvature);
	if (ps_ptr->N_additional_points())
	   then {
		const int np = ps_ptr->N_grid_points();
		additional_save_Theta = ps_ptr->gridfn_data(gfns::gfn__Theta)[np];
		}
	ps_ptr->add_to_ghosted_gridfn(-epsilon/2, gfns::gfn__h);
	}
const
  enum expansion_status status = expansion(ps_ptr,
                                           this_compute_info,
					   cgi, gi,
					   error_info, initial_flag);
if (status != expansion_success)
   then {
        expansion(NULL,
                  this_compute_info,
                  cgi, gi,
                  error_info, false);
        return status;					// *** ERROR RETURN ***
  }

// compute Theta(h+epsilon/2)
if (active_flag)
   then {
	ps_ptr->gridfn_copy(gfns::gfn__Theta, gfns::gfn__old_Theta);
	ps_ptr->add_to_ghosted_gridfn(epsilon, gfns::gfn__h);
	}
const
  enum expansion_status status2 = expansion(ps_ptr,
                                            this_compute_info,
					    cgi, gi,
					    error_info, initial_flag);
if (status2 != expansion_success)
   then return status2;					// *** ERROR RETURN ***

if (active_flag)
   then {
	    for (int pn = 0 ; pn < ps_ptr->N_patches() ; ++pn)
	    {
	    patch& p = ps_ptr->ith_patch(pn);
		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		const int II = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
		const fp old_Theta = p.gridfn(gfns::gfn__old_Theta,
					      irho,isigma);
		const fp new_Theta = p.gridfn(gfns::gfn__Theta,
					      irho,isigma);
		const fp d_dr_term = (new_Theta - old_Theta) / epsilon;
		Jac_ptr->sum_into_element(II,II, d_dr_term);
		}
		}
	    }

	// restore h and Theta
	ps_ptr->add_to_ghosted_gridfn(-epsilon/2, gfns::gfn__h);
	ps_ptr->gridfn_copy(gfns::gfn__save_Theta, gfns::gfn__Theta);
	ps_ptr->gridfn_copy(gfns::gfn__save_mean_curvature, gfns::gfn__mean_curvature);
	if (ps_ptr->N_additional_points())
	   then {
		const int np = ps_ptr->N_grid_points();
		ps_ptr->gridfn_data(gfns::gfn__Theta)[np] = additional_save_Theta;
		}
	}

return expansion_success;				// *** NORMAL RETURN ***
}
	  }

//******************************************************************************

	  }	// namespace AHFinderDirect
