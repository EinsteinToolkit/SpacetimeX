// test_patch_system.cc -- test driver for patch_system::
// $Header$

//
// <<<misc global data>>>
// <<<prototypes>>>
//
// driver - Cactus interface
/// test_ghost_zone_Jacobians - test Jacobian of ghost zone synchronization ops
//
/// verify_gpn - verify gpn <--> (patch,irho,isigma) conversions
///
/// setup_sym_fn_xyz - set up symmetrized test function fn(global_[xyz])
/// setup_fn_rho_sigma - set up test function fn(rho,sigma)
/// finite_diff - compute linear combination of finite differences
/// analytic_derivs - compute linear combination of analytic derivatives
///
/// sym_fn_xyz - symmetrized test function fn(x,y,z) + fn(-y,x,z) + ...
/// fn_xyz - test function fn(x,y,z)
///
/// gridfn_minus - compute gridfn x - gridfn y --> gridfn z
/// ghosted_gridfn_minus - compute [ghosted] gridfn x - gridfn y --> gridfn z
///
/// fn_rho_sigma - test function fn(rho,sigma)
/// finite_diff_fn - finite differences of fn(rho,sigma)
/// analytic_deriv_fn - analytical derivs of fn(rho,sigma) (via Maple)
//

#include <stdio.h>
#include <assert.h>
#include <math.h>

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

#include "coords.hh"
#include "grid.hh"
#include "fd_grid.hh"
#include "patch.hh"
#include "patch_edge.hh"
#include "patch_interp.hh"
#include "ghost_zone.hh"
#include "patch_system.hh"

using namespace AHFinderDirect;

//******************************************************************************

//
// misc global data
//

// which test are we going to do?
static const int which_deriv_fn          = 0x1;
static const int which_deriv_rho         = 0x2;
static const int which_deriv_sigma       = 0x4;
static const int which_deriv_rho_rho     = 0x8;
static const int which_deriv_rho_sigma   = 0x10;
static const int which_deriv_sigma_sigma = 0x20;

// when testing multiple derivatives, we combine them together
// with these "quasi-random" weights (chosen from digits of pi)
static const fp deriv_weight_fn          = 3.14,
		deriv_weight_rho         = 1.59,
		deriv_weight_sigma       = 2.65,
		deriv_weight_rho_rho     = 3.58,
		deriv_weight_rho_sigma   = 9.79,
		deriv_weight_sigma_sigma = 3.23;

// n.b. nominal gfns must all be < 0, ghosted > 0
// nominal gridfns
static const int FD_derivs_gfn = -1;
static const int analytic_derivs_gfn = -2;
static const int nominal_error_gfn = -3;
static const int nominal_min_gfn = -3;
static const int nominal_max_gfn = -1;

// ghosted gridfns
static const int test_fn_gfn = 1;
static const int test_fn_copy_gfn = 2;
static const int ghosted_error_gfn = 3;
static const int ghosted_min_gfn = 1;
static const int ghosted_max_gfn = 3;

//******************************************************************************

//
// ***** prototypes *****
//

extern "C"
  void test_patch_system(CCTK_ARGUMENTS);

namespace {
void test_ghost_zone_Jacobians(patch_system& ps,
				int test_gfn, int NP_test_gfn,
				fp perturbation_amplitude,
				bool perturb_all_y_patch_points,
				const char Jacobian_file_name[]);

void verify_gpn(const patch_system& ps);

void setup_sym_fn_xyz(patch_system& ps, int ghosted_gfn, bool want_ghost_zones);
void setup_fn_rho_sigma(patch_system& ps, int ghosted_gfn);
void finite_diff(patch_system& ps,
		 int ghosted_gfn_src, int gfn_dst,
		 int which_derivs);
void analytic_derivs(patch_system& ps, int gfn_dst, int which_derivs);

void gridfn_minus(patch_system& ps, 
		  int gfn_x, int gfn_y, int gfn_dst);
void ghosted_gridfn_minus(patch_system& ps, 
			  int ghosted_gfn_x, int ghosted_gfn_y,
			  int ghosted_gfn_dst);

fp sym_fn_xyz(enum patch_system::patch_system_type type, fp x, fp y, fp z);
fp fn_xyz(fp x, fp y, fp z);

fp fn_rho_sigma(fp rho, fp sigma);
fp finite_diff_fn(const patch& p,
		  int ghosted_gfn_src, int which_derivs,
		  int irho, int isigma);
fp analytic_deriv_fn(fp rho, fp sigma, int which_derivs);
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function is the Cactus interface for the test driver.
//
extern "C"
  void test_patch_system(CCTK_ARGUMENTS)
{
DECLARE_CCTK_ARGUMENTS
DECLARE_CCTK_PARAMETERS

//
// set up the interpatch interpolator
//
CCTK_VInfo(CCTK_THORNSTRING, "setting up interpatch interpolator");
const int interp_handle = CCTK_InterpHandle(interpatch_interpolator_name);
if (interp_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "couldn't find interpolator \"%s\"!",
		   interpatch_interpolator_name);		/*NOTREACHED*/
const int interp_par_table_handle
	= Util_TableCreateFromString(interpatch_interpolator_pars);
if (interp_par_table_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "bad interpatch-interpolator parameter(s) \"%s\"!",
		   interpatch_interpolator_pars);		/*NOTREACHED*/

//
// create the patch system
//
CCTK_VInfo(CCTK_THORNSTRING, "about to create patch system");
patch_system ps(origin_x, origin_y, origin_z,
		patch_system::type_of_name(patch_system_type),
		ghost_zone_width, patch_overlap_width, N_zones_per_right_angle,
		nominal_min_gfn, nominal_max_gfn,
		ghosted_min_gfn, ghosted_max_gfn,
		interp_handle, interp_par_table_handle);
CCTK_VInfo(CCTK_THORNSTRING, "patch system created ok");

//
// do the actual tests
//
if      (STRING_EQUAL(which_test, "gridfn"))
   then {
	verify_gpn(ps);
	setup_sym_fn_xyz(ps, test_fn_gfn, true);
	ps.print_ghosted_gridfn(test_fn_gfn, "test_fn.dat");
	}

else if (STRING_EQUAL(which_test, "read gridfn"))
   then {
	verify_gpn(ps);
	ps.read_ghosted_gridfn(test_fn_gfn, "test_fn.dat");
	ps.print_ghosted_gridfn(test_fn_gfn, "test_fn2.dat");
	}

else if (STRING_EQUAL(which_test, "synchronize"))
   then {
	verify_gpn(ps);
	setup_sym_fn_xyz(ps, test_fn_gfn, false);
	ps.print_ghosted_gridfn(test_fn_gfn, "test_fn_init.dat");

	setup_sym_fn_xyz(ps, test_fn_copy_gfn, true);
	ps.print_ghosted_gridfn(test_fn_copy_gfn, "test_fn_copy.dat");

	ps.synchronize(test_fn_gfn, test_fn_gfn);
	ps.print_ghosted_gridfn(test_fn_gfn, "test_fn_sync.dat");

	ghosted_gridfn_minus(ps,
			     test_fn_gfn, test_fn_copy_gfn,
			     ghosted_error_gfn);
	ps.print_ghosted_gridfn(ghosted_error_gfn, "ghosted_error.dat");
	}

else if (STRING_EQUAL(which_test, "ghost zone Jacobian"))
   then {
	verify_gpn(ps);
	test_ghost_zone_Jacobians(ps,
				   test_fn_gfn, test_fn_copy_gfn,
				   NP_Jacobian__perturbation_amplitude,
				 (NP_Jacobian__perturb_all_y_patch_points != 0),
				   Jacobian_file_name);
	}

else if (STRING_EQUAL(which_test, "derivatives"))
   then {
	verify_gpn(ps);
	setup_fn_rho_sigma(ps, test_fn_gfn);
	ps.print_ghosted_gridfn(test_fn_gfn, "test_fn.dat");
	finite_diff(ps, test_fn_gfn, FD_derivs_gfn, which_derivs);
	ps.print_gridfn(FD_derivs_gfn, "FD_derivs.dat");
	analytic_derivs(ps, analytic_derivs_gfn, which_derivs);
	ps.print_gridfn(analytic_derivs_gfn, "analytic_derivs.dat");
	gridfn_minus(ps,
		     FD_derivs_gfn, analytic_derivs_gfn,
		     nominal_error_gfn);
	ps.print_gridfn(nominal_error_gfn, "nominal_error.dat");
	}

else	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "unknown which_test=\"%s\"!",
		   which_test);					/*NOTREACHED*/

CCTK_VInfo(CCTK_THORNSTRING, "destroying patch system");
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function tests the computation of the Jacobian of the
// ghost_zone::synchronize() operation.  In outline, it does the following:
//
//	set up a test function in test_gfn on the nominal grid 
//	synchronize test_gfn
//	print ths synchronized function
//	set up the same test function in NP_test_gfn on the nominal grid
//		for each patch xp and ghost zone xgz
//		{
//		compute the synchronize() Jacobian for this ghost zone
//		   (via  xgz.compute_Jacobian()  et al)
//			for each point x in (p,xgz)
//			{
//				for each point y in xgz.Jacobian_y_patch()
//				{
//				if ( !perturb_all_y_patch_points
//				     && (y iperp != xgz Jacobian y iperp) )
//				   then continue;	// *** LOOP CONTROL ***
//				const fp save_y_gridfn = NP_test_gfn at y
//				NP_test_gfn at y += perturbation_amplitude
//				synchronize NP_test_gfn
//				NP_Jacobian = (NP_test_gfn at x - test_gfn at x)
//					      / perturbation_amplitude
//				NP_test_gfn at y = save_y_gridfn
//				Jacobian = xgz.Jacobian(...)
//				print Jacobians to output file
//				if (the Jacobians differ)
//				   then {
//					print lots of debugging info to stdout
//					goto done;
//					}
//				}
//			}
//		}
//	done:
//
// Note that this test can be *very* slow to run!  Eg. on a 731 MHz PIII
// laptop, a full-sphere patch system with 5 degree resolution, a ghost zone
// iperp width of 2 points, and perturb_all_y_patch_points set to true,
// took around 25 minutes of cpu time, and produced an 20 MB Jacobian
// output file.  Even with perturb_all_y_patch_points set to false, it
// took 1 minute 15 seconds cpu time, and produced a 1 MB Jacobian output
// file.
//
// Arguments:
// ps = (in out) THe patch system in/on which to do the computations.
// test_gfn, NP_test_gfn = The gfns of two ghosted test gridfns.
// perturbation_amplitude = The perturbation amplitude for the NP Jacobian.
// perturb_all_y_patch_points
//	= true  ==> When computing the NP Jacobian, try perturbing at
//		    every point in the y patch.  This gives 
//	  false ==> Only try perturbing at y_iper == Jacobian_y_iperp.
// Jacobian_file_name = The name of the output file to which both Jacobians
//			should be written.
//
namespace {
void test_ghost_zone_Jacobians(patch_system& ps,
				int test_gfn, int NP_test_gfn,
				fp perturbation_amplitude,
				bool perturb_all_y_patch_points,
				const char Jacobian_file_name[])
{
CCTK_VInfo(CCTK_THORNSTRING, "testing ghost zone synchronize() Jacobian...");

setup_sym_fn_xyz(ps, test_gfn, false);
ps.synchronize(test_fn_gfn, test_fn_gfn);
ps.print_ghosted_gridfn(test_fn_gfn, "test_fn.dat");

setup_sym_fn_xyz(ps, NP_test_gfn, false);

FILE *fileptr = fopen(Jacobian_file_name, "w");
if (fileptr == NULL)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"test_ghost_zone_Jacobians(): can't open output plot file \"%s\"!",
		   Jacobian_file_name);				/*NOTREACHED*/
fprintf(fileptr, "# column  1 = x patch number\n");
fprintf(fileptr, "# column  2 = x ghost_zone is_min()\n");
fprintf(fileptr, "# column  3 = x ghost_zone is_rho()\n");
fprintf(fileptr, "# column  4 = x_iperp\n");
fprintf(fileptr, "# column  5 = x_ipar\n");
fprintf(fileptr, "# column  6 = x_irho\n");
fprintf(fileptr, "# column  7 = x_isigma\n");
fprintf(fileptr, "# column  8 = y patch number\n");
fprintf(fileptr, "# column  9 = y ghost_zone is_min()\n");
fprintf(fileptr, "# column 10 = y ghost_zone is_rho()\n");
fprintf(fileptr, "# column 11 = y_iperp\n");
fprintf(fileptr, "# column 12 = y_ipar\n");
fprintf(fileptr, "# column 13 = y_irho\n");
fprintf(fileptr, "# column 14 = y_isigma\n");
fprintf(fileptr, "# column 15 = Jacobian\n");
fprintf(fileptr, "# column 16 = NP_Jacobian\n");
fprintf(fileptr, "# column 17 = Jacobian error\n");

    //*** for each patch p and ghost zone xgz
    for (int xpn = 0 ; xpn < ps.N_patches() ; ++xpn)
    {
    patch& xp = ps.ith_patch(xpn);

    // n.b. these loops must use _int_ variables for the loop
    //      to terminate!
    for (int want_min = false ; want_min <= true ; ++want_min)
    {
    for (int want_rho = false ; want_rho <= true ; ++want_rho)
    {
    const patch_edge& xe = xp.minmax_ang_patch_edge(want_min, want_rho);
    ghost_zone& xgz      = xp.minmax_ang_ghost_zone(want_min, want_rho);

    patch&            yp = xgz.Jacobian_y_patch();
    const patch_edge& ye = xgz.Jacobian_y_edge();

    CCTK_VInfo(CCTK_THORNSTRING,
	       "   testing x patch %s, edge %s",
	       xp.name(), xe.name());
    CCTK_VInfo(CCTK_THORNSTRING,
	       "           y patch %s, edge %s",
	       yp.name(), ye.name());

    xgz.compute_Jacobian(test_gfn, test_gfn);
    const int Jacobian_min_y_ipar_m = xgz.Jacobian_min_y_ipar_m();
    const int Jacobian_max_y_ipar_m = xgz.Jacobian_max_y_ipar_m();

	//*** for each point x in (p,xgz)
	for (int x_iperp = xgz.min_iperp() ;
	     x_iperp <= xgz.max_iperp() ;
	     ++x_iperp)
	{
	for (int x_ipar = xgz.min_ipar(x_iperp) ;
	     x_ipar <= xgz.max_ipar(x_iperp) ;
	     ++x_ipar)
	{
	const int x_irho   = xe.  irho_of_iperp_ipar(x_iperp, x_ipar);
	const int x_isigma = xe.isigma_of_iperp_ipar(x_iperp, x_ipar);

	const int Jacobian_y_iperp = xgz.Jacobian_y_iperp(x_iperp);
	const int Jacobian_y_ipar_posn
		= xgz.Jacobian_y_ipar_posn(x_iperp, x_ipar);

	    //*** for each point y in gz.Jacobian_patch()
	    for (int y_irho = yp.min_irho() ;
		     y_irho <= yp.max_irho() ;
		     ++y_irho)
	    {
	    for (int y_isigma = yp.min_isigma() ;
		 y_isigma <= yp.max_isigma() ;
		 ++y_isigma)
	    {
	    const int y_iperp = ye.iperp_of_irho_isigma(y_irho,y_isigma);
	    const int y_ipar  = ye. ipar_of_irho_isigma(y_irho,y_isigma);

	    if ( !perturb_all_y_patch_points && (y_iperp != Jacobian_y_iperp) )
	       then continue;				// *** LOOP CONTROL ***

	    // compute the NP Jacobian
	    const fp save_y_gridfn
		= yp.ghosted_gridfn(NP_test_gfn, y_irho,y_isigma);
	    yp.ghosted_gridfn(NP_test_gfn, y_irho,y_isigma)
		+= perturbation_amplitude;

	    ps.synchronize(NP_test_gfn, NP_test_gfn);
	    const fp NP_Jacobian
		= (   xp.ghosted_gridfn(NP_test_gfn, x_irho,x_isigma)
		    - xp.ghosted_gridfn(   test_gfn, x_irho,x_isigma) )
		  / perturbation_amplitude;
	    yp.ghosted_gridfn(NP_test_gfn, y_irho,y_isigma) = save_y_gridfn;

	    // compute the query Jacobian
	    const int y_ipar_m = y_ipar - Jacobian_y_ipar_posn;
	    const bool m_in_molecule
		= (y_iperp == Jacobian_y_iperp)
		  && (y_ipar_m >= Jacobian_min_y_ipar_m)
		  && (y_ipar_m <= Jacobian_max_y_ipar_m);
	    const fp Jacobian = m_in_molecule
				? xgz.Jacobian(x_iperp, x_ipar, y_ipar_m)
				: 0.0;

	    // print the results
	    const fp error = Jacobian - NP_Jacobian;
	    fprintf(fileptr,
"%d %d %d\t%d %d\t%d %d\t%d %d %d\t%d %d\t%d %d\t%.10g\t%.10g\t%e\n",
		    xp.patch_number(), xe.is_min(), xe.is_rho(),
		    x_iperp, x_ipar, x_irho, x_isigma,
		    yp.patch_number(), ye.is_min(), ye.is_rho(),
		    y_iperp, y_ipar, y_irho, y_isigma,
		    double(Jacobian), double(NP_Jacobian),
		    double(error));

	    // debugging code in case the Jacobian is wrong :(
	    if (jtutil::abs(error) > 1.0e-6)
	       then {
		    printf("### large Jacobian error!\n");

		    printf("x: p patch %s, edge %s\n", xp.name(), xe.name());
		    printf("y: q patch %s, edge %s\n", yp.name(), ye.name());

		    const fp x_rho   = xp.rho_of_irho    (x_irho);
		    const fp x_sigma = xp.sigma_of_isigma(x_isigma);
		    const fp y_rho   = yp.rho_of_irho    (y_irho);
		    const fp y_sigma = yp.sigma_of_isigma(y_isigma);

		    const fp x_drho   = jtutil::degrees_of_radians(x_rho);
		    const fp x_dsigma = jtutil::degrees_of_radians(x_sigma);
		    const fp y_drho   = jtutil::degrees_of_radians(y_rho);
		    const fp y_dsigma = jtutil::degrees_of_radians(y_sigma);

		    printf(
"x iperp=%d ipar=%d   irho=%d isigma=%d   drho=%g dsigma=%g\n",
			   x_iperp, x_ipar, x_irho, x_isigma, x_drho, x_dsigma);
		    printf(
"y iperp=%d ipar=%d   irho=%d isigma=%d   drho=%g dsigma=%g\n",
			   y_iperp, y_ipar, y_irho, y_isigma, y_drho, y_dsigma);

		    const fp x_mu  = xp.mu_of_rho_sigma (x_rho, x_sigma);
		    const fp x_nu  = xp.nu_of_rho_sigma (x_rho, x_sigma);
		    const fp x_phi = xp.phi_of_rho_sigma(x_rho, x_sigma);
		    const fp y_mu  = yp.mu_of_rho_sigma (y_rho, y_sigma);
		    const fp y_nu  = yp.nu_of_rho_sigma (y_rho, y_sigma);
		    const fp y_phi = yp.phi_of_rho_sigma(y_rho, y_sigma);

		    const fp x_dmu  = jtutil::degrees_of_radians(x_mu);
		    const fp x_dnu  = jtutil::degrees_of_radians(x_nu);
		    const fp x_dphi = jtutil::degrees_of_radians(x_phi);
		    const fp y_dmu  = jtutil::degrees_of_radians(y_mu);
		    const fp y_dnu  = jtutil::degrees_of_radians(y_nu);
		    const fp y_dphi = jtutil::degrees_of_radians(y_phi);
		
		    printf("x dmu=%g dnu=%g dphi=%g\n", x_dmu, x_dnu, x_dphi);
		    printf("y dmu=%g dnu=%g dphi=%g\n", y_dmu, y_dnu, y_dphi);

		    printf("Jacobian=%.10g\tNP_Jacobian=%.10g\terror=%e\n",
			   double(Jacobian), double(NP_Jacobian),
			   double(error));

		    printf("Jacobian_y_[min,max]_ipar_m=[%d,%d]\n",
			   Jacobian_min_y_ipar_m, Jacobian_max_y_ipar_m);
		    printf("Jacobian_y_iperp=%d Jacobian_y_ipar_posn=%d\n",
			   Jacobian_y_iperp, Jacobian_y_ipar_posn);
		    printf("y_ipar_m=%d\n", y_ipar_m);

		    printf("###\n");
		    goto done;
		    }
	    }
	    }
	}
	}
    }
    }

    }

done:
fclose(fileptr);
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function verifies the nominal- and ghosted-grid
// gpn <--> (patch,irho,isigma) conversions.
//
namespace {
void verify_gpn(const patch_system& ps)
{
CCTK_VInfo(CCTK_THORNSTRING,
	   "testing gpn <--> (patch,irho,isigma) conversions");
int gpn = 0;
int ghosted_gpn = 0;

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	const patch& p = ps.ith_patch(pn);

	// nominal grid
		  {
		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		assert( ps.gpn_of_patch_irho_isigma(p,irho,isigma) == gpn );
		int gpn_irho, gpn_isigma;
		const patch& gpn_p
		   = ps.patch_irho_isigma_of_gpn(gpn, gpn_irho, gpn_isigma);
		assert( gpn_p      == p      );
		assert( gpn_irho   == irho   );
		assert( gpn_isigma == isigma );

		++gpn;
		}
		}
		  }

	// ghosted grid
		  {
		for (int irho = p.ghosted_min_irho() ;
		     irho <= p.ghosted_max_irho() ;
		     ++irho)
		{
		for (int isigma = p.ghosted_min_isigma() ;
		     isigma <= p.ghosted_max_isigma() ;
		     ++isigma)
		{
		assert( ps.ghosted_gpn_of_patch_irho_isigma(p,irho,isigma)
			== ghosted_gpn );
		int gpn_irho, gpn_isigma;
		const patch& gpn_p
		    = ps.ghosted_patch_irho_isigma_of_gpn(ghosted_gpn,
							  gpn_irho, gpn_isigma);
		assert( gpn_p      == p      );
		assert( gpn_irho   == irho   );
		assert( gpn_isigma == isigma );

		++ghosted_gpn;
		}
		}
		  }
	}

assert(         gpn == ps.        N_grid_points() );
assert( ghosted_gpn == ps.ghosted_N_grid_points() );
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function sets up the ghosted test function for the gridfn and
// synchronize tests, symmetrizing the test function to match the
// symmetry of the patch system.
//
// Arguments:
// ps = The patch system.
// ghosted_gfn = Specifies the gridfn to set up.
// want_ghost_zones = true ==> Set up on ghosted part of ghosted grid
//		      false ==> Set up on nominal part of ghosted grid
//
namespace {
void setup_sym_fn_xyz(patch_system& ps, int ghosted_gfn, bool want_ghost_zones)
{
CCTK_VInfo(CCTK_THORNSTRING, "setting up ghosted test fn(x,y,z)");
CCTK_VInfo(CCTK_THORNSTRING, "   on %s part of ghosted grid",
	   want_ghost_zones ? "ghosted" : "nominal");

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.effective_min_irho(want_ghost_zones) ;
		     irho <= p.effective_max_irho(want_ghost_zones) ;
		     ++irho)
		{
		for (int isigma = p.effective_min_isigma(want_ghost_zones) ;
		     isigma <= p.effective_max_isigma(want_ghost_zones) ;
		     ++isigma)
		{
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);

		fp local_x, local_y, local_z;
		p.xyz_of_r_rho_sigma(1.0, rho, sigma,
				     local_x, local_y, local_z);

		p.ghosted_gridfn(ghosted_gfn, irho,isigma)
			= sym_fn_xyz(ps.type(), local_x, local_y, local_z);
		}
		}
	}
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function sets up the test function for the finite differencing
// tests.
//
// Arguments:
// ps = The patch system.
// ghosted_gfn = Specifies the gridfn to set up.
//
namespace {
void setup_fn_rho_sigma(patch_system& ps, int ghosted_gfn)
{
CCTK_VInfo(CCTK_THORNSTRING, "setting up ghosted test fn(rho,sigma)");

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.ghosted_min_irho() ;
		     irho <= p.ghosted_max_irho() ;
		     ++irho)
		{
		for (int isigma = p.ghosted_min_isigma() ;
		     isigma <= p.ghosted_max_isigma() ;
		     ++isigma)
		{
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);

		p.ghosted_gridfn(ghosted_gfn, irho,isigma)
			= fn_rho_sigma(rho, sigma);
		}
		}
	}
}
	  }

//******************************************************************************

//
// This function computes (on the nominal grid only) the specified
// linear combination of finite derivatives of a test function.
//
// Arguments:
// ps = The patch system.
// gfn_src = Specifies the gridfn to finite difference.
// gfn_dst = Specifies the gridfn in which to store the result.
// which = Specifies which finite derivatives to include in the test.
//
namespace {
void finite_diff(patch_system& ps,
		 int ghosted_gfn_src, int gfn_dst,
		 int which_derivs)
{
CCTK_VInfo(CCTK_THORNSTRING,
	   "finite differencing ghosted_gfn=%d --> gfn=%d (which_derivs=%d)",
	   ghosted_gfn_src, gfn_dst, which_derivs);

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

	for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
	{
	for (int isigma = p.min_isigma() ; isigma <= p.max_isigma() ; ++isigma)
	{
	p.gridfn(gfn_dst, irho,isigma)
		= finite_diff_fn(p,
				 ghosted_gfn_src, which_derivs,
				 irho,isigma);
	}
	}

	}
}
	  }

//******************************************************************************

//
// This function computes the specified linear combination of
// derivatives of the test function, everywhere on the nominal grid.
//
namespace {
void analytic_derivs(patch_system& ps,
		     int gfn_dst,
		     int which_derivs)
{
CCTK_VInfo(CCTK_THORNSTRING,
	   "computing analytic derivatives (which_derivs=%d)",
	   which_derivs);

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

	for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
	{
	for (int isigma = p.min_isigma() ; isigma <= p.max_isigma() ; ++isigma)
	{
	const fp rho = p.rho_of_irho(irho);
	const fp sigma = p.sigma_of_isigma(isigma);
	p.gridfn(gfn_dst, irho,isigma)
		= analytic_deriv_fn(rho,sigma, which_derivs);
	}
	}

	}
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function computes [nominal]  gridfn_x - gridfn_y --> gridfn_z .
//
namespace {
void gridfn_minus(patch_system& ps, 
		  int gfn_x, int gfn_y, int gfn_dst)
{
CCTK_VInfo(CCTK_THORNSTRING,
	   "[nominal] gfn=%d - gfn=%d --> gfn=%d",
	   gfn_x, gfn_y, gfn_dst);

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

	for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
	{
	for (int isigma = p.min_isigma() ; isigma <= p.max_isigma() ; ++isigma)
	{
	p.gridfn(gfn_dst, irho,isigma)
		=   p.gridfn(gfn_x, irho,isigma)
		  - p.gridfn(gfn_y, irho,isigma);
	}
	}

	}
}
	  }

//******************************************************************************

//
// This function computes [ghosted]  gridfn_x - gridfn_y --> gridfn_z .
//
namespace {
void ghosted_gridfn_minus(patch_system& ps, 
			  int ghosted_gfn_x, int ghosted_gfn_y,
			  int ghosted_gfn_dst)
{
CCTK_VInfo(CCTK_THORNSTRING,
	   "[ghosted] gfn=%d - gfn=%d --> gfn=%d",
	   ghosted_gfn_x, ghosted_gfn_y, ghosted_gfn_dst);

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.ghosted_min_irho() ;
		     irho <= p.ghosted_max_irho() ;
		     ++irho)
		{
		for (int isigma = p.ghosted_min_isigma() ;
		     isigma <= p.ghosted_max_isigma() ;
		     ++isigma)
		{
		p.ghosted_gridfn(ghosted_gfn_dst, irho,isigma)
			=   p.ghosted_gridfn(ghosted_gfn_x, irho,isigma)
			  - p.ghosted_gridfn(ghosted_gfn_y, irho,isigma);
		}
		}

	}
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function symmetrizes  fn_xyz()  (about the origin) to match
// the patch system's symmetries.
//
// To rotate f(x,y) by 90, 180, or 270 degrees:
//
//	      (-y,x) |
//	             |
//	             |
//	             |
//	             |        (x,y)
//	             |
//	-------------+-------------
//	             |
//	(-x,-y)      |
//	             |
//	             |
//	             |
//	             | (y,-x)
//
namespace {
fp sym_fn_xyz(enum patch_system::patch_system_type type, fp x, fp y, fp z)
{
switch	(type)
	{
case patch_system::full_sphere_patch_system:
	return fn_xyz(x,y,z);
	break;
case patch_system::plus_z_hemisphere_patch_system:
	return fn_xyz(x,y,+z) + fn_xyz(x,y,-z);
	break;
case patch_system::plus_xy_quadrant_patch_system:
	return   fn_xyz(+x,+y,z)
	       + fn_xyz(-y,+x,z)
	       + fn_xyz(-x,-y,z)
	       + fn_xyz(+y,-x,z);
	break;
case patch_system::plus_xz_quadrant_patch_system:
	return   fn_xyz(+x,+y,+z) + fn_xyz(-x,-y,+z)
	       + fn_xyz(+x,+y,-z) + fn_xyz(-x,-y,-z);
	break;
case patch_system::plus_xyz_octant_patch_system:
	return   fn_xyz(+x,+y,+z) + fn_xyz(+x,+y,-z)
	       + fn_xyz(-y,+x,+z) + fn_xyz(-y,+x,-z)
	       + fn_xyz(-x,-y,+z) + fn_xyz(-x,-y,-z)
	       + fn_xyz(+y,-x,+z) + fn_xyz(+y,-x,-z);
	break;
default:
	error_exit(PANIC_EXIT,
"***** sym_fn_xyz(): impossible type=(int)%d!\n",
		   int(type));					/*NOTREACHED*/
			}
}
	  }

//******************************************************************************

//
// This is the underlying test function for our function and ghost-zone tests.
//
namespace {
fp fn_xyz(fp x, fp y, fp z)
{
return (x*(x+0.238) + 2.417*y*(y-0.917) + 1.38*z*(z-0.472))
       * tanh(jtutil::pow3(cos(z)));
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This is the underlying test function for our finite differencing
// tests.
//
namespace {
fp fn_rho_sigma(fp rho, fp sigma)
{
return exp(sin(1.38*rho)) * tanh(0.17+0.83*jtutil::pow2(sin(sigma)));
}
	  }

//******************************************************************************

//
// This function computes the sum of our linear combination of various
// finite difference approximations to derivatives of fn_rho_sigma().
//
// Arguments:
// ghosted_gfn_src = Specifies the gridfn to finite difference.
//
namespace {
fp finite_diff_fn(const patch& p,
		  int ghosted_gfn_src, int which_derivs,
		  int irho, int isigma)
{
fp sum = 0.0;

if (which_derivs & which_deriv_fn)
   then sum += deriv_weight_fn
	       * p.ghosted_gridfn(ghosted_gfn_src, irho,isigma);

if (which_derivs & which_deriv_rho)
   then sum += deriv_weight_rho
	       * p.partial_rho(ghosted_gfn_src, irho,isigma);
if (which_derivs & which_deriv_sigma)
   then sum += deriv_weight_sigma
	       * p.partial_sigma(ghosted_gfn_src, irho,isigma);

if (which_derivs & which_deriv_rho_rho)
   then sum += deriv_weight_rho_rho
	       * p.partial_rho_rho(ghosted_gfn_src, irho,isigma);
if (which_derivs & which_deriv_rho_sigma)
   then sum += deriv_weight_rho_sigma
	       * p.partial_rho_sigma(ghosted_gfn_src, irho,isigma);
if (which_derivs & which_deriv_sigma_sigma)
   then sum += deriv_weight_sigma_sigma
	       * p.partial_sigma_sigma(ghosted_gfn_src, irho,isigma);

return sum;
}
	  }

//******************************************************************************

//
// This function computes the sum of a specified linear combination of
// various (analytical) derivatives of fn_rho_sigma().
//
// The derivatives were machine-generated via Maple's codegen[C]() function:
// "deriv_patch_system.maple" is the Maple input
// "deriv_patch_system.out" is the Maple input; code here is cut-n-pasted
//			   from the Maple codegen[C]() output there
//
namespace {
fp analytic_deriv_fn(fp rho, fp sigma, int which_derivs)
{
fp sum = 0.0;

if (which_derivs & which_deriv_fn)
   then sum += deriv_weight_fn
	       * fn_rho_sigma(rho,sigma);

if (which_derivs & which_deriv_rho)
   then sum += deriv_weight_rho
	       * 0.138E1*cos(0.138E1*rho)*exp(sin(0.138E1*rho))
		 *tanh(0.17+0.83*pow(sin(sigma),2.0));
if (which_derivs & which_deriv_sigma)
   then sum += deriv_weight_sigma
	       * 0.166E1*exp(sin(0.138E1*rho))
		 *(1.0-pow(tanh(0.17+0.83*pow(sin(sigma),2.0)),2.0))
		 *sin(sigma)*cos(sigma);

if (which_derivs & which_deriv_rho_rho)
   then sum += deriv_weight_rho_rho
	       * (
		 -0.19044E1*sin(0.138E1*rho)*exp(sin(0.138E1*rho))
		  *tanh(0.17+0.83*pow(sin(sigma),2.0))
		 +0.19044E1*pow(cos(0.138E1*rho),2.0)*exp(sin(0.138E1*rho))
		  *tanh(0.17+0.83*pow(sin(sigma),2.0))
		 );
if (which_derivs & which_deriv_rho_sigma)
   then sum += deriv_weight_rho_sigma
	       * 0.22908E1*cos(0.138E1*rho)*exp(sin(0.138E1*rho))
		 *(1.0-pow(tanh(0.17+0.83*pow(sin(sigma),2.0)),2.0))
		 *sin(sigma)*cos(sigma);
if (which_derivs & which_deriv_sigma_sigma)
   then sum += deriv_weight_sigma_sigma
	       * (
		 -0.55112E1*exp(sin(0.138E1*rho))
		  *tanh(0.17+0.83*pow(sin(sigma),2.0))
		  *(1.0-pow(tanh(0.17+0.83*pow(sin(sigma),2.0)),2.0))
		  *pow(sin(sigma),2.0)*pow(cos(sigma),2.0)
		 +0.166E1*exp(sin(0.138E1*rho))
		  *(1.0-pow(tanh(0.17+0.83*pow(sin(sigma),2.0)),2.0))
		  *pow(cos(sigma),2.0)
		 -0.166E1*exp(sin(0.138E1*rho))
		  *(1.0-pow(tanh(0.17+0.83*pow(sin(sigma),2.0)),2.0))
		  *pow(sin(sigma),2.0)
		 );

return sum;
}
	  }
