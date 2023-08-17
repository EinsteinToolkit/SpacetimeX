// initial_guess.cc -- set up the initial guess
// $Header$
//
// <<<access to persistent data>>>
// <<<prototypes for functions local to this file>>>
// setup_initial_guess - set up initial guess in h
// decode_initial_guess_method - decode the  initial_guess_method  parameter
/// setup_Kerr_horizon - set up Kerr horizon in h (Kerr or Kerr-Schild coords)
/// setup_coord_ellipsoid - setup up a coordinate ellipsoid in h
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
// ***** prototypes for functions local to this file *****
//

namespace {
void setup_Kerr_horizon(patch_system& ps,
			fp x_posn, fp y_posn, fp z_posn,
			fp m, fp a,
			bool Kerr_Schild_flag,
			const struct verbose_info& verbose_info);
void setup_coord_ellipsoid(patch_system& ps,
			   fp x_center, fp y_center, fp z_center,
			   fp x_radius, fp y_radius, fp z_radius,
			   bool print_msg_flag);
	  }

//******************************************************************************

//
// This function sets up the initial guess for a single apparent horizon,
// in the  h  angular gridfn.
//
void setup_initial_guess(patch_system& ps,
			 const struct initial_guess_info& igi,
			 const struct IO_info& IO_info,
			 int hn, int N_horizons,
			 const struct verbose_info& verbose_info)
{
if (verbose_info.print_algorithm_highlights)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "setting initial guess for horizon %d/%d",
		   hn, N_horizons);

switch	(igi.method)
	{
case initial_guess__read_from_named_file:
	input_gridfn__explicit_name(ps, gfns::gfn__h,
				    IO_info,
				    igi.read_from_named_file_info.file_name,
				    verbose_info.print_algorithm_highlights);
	break;

case initial_guess__read_from_h_file:
	input_gridfn(ps, gfns::gfn__h,
		     IO_info, IO_info.h_base_file_name, IO_info.h_min_digits,
		     hn, verbose_info.print_algorithm_highlights);
	break;

case initial_guess__Kerr_Kerr:
	setup_Kerr_horizon(ps,
			   igi.Kerr_Kerr_info.x_posn,
			   igi.Kerr_Kerr_info.y_posn,
			   igi.Kerr_Kerr_info.z_posn,
			   igi.Kerr_Kerr_info.mass,
			   igi.Kerr_Kerr_info.spin,
			   false,		// Kerr coordinates
			   verbose_info);
	break;

case initial_guess__Kerr_KerrSchild:
	setup_Kerr_horizon(ps,
			   igi.Kerr_KerrSchild_info.x_posn,
			   igi.Kerr_KerrSchild_info.y_posn,
			   igi.Kerr_KerrSchild_info.z_posn,
			   igi.Kerr_KerrSchild_info.mass,
			   igi.Kerr_KerrSchild_info.spin,
			   true,		// Kerr-Schild coordinates
			   verbose_info);
	break;

case initial_guess__coord_sphere:
	setup_coord_ellipsoid(ps,
			      igi.coord_sphere_info.x_center,
			      igi.coord_sphere_info.y_center,
			      igi.coord_sphere_info.z_center,
			      igi.coord_sphere_info.radius,
			      igi.coord_sphere_info.radius,
			      igi.coord_sphere_info.radius,
			      verbose_info.print_algorithm_highlights);
	break;

case initial_guess__coord_ellipsoid:
	setup_coord_ellipsoid(ps,
			      igi.coord_ellipsoid_info.x_center,
			      igi.coord_ellipsoid_info.y_center,
			      igi.coord_ellipsoid_info.z_center,
			      igi.coord_ellipsoid_info.x_radius,
			      igi.coord_ellipsoid_info.y_radius,
			      igi.coord_ellipsoid_info.z_radius,
			      verbose_info.print_algorithm_highlights);
	break;

default:
	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "unknown initial guess method=(int)%d!",
		   int(igi.method));				/*NOTREACHED*/
	}
}

//******************************************************************************

//
// This function decodes the  initial_guess_method  parameter (string)
// into an internal enum for future use.
//
enum initial_guess_method
  decode_initial_guess_method(const char initial_guess_method_string[])
{
if	(STRING_EQUAL(initial_guess_method_string, "read from named file"))
   then return initial_guess__read_from_named_file;
else if (STRING_EQUAL(initial_guess_method_string, "read from h file"))
   then return initial_guess__read_from_h_file;
else if (STRING_EQUAL(initial_guess_method_string, "Kerr/Kerr"))
   then return initial_guess__Kerr_Kerr;
else if (STRING_EQUAL(initial_guess_method_string, "Kerr/Kerr-Schild"))
   then return initial_guess__Kerr_KerrSchild;
else if (STRING_EQUAL(initial_guess_method_string, "coordinate sphere"))
   then return initial_guess__coord_sphere;
else if (STRING_EQUAL(initial_guess_method_string, "coordinate ellipsoid"))
   then return initial_guess__coord_ellipsoid;
else	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   decode_initial_guess_method():\n"
"        unknown initial_guess_method_string=\"%s\"!",
		   initial_guess_method_string);		/*NOTREACHED*/
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function sets up the horizon of a Kerr black hole in Kerr or
// Kerr-Schild coordinates, on the nominal grid, in the  h  gridfn.
//
// Kerr-Schild coordinates are described in MTW Exercise 33.8, page 903,
// and the horizon is worked out on page 13.2 of my AHFinderDirect notes.
//
// Arguments:
// [xyz]_posn = The global-coordinate position of the Kerr black hole.
// (m,a) = Describe the Kerr black hole.  Note that my convention has
//	   a=J/m^2 dimensionless, while MTW take a=J/m=m*(my a).
// Kerr_Schild_flag = false to use Kerr coordinates,
//		      true to use Kerr-Schild coordinates
// verbose_info = controls what messages we print
//
namespace {
void setup_Kerr_horizon(patch_system& ps,
			fp x_posn, fp y_posn, fp z_posn,
			fp m, fp a,
			bool Kerr_Schild_flag,
			const struct verbose_info& verbose_info)
{
const char* const name = Kerr_Schild_flag ? "Kerr-Schild" : "Kerr";

if (verbose_info.print_algorithm_highlights)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "   setting Kerr/%s horizon shape",
		   name);
	CCTK_VInfo(CCTK_THORNSTRING,
		   "           posn=(%g,%g,%g) mass=%g spin=J/m^2=%g",
		   double(x_posn), double(y_posn), double(z_posn),
		   double(m), double(a));
	}

// horizon in Kerr coordinates is coordinate sphere
const fp r = m * (1.0 + sqrt(1.0 - a*a));

// horizon in Kerr-Schild coordinates is coordinate ellipsoid
const fp  z_radius = r;
const fp xy_radius = Kerr_Schild_flag ? r * sqrt(1.0 + a*a*m*m/(r*r)) : r;

if (verbose_info.print_algorithm_details)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   setting coordinate %s",
		   Kerr_Schild_flag ? "ellipsoid" : "sphere");
setup_coord_ellipsoid(ps,
		      x_posn, y_posn, z_posn,
		      xy_radius, xy_radius, z_radius,
		      verbose_info.print_algorithm_details);
}
	  }

//******************************************************************************

//
// This function sets up an ellipsoid initial guess, using the formulas
// in "ellipsoid.maple" and the Maple-generated C code in "ellipsoid.c":
//
// ellipsoid has global-coordinates center (A,B,C), radius (a,b,c)
// angular coordinate system has center (U,V,W)
//
// direction cosines wrt angular coordinate center are (xcos,ycos,zcos)
// i.e. a point has coordinates (U+xcos*r, V+ycos*r, W+zcos*r)
//
// then the equation of the ellipsoid is
//	(U+xcos*r - A)^2     (V+ycos*r - B)^2     (W+zcos*r - C)^2
//	-----------------  +  ----------------  +  -----------------  =  1
//	        a^2                  b^2                   c^2
//
// to solve this, we introduce intermediate variables
//	AU = A - U
//	BV = B - V
//	CW = C - W
//
namespace {
void setup_coord_ellipsoid(patch_system& ps,
			   fp x_center, fp y_center, fp z_center,
			   fp x_radius, fp y_radius, fp z_radius,
			   bool print_msg_flag)
{
if (print_msg_flag)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "   setting ellipsoid: center=(%g,%g,%g)",
		   double(x_center), double(y_center), double(z_center));
	CCTK_VInfo(CCTK_THORNSTRING,
		   "                      radius=(%g,%g,%g)",
		   double(x_radius), double(y_radius), double(z_radius));
	}

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);
		fp xcos, ycos, zcos;
		p.xyzcos_of_rho_sigma(rho,sigma, xcos,ycos,zcos);

		// set up variables used by Maple-generated code
		const fp AU = x_center - ps.origin_x();
		const fp BV = y_center - ps.origin_y();
		const fp CW = z_center - ps.origin_z();
		const fp a = x_radius;
		const fp b = y_radius;
		const fp c = z_radius;

		// compute the solutions r_plus and r_minus
		fp r_plus, r_minus;
		#include "ellipsoid.c"

		// exactly one of the solutions (call it r) should be positive
		fp r;
		if      ((r_plus > 0.0) && (r_minus < 0.0))
		   then r = r_plus;
		else if ((r_plus < 0.0) && (r_minus > 0.0))
		   then r = r_minus;
		else    CCTK_VWarn(FATAL_ERROR,
				   __LINE__, __FILE__, CCTK_THORNSTRING,
				   "\n"
"   setup_coord_ellipsoid():\n"
"        expected exactly one r>0 solution to quadratic, got 0 or 2!\n"
"        %s patch (irho,isigma)=(%d,%d) ==> (rho,sigma)=(%g,%g)\n"
"        direction cosines (xcos,ycos,zcos)=(%g,%g,%g)\n"
"        r_plus=%g r_minus=%g\n"
"        ==> this probably means the initial guess surface doesn't contain\n"
"            the local origin point, or more generally that the initial\n"
"            guess surface isn't a Strahlkoerper (\"star-shaped region\")\n"
"            with respect to the local origin point\n"
				   ,
				   p.name(), irho, isigma,
				   double(rho), double(sigma),
				   double(xcos), double(ycos), double(zcos),
				   double(r_plus), double(r_minus));
		   						/*NOTREACHED*/

		// r = horizon radius at this grid point
		p.ghosted_gridfn(gfns::gfn__h, irho,isigma) = r;
		}
		}
	}
}
	  }

//******************************************************************************

	  }	// namespace AHFinderDirect
