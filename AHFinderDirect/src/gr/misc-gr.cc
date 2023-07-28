// misc-gr.cc -- misc support routines
// $Header$
//
// expansion_status_string - string describing expansion_status
//

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

//
// This function decodes the  Jacobian_compute_method  parameter (string) into
// an internal enum for future use.
//
enum Jacobian_compute_method
  decode_Jacobian_compute_method(const char Jacobian_compute_method_string[])
{
if	(STRING_EQUAL(Jacobian_compute_method_string,
		      "numerical perturbation"))
   then return Jacobian__numerical_perturbation;
else if (STRING_EQUAL(Jacobian_compute_method_string,
		      "symbolic differentiation with finite diff d/dr"))
   then return Jacobian__symbolic_diff_with_FD_dr;
else if (STRING_EQUAL(Jacobian_compute_method_string,
		      "symbolic differentiation"))
   then return Jacobian__symbolic_diff;
else	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   decode_Jacobian_compute_method():\n"
"        unknown Jacobian_compute_method_string=\"%s\"!",
		   Jacobian_compute_method_string);		/*NOTREACHED*/
}

//******************************************************************************

//
// This function returns (a pointer to) a C-style string describing
// an  expansion_status  value.
//
const char* expansion_status_string(enum expansion_status status)
{
switch	(status)
	{
case expansion_success:
	return "success";
	break;
case expansion_failure__surface_nonfinite:
	return "infinity/NaN in surface shape!";
	break;
case expansion_failure__surface_too_large:
	return "surface too large";
	break;
case expansion_failure__surface_outside_grid:
	return "surface outside grid";
	break;
case expansion_failure__surface_in_excised_region:
	return "surface in excised region";
	break;
case expansion_failure__geometry_nonfinite:
	return "infinity/NaN in 3-geometry!";
	break;
case expansion_failure__gij_not_positive_definite:
	return "g_ij not positive definite!";
	break;
default:
	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
		    "expansion_status_string(): unknown status=(int)%d!",
		    status);					/*NOTREACHED*/
	}
}

//******************************************************************************

	  }	// namespace AHFinderDirect
