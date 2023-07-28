// misc-driver.cc -- misc support routines
// $Header$
//
// Cactus_gridfn_varindex - get Cactus gridfn variable index
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
// This function gets the Cactus variable index of a given gridfn, and
// checks to make sure this is valid (i.e. that it's not an error code).
//
int Cactus_gridfn_varindex(const char gridfn_name[])
{
const int varindex = CCTK_VarIndex(gridfn_name);
if (varindex < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   Cactus_gridfn_varindex(): error return from CCTK_VarIndex()\n"
"                             for Cactus gridfn!\n"
"                             name=\"%s\" status=%d"
		   ,
		   gridfn_name, varindex);			/*NOTREACHED*/

return varindex;
}

//******************************************************************************

	  }	// namespace AHFinderDirect
