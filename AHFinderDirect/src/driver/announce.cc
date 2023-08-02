// announce.cc -- annnounce apparent horizon info to other thorns
// $Header$
//
// <<<access to persistent data>>>
// AHFinderDirect_announce - top-level driver for announce stuff
//

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "util_Table.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

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
// This function is called by the Cactus scheduler, to announce any
// desired apparent horizon info to any other thorns that may be interested.
// At present the only info we announce is the centroid position of a
// single selected apparent horizon; if the SetAHCentroid() aliased
// function has been defined then we announce by calling that.
//
extern "C"
  void AHFinderDirect_announce(CCTK_ARGUMENTS)
{
DECLARE_CCTK_ARGUMENTS_AHFinderDirect_announce
DECLARE_CCTK_PARAMETERS

const struct verbose_info& verbose_info = state.verbose_info;

// which horizon to announce?
const int hn = which_horizon_to_announce_centroid;
if (hn == 0)
   then return;						// *** NO-OP RETURN ***

if (! ((hn >= 1) && (hn <= N_horizons)) )
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   AHFinderDirect_announce():\n"
"        invalid horizon number %d to announce\n"
"        (valid range is [1,N_horizons=%d])!\n"
		   ,
		   hn, int(N_horizons));			/*NOTREACHED*/

assert(state.AH_data_array[hn] != NULL);
const struct AH_data& AH_data = *state.AH_data_array[hn];

// only try to announce AH info if we've found AHs at this time level
if (! AH_data.search_flag)
   then return;						// *** NO-OP RETURN ***

// did we actually *find* this horizon?
if (! AH_data.found_flag)
   then return;						// *** NO-OP RETURN ***

// is there anyone to announce it to?
if (CCTK_IsFunctionAliased("SetDriftCorrectPosition"))
   then {
	const struct BH_diagnostics& BH_diagnostics = AH_data.BH_diagnostics;
	const CCTK_REAL xx = BH_diagnostics.centroid_x;
	const CCTK_REAL yy = BH_diagnostics.centroid_y;
	const CCTK_REAL zz = BH_diagnostics.centroid_z;
	if (verbose_info.print_physics_details)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "horizon %d centroid (%g,%g,%g) --> DriftCorrect",
			   hn, double(xx), double(yy), double(zz));
	SetDriftCorrectPosition(cctkGH, xx, yy, zz);
	}
}

//******************************************************************************

//
// This function is called by the Cactus scheduler, to copy any
// desired apparent horizon info to Cactus variables.
//
extern "C"
  void AHFinderDirect_store(CCTK_ARGUMENTS)
{
DECLARE_CCTK_ARGUMENTS_AHFinderDirect_store
DECLARE_CCTK_PARAMETERS

for (int hn = 1; hn <= N_horizons; ++ hn)
  {

  // Store in spherical surface
  const int sn = sf_IdFromName(which_surface_to_store_info[hn], 
                               which_surface_to_store_info_by_name[hn]);
  if (sn == -1)
    then continue;

  if (sn < 0 || sn >= nsurfaces)
    then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   AHFinderDirect_store():\n"
"        invalid surface number %d for horizon number %d\n"
"        (valid range is [0,nsurfaces-1=%d])!\n"
		   ,
		   sn, hn,
                   int(nsurfaces-1));			/*NOTREACHED*/
  
  const struct AH_data& AH_data = *state.AH_data_array[hn];
  const struct BH_diagnostics& BH_diagnostics = AH_data.BH_diagnostics;
  BH_diagnostics.store(cctkGH, hn, sn);

  }
}

//******************************************************************************

//
// This function is called by the Cactus scheduler, to copy any
// desired apparent horizon info to Cactus variables.
//
extern "C"
  void AHFinderDirect_save(CCTK_ARGUMENTS)
{
DECLARE_CCTK_ARGUMENTS_AHFinderDirect_save
DECLARE_CCTK_PARAMETERS

for (int hn = 1; hn <= N_horizons; ++ hn)
  {

  const struct AH_data& AH_data = *state.AH_data_array[hn];
  const struct BH_diagnostics& BH_diagnostics = AH_data.BH_diagnostics;

  // Save in grid array
  BH_diagnostics.save(cctkGH, hn);

  }
}

//******************************************************************************

//
// This function is called by the Cactus scheduler, to copy any
// desired apparent horizon info from Cactus variables.
//
extern "C"
  void AHFinderDirect_recover(CCTK_ARGUMENTS)
{
DECLARE_CCTK_ARGUMENTS_AHFinderDirect_recover
DECLARE_CCTK_PARAMETERS

for (int hn = 1; hn <= N_horizons; ++ hn)
  {

  struct AH_data& AH_data = *state.AH_data_array[hn];
  struct BH_diagnostics& BH_diagnostics = AH_data.BH_diagnostics;

  // Load from grid array
  BH_diagnostics.load(cctkGH, hn);

  }
}

//******************************************************************************

	  }	// namespace AHFinderDirect
