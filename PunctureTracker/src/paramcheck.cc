#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

extern "C" void PunctureTracker_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  // do this only on processor 0
  if (CCTK_MyProc(cctkGH) != 0)
    return;

  // for (int i = 0; i < 10; ++i) {
  //   if (which_surface_to_store_info[i] != -1) {
  //     if (which_surface_to_store_info[i] >= nsurfaces)
  //       CCTK_PARAMWARN("You assigned a greater surface index than there are "
  //                      "spherical surfaces!");
  //   }
  // }
}
