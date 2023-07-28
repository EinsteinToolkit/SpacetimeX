#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



void SphericalSurface_CheckState (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_SphericalSurface_CheckState;
  DECLARE_CCTK_PARAMETERS;
  
  int n;
  
  for (n=0; n<nsurfaces; ++n) {
    
    if (sf_valid[n] > 0 && sf_active[n] == 0) {
      
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Surface #%d has sf_valid set to a positive value, but does "
                  "not have sf_active set.  This is an error in the thorn "
                  "which calculated this surface", n);
      
    }
    
  }
}
