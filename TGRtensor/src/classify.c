/* $Header$ */

#include <math.h>

#include "cctk.h"

int CCTK_FCALL CCTK_FNAME(TAT_isnan) (const CCTK_REAL * restrict const x)
{
#ifdef HAVE_ISNAN
  return isnan(*x);
#else
  return 0;
#endif
}

int CCTK_FCALL CCTK_FNAME(TAT_isinf) (const CCTK_REAL * restrict const x)
{
#ifdef HAVE_ISINF
  return isinf(*x);
#else
  return 0;
#endif
}

int CCTK_FCALL CCTK_FNAME(TAT_finite) (const CCTK_REAL * restrict const x)
{
#ifdef HAVE_ISFINITE
  return isfinite(*x);
#else
  return 1;
#endif
}

CCTK_REAL CCTK_FCALL CCTK_FNAME(TAT_nan) (void)
{
  return 0.0 / 0.0;
}

CCTK_REAL CCTK_FCALL CCTK_FNAME(TAT_inf) (void)
{
  return 1.0 / 0.0;
}
