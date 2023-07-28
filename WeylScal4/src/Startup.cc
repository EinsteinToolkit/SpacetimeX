/*  File produced by Kranc */

#include "cctk.h"

extern "C" int WeylScal4_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "WeylScal4";
  CCTK_RegisterBanner(banner);
  return 0;
}
