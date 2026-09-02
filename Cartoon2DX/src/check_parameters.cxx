#include <cartoon2dx.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Cartoon2DX {

static const char *carpetx_keyword(const char *name) {
  int type = -1;
  const void *const ptr = CCTK_ParameterGet(name, "CarpetX", &type);
  if (!ptr)
    CCTK_VERROR("Cartoon2DX: CarpetX::%s not found", name);
  return *static_cast<const char *const *>(ptr);
}

extern "C" void Cartoon2DX_CheckParameters(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Cartoon2DX_CheckParameters;
  DECLARE_CCTK_PARAMETERS;

  const char *const boundary_y = carpetx_keyword("boundary_y");
  const char *const boundary_upper_y = carpetx_keyword("boundary_upper_y");
  const char *const boundary_x = carpetx_keyword("boundary_x");
  if (!CCTK_EQUALS(boundary_y, "cartoon") ||
      !CCTK_EQUALS(boundary_upper_y, "cartoon"))
    CCTK_VERROR("Cartoon2DX requires CarpetX::boundary_y and "
                "boundary_upper_y = \"cartoon\"");
  if (fill_negative_x && !CCTK_EQUALS(boundary_x, "cartoon"))
    CCTK_VERROR("Cartoon2DX::fill_negative_x requires CarpetX::boundary_x = "
                "\"cartoon\"");
  if (!fill_negative_x && CCTK_EQUALS(boundary_x, "cartoon"))
    CCTK_VERROR("CarpetX::boundary_x = \"cartoon\" requires "
                "Cartoon2DX::fill_negative_x = yes");

  if (verbose) {
    CCTK_VINFO("Cartoon2DX: interpolation order %d (%d-point stencil)%s", order,
               order + 1,
               fill_negative_x ? ", fill x<0 by rotation from x>=0" : "");
  }
}

} // namespace Cartoon2DX
