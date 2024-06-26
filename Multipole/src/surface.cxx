/* surface.cxx */
/* (c) Liwei Ji 06/2024 */

#include "surface.hxx"

#include <loop_device.hxx>

#include <util_Table.h>

namespace Multipole {

void Surface::interp(CCTK_ARGUMENTS, int real_idx, int imag_idx) {
  DECLARE_CCTK_PARAMETERS;

  const CCTK_INT npoints =
      CCTK_MyProc(cctkGH) == 0 ? (nTheta_ + 1) * (nPhi_ + 1) : 0;

  const void *interp_coords[Loop::dim] = {(const void *)x_.data(),
                                          (const void *)y_.data(),
                                          (const void *)z_.data()};

  CCTK_INT N_input_arrays = imag_idx == -1 ? 1 : 2;
  CCTK_INT N_output_arrays = imag_idx == -1 ? 1 : 2;

  const CCTK_INT input_array_indices[2] = {real_idx, imag_idx};

  // Interpolation result
  CCTK_POINTER output_arrays[2];
  output_arrays[0] = real_.data();
  output_arrays[1] = imag_.data();

  /* DriverInterpolate arguments that aren't currently used */
  const int coord_system_handle = 0;
  CCTK_INT const interp_coords_type_code = 0;
  CCTK_INT const output_array_types[1] = {0};

  int interp_handle = CCTK_InterpHandle("CarpetX");
  if (interp_handle < 0) {
    CCTK_VERROR("Could not obtain inteprolator handle for built-in 'CarpetX' "
                "interpolator: %d",
                interp_handle);
  }

  // Interpolation parameter table
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);

  int ierr = Util_TableSetFromString(param_table_handle, interpolator_pars);

  ierr = DriverInterpolate(
      cctkGH, Loop::dim, interp_handle, param_table_handle, coord_system_handle,
      npoints, interp_coords_type_code, interp_coords, N_input_arrays,
      input_array_indices, N_output_arrays, output_array_types, output_arrays);

  if (ierr < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "CCTK_InterpGridArrays returned error code %d", ierr);
  }

  if (imag_idx == -1) {
    for (int i = 0; i < (nTheta_ + 1) * (nPhi_ + 1); i++) {
      imag_[i] = 0;
    }
  }

  Util_TableDestroy(param_table_handle);
}

} // namespace Multipole
