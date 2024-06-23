#include <stdio.h>

#include "interpolate.hxx"

static void report_interp_error(int ierr) {
  if (ierr < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "CCTK_InterpGridArrays returned error code %d", ierr);
  }
}

void Multipole_Interp(CCTK_ARGUMENTS, CCTK_REAL xs[], CCTK_REAL ys[],
                      CCTK_REAL zs[], int real_idx, int imag_idx,
                      CCTK_REAL sphere_real[], CCTK_REAL sphere_imag[]) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_INT N_dims = 3;
  const CCTK_INT npoints =
      CCTK_MyProc(cctkGH) == 0 ? (ntheta + 1) * (nphi + 1) : 0;
  const void *interp_coords[N_dims] = {(const void *)xs, (const void *)ys,
                                       (const void *)zs};

  CCTK_INT N_input_arrays = imag_idx == -1 ? 1 : 2;
  CCTK_INT N_output_arrays = imag_idx == -1 ? 1 : 2;
  const CCTK_INT input_array_indices[2] = {real_idx, imag_idx};
  // Interpolation result
  CCTK_POINTER output_arrays[2];
  output_arrays[0] = sphere_real;
  output_arrays[1] = sphere_imag;

  int ierr = -1;
  /* DriverInterpolate arguments that aren't currently used */
  const int coord_system_handle = 0;
  CCTK_INT const interp_coords_type_code = 0;
  CCTK_INT const output_array_types[1] = {0};
  int interp_handle = CCTK_InterpHandle("CarpetX");

  if (interp_handle < 0)
    CCTK_VERROR("Could not obtain inteprolator handle for built-in 'CarpetX' "
                "interpolator: %d",
                interp_handle);

  // Interpolation parameter table
  CCTK_INT operations[1][N_dims];
  for (int var = 0; var < npoints; var++) {
    operations[0][var] = 0;
  }
  int operands[1][N_dims];
  for (int var = 0; var < npoints; var++) {
    operands[0][var] = var;
  }
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (param_table_handle < 0)
    CCTK_VERROR("Can't create parameter table: %d", param_table_handle);
  if ((ierr = Util_TableSetInt(param_table_handle, 1, "order")) < 0)
    CCTK_VERROR("Can't set order in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, N_input_arrays,
                                    (int const *const)operands,
                                    "operand_indices")) < 0)
    CCTK_VERROR("Can't set operand_indices array in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, N_input_arrays,
                                    (int const *const)operations,
                                    "operation_codes")) < 0)
    CCTK_VERROR("Can't set operation_codes array in parameter table: %d", ierr);

  ierr = DriverInterpolate(
      cctkGH, N_dims, interp_handle, param_table_handle, coord_system_handle,
      npoints, interp_coords_type_code, interp_coords, N_input_arrays,
      input_array_indices, N_output_arrays, output_array_types, output_arrays);

  report_interp_error(ierr);

  if (imag_idx == -1) {
    for (int i = 0; i < (ntheta + 1) * (nphi + 1); i++) {
      sphere_imag[i] = 0;
    }
  }

  Util_TableDestroy(param_table_handle);
}
