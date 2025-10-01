#include <cctk.h>
#include <loop_device.hxx>

namespace NewRadX {

/**
 * @brief Applies radiative boundary condition to the RHS of a state variable
 *
 * Adapted from NewRad thorn by E. Schnetter, used with Carpet. Original code
 * adapted from BSSN_MoL's files NewRad.F and newrad.h. This code was probably
 * originally written by Miguel Alcubierre.
 *
 * @param cctkGH Pointer to Cactus grid hierarchy struct.
 * @param var State variable which will have boundary conditions applied to it.
 * @param rhs RHS of the evolution equation for @param var
 * @param var0 Value at infinity.
 * @param v0 Propagation speed.
 * @param radpower Radial fall-off exponent
 */
void NewRadX_Apply(const cGH *restrict const cctkGH,
                   const Loop::GF3D2<const CCTK_REAL> &var,
                   const Loop::GF3D2<CCTK_REAL> &rhs, const CCTK_REAL var0,
                   const CCTK_REAL v0, const CCTK_REAL radpower);

#define NEWRADX_MULTIPATCH_QUANTITIES                                          \
  vcoordx, vcoordy, vcoordz, vJ_da_dx, vJ_da_dy, vJ_da_dz, vJ_db_dx, vJ_db_dy, \
      vJ_db_dz, vJ_dc_dx, vJ_dc_dy, vJ_dc_dz

/**
 * @brief Applies radiative boundary condition to the RHS of a state variable.
 * Assumes that:
 *   1. Using multiple patches
 *   2. Patch 0 is cartesian
 *   3. Patches != 0 are spherical-like
 *   4. The local c coordinate is radial and points outward
 *
 * @param cctkGH Pointer to Cactus grid hierarchy struct.
 * @param var State variable which will have boundary conditions applied to it.
 * @param rhs RHS of the evolution equation for @param var
 * @param vcoordx x coordinates grid function.
 * @param vcoordy y coordinates grid function.
 * @param vcoordz z coordinates grid function.
 * @param var0 Value at infinity.
 * @param v0 Propagation speed.
 * @param radpower Radial fall-off exponent
 */
void NewRadX_Apply(const cGH *restrict const cctkGH,
                   const Loop::GF3D2<const CCTK_REAL> &var,
                   const Loop::GF3D2<CCTK_REAL> &rhs,
                   const Loop::GF3D2<const CCTK_REAL> &vcoordx,
                   const Loop::GF3D2<const CCTK_REAL> &vcoordy,
                   const Loop::GF3D2<const CCTK_REAL> &vcoordz,
                   const Loop::GF3D2<const CCTK_REAL> &vJ_da_dx,
                   const Loop::GF3D2<const CCTK_REAL> &vJ_da_dy,
                   const Loop::GF3D2<const CCTK_REAL> &vJ_da_dz,
                   const Loop::GF3D2<const CCTK_REAL> &vJ_db_dx,
                   const Loop::GF3D2<const CCTK_REAL> &vJ_db_dy,
                   const Loop::GF3D2<const CCTK_REAL> &vJ_db_dz,
                   const Loop::GF3D2<const CCTK_REAL> &vJ_dc_dx,
                   const Loop::GF3D2<const CCTK_REAL> &vJ_dc_dy,
                   const Loop::GF3D2<const CCTK_REAL> &vJ_dc_dz,
                   const CCTK_REAL var0, const CCTK_REAL v0,
                   const CCTK_REAL radpower);

} // namespace NewRadX
