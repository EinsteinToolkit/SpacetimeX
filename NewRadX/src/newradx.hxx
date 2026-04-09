#ifndef NEWRADX_HXX
#define NEWRADX_HXX

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
 * @brief Applies radiative boundary condition to the RHS of a state variable
 * on multipatch grids. See above for original credits.
 *
 * @param cctkGH Pointer to Cactus grid hierarchy struct.
 * @param var State variable which will have boundary conditions applied to it.
 * @param rhs RHS of the evolution equation for @param var
 * @param vcoordx Global x vertex coordinates grid function. Providade by
 * CoordinatesX
 * @param vcoordy Global y vertex coordinates grid function. Providade by
 * CoordinatesX
 * @param vcoordz Global z vertex coordinates grid function. Providade by
 * CoordinatesX
 * @param vJ_da_dx Vertex centered coordinate transformation Jacobian da/dx.
 * Provided by any thorn implementing the Multipatch interface.
 * @param vJ_da_dy Vertex centered coordinate transformation Jacobian da/dy.
 * Provided by any thorn implementing the Multipatch interface.
 * @param vJ_da_dz Vertex centered coordinate transformation Jacobian da/dz.
 * Provided by any thorn implementing the Multipatch interface.
 * @param vJ_db_dx Vertex centered coordinate transformation Jacobian db/dx.
 * Provided by any thorn implementing the Multipatch interface.
 * @param vJ_db_dy Vertex centered coordinate transformation Jacobian db/dy.
 * Provided by any thorn implementing the Multipatch interface.
 * @param vJ_db_dz Vertex centered coordinate transformation Jacobian db/dz.
 * Provided by any thorn implementing the Multipatch interface.
 * @param vJ_dc_dx Vertex centered coordinate transformation Jacobian dc/dx.
 * Provided by any thorn implementing the Multipatch interface.
 * @param vJ_dc_dy Vertex centered coordinate transformation Jacobian dc/dy.
 * Provided by any thorn implementing the Multipatch interface.
 * @param vJ_dc_dz Vertex centered coordinate transformation Jacobian dc/dz.
 * Provided by any thorn implementing the Multipatch interface.
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

#endif // NEWRADX_HXX