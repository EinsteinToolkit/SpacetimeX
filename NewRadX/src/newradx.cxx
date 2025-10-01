#include <cctk.h>

#include <loop_device.hxx>
#include <driver.hxx>

#include <global_derivatives.hxx>

#include "newradx.hxx"

#include <cmath>

namespace NewRadX {

using namespace Loop;
using namespace std;

template <typename T> static inline constexpr T pow2(const T x) {
  return x * x;
}

template <std::size_t dir>
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
c2o(const Loop::PointDesc &p, const vect<int, dim> &pI,
    const Loop::GF3D2<const CCTK_REAL> &gf) noexcept {
  const auto num{gf(pI + p.DI[dir]) - gf(pI - p.DI[dir])};
  const auto den{1.0 / (2.0 * p.DX[dir])};
  return num * den;
}

template <std::size_t dir>
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
l2o(const Loop::PointDesc &p, const vect<int, dim> &pI,
    const Loop::GF3D2<const CCTK_REAL> &gf) noexcept {
  const auto num{-3.0 * gf(pI) + 4.0 * gf(pI + p.DI[dir]) -
                 gf(pI + 2 * p.DI[dir])};
  const auto den{1.0 / (2.0 * p.DX[dir])};
  return num * den;
}

template <std::size_t dir>
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
r2o(const Loop::PointDesc &p, const vect<int, dim> &pI,
    const Loop::GF3D2<const CCTK_REAL> &gf) noexcept {
  const auto num{3.0 * gf(pI) - 4.0 * gf(pI - p.DI[dir]) +
                 gf(pI - 2 * p.DI[dir])};
  const auto den{1.0 / (2.0 * p.DX[dir])};
  return num * den;
}

template <std::size_t dir>
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL calc_deriv(
    const Loop::PointDesc &p, const Loop::GF3D2<const CCTK_REAL> &gf) noexcept {
  if (p.NI[dir] == 0) {
    /* interior
     * in direction parallel to face/edge, apply symmetric stencil
     * currently second-order accurate
     * TODO: apply same finite-difference order as interior
     */
    return c2o<dir>(p, p.I, gf);

  } else if (p.NI[dir] == +1) {
    /* upper boundary
     * in direction normal to face/edge/corner, apply asymmetric stencil
     * currently second-order accurate
     * TODO: apply same finite-difference order as interior
     */
    return r2o<dir>(p, p.I, gf);

  } else if (p.NI[dir] == -1) {
    /* lower boundary
     * in direction normal to face/edge/corner, apply asymmetric stencil
     * currently second-order accurate
     * TODO: apply same finite-difference order as interior
     */
    return l2o<dir>(p, p.I, gf);

  } else {
    assert(0);
  }
}

void NewRadX_Apply(const cGH *restrict const cctkGH,
                   const Loop::GF3D2<const CCTK_REAL> &var,
                   const Loop::GF3D2<CCTK_REAL> &rhs, const CCTK_REAL var0,
                   const CCTK_REAL v0, const CCTK_REAL radpower) {
  DECLARE_CCTK_ARGUMENTS;

  const auto symmetries = CarpetX::ghext->patchdata.at(cctk_patch).symmetries;
  const vect<vect<bool, Loop::dim>, 2> is_sym_bnd{
      {symmetries[0][0] != CarpetX::symmetry_t::none,
       symmetries[0][1] != CarpetX::symmetry_t::none,
       symmetries[0][2] != CarpetX::symmetry_t::none},
      {symmetries[1][0] != CarpetX::symmetry_t::none,
       symmetries[1][1] != CarpetX::symmetry_t::none,
       symmetries[1][2] != CarpetX::symmetry_t::none}};

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_outermost_int_device<0, 0, 0>(
      grid.nghostzones, is_sym_bnd,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // The main part of the boundary condition assumes that we have an
        // outgoing radial wave with some speed v0:
        //
        //    var  =  var0 + u(r-v0*t)/r
        //
        // This implies the following differential equation:
        //
        //    d_t var  =  - v^i d_i var  -  v0 (var - var0) / r
        //
        // where  vi = v0 xi/r

        // coordinate radius at p.I
        const auto r = sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z));

        // Find local wave speeds at radiative boundary point p.I
        const auto vx = v0 * p.x / r;
        const auto vy = v0 * p.y / r;
        const auto vz = v0 * p.z / r;

        // Derivatives
        const auto varx = calc_deriv<0>(p, var);
        const auto vary = calc_deriv<1>(p, var);
        const auto varz = calc_deriv<2>(p, var);

        // radiative rhs
        rhs(p.I) =
            -vx * varx - vy * vary - vz * varz - v0 * (var(p.I) - var0) / r;

        if (radpower > 0) {
          // When solution is known to have a Coulomb-like component
          // asymptotically, this is estimated and extrapolated from
          // physical interior points to the radiative boundary. I.e.:
          //
          //    var  =  var0 + u(r-v0*t)/r + h(t)/r^n
          //
          // This implies the following differential equation:
          //
          //    d_t var  =  - v^i d_i var  -  v0 (var - var0) / r
          //             + v0 (1 - n) h / r^n+1 + d_t h / r^n
          //             ~  - v^i d_i var  -  v0 (var - var0) / r
          //             + d_t h / r^n + O(1/r^n+1)
          //
          // where  vi = v0 xi/r
          //
          // the Coulomb term, d_t h / r ^n , is estimated by extrapolating
          // from the nearest interior points
          //
          // Displacement to get from p.I to interior point placed
          // nghostpoints away
          const vect<int, dim> displacement{grid.nghostzones[0] * p.NI[0],
                                            grid.nghostzones[1] * p.NI[1],
                                            grid.nghostzones[2] * p.NI[2]};
          const vect<int, dim> intp = p.I - displacement;

          assert(intp[0] >= grid.nghostzones[0]);
          assert(intp[1] >= grid.nghostzones[1]);
          assert(intp[2] >= grid.nghostzones[2]);
          assert(intp[0] <= grid.lsh[0] - grid.nghostzones[0] - 1);
          assert(intp[1] <= grid.lsh[1] - grid.nghostzones[1] - 1);
          assert(intp[2] <= grid.lsh[2] - grid.nghostzones[2] - 1);

          // coordinates at p.I-displacement
          const auto xint = p.x - displacement[0] * p.DX[0];
          const auto yint = p.y - displacement[1] * p.DX[1];
          const auto zint = p.z - displacement[2] * p.DX[2];
          const auto rint = sqrt(pow2(xint) + pow2(yint) + pow2(zint));

          // Find local wave speeds at physical point p.I-displacement
          const auto vxint = v0 * xint / rint;
          const auto vyint = v0 * yint / rint;
          const auto vzint = v0 * zint / rint;

          // Derivatives at physical point p.I-displacement
          const auto varxint = c2o<0>(p, intp, var);
          const auto varyint = c2o<1>(p, intp, var);
          const auto varzint = c2o<2>(p, intp, var);

          // Eextrapolate Coulomb component, rescale to account for radial
          // fall-off
          const auto rad = -vxint * varxint - vyint * varyint -
                           vzint * varzint - v0 * (var(intp) - var0) / rint;
          const auto aux = (rhs(intp) - rad) * pow(rint / r, radpower);

          // Radiative rhs with extrapolated Coulomb correction
          rhs(p.I) += aux;
        }
      });
}

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
                   const CCTK_REAL radpower) {
  using namespace CapyrX::MultiPatch::GlobalDerivatives;

  DECLARE_CCTK_ARGUMENTS;

  const auto symmetries = CarpetX::ghext->patchdata.at(cctk_patch).symmetries;
  const vect<vect<bool, Loop::dim>, 2> is_sym_bnd{
      {symmetries[0][0] != CarpetX::symmetry_t::none,
       symmetries[0][1] != CarpetX::symmetry_t::none,
       symmetries[0][2] != CarpetX::symmetry_t::none},
      {symmetries[1][0] != CarpetX::symmetry_t::none,
       symmetries[1][1] != CarpetX::symmetry_t::none,
       symmetries[1][2] != CarpetX::symmetry_t::none}};

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_outermost_int_device<0, 0, 0>(
      grid.nghostzones, is_sym_bnd,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Make sure we are always on the outer boundary of a patch system
        assert(p.patch != 0);
        assert(p.NI[2] >= 0);

        // Find local wave speeds at radiative boundary point p.I
        const auto x = vcoordx(p.I);
        const auto y = vcoordy(p.I);
        const auto z = vcoordz(p.I);
        const auto r = sqrt(pow2(x) + pow2(y) + pow2(z));

        const auto vx = v0 * vcoordx(p.I) / r;
        const auto vy = v0 * vcoordy(p.I) / r;
        const auto vz = v0 * vcoordz(p.I) / r;

        // Local derivatives
        const LocalFirstDerivatives l_dvar{.da = r2o<0>(p, p.I, var),
                                           .db = r2o<1>(p, p.I, var),
                                           .dc = r2o<2>(p, p.I, var)};

        // Global derivatives
        const Jacobians jac{VERTEX_JACOBIANS(p)};
        const auto g_dvar{project_first(l_dvar, jac)};

        // radiative rhs
        rhs(p.I) = -vx * g_dvar.dx - vy * g_dvar.dy - vz * g_dvar.dz -
                   v0 * (var(p.I) - var0) / r;

        if (radpower > 0.0) {
          // TODO
        }
      });
}

} // namespace NewRadX
