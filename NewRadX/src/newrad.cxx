#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cmath>
#include "newrad.hxx"
#include "loop.hxx"
#include "loop_device.hxx"

namespace NewRadX {

using namespace Loop;
using namespace std;


namespace {
template <typename T> constexpr T pow2(const T x) { return x * x; }
} // namespace

// Adapted from NewRad thorn by E. Schnetter, used with Carpet.
// Original code adapted from BSSN_MoL's files NewRad.F and newrad.h.
// This code was probably originally written by Miguel Alcubierre.
void newrad(const cGH *restrict const cctkGH,
            const Loop::GF3D2<const CCTK_REAL> var, Loop::GF3D2<CCTK_REAL> rhs,
            const CCTK_REAL var0, //!< value at infinity
            const CCTK_REAL v0,    //!< propagation speed
            const CCTK_REAL radpower //!< exponent in radial fall-off
) {
  DECLARE_CCTK_ARGUMENTS;

  constexpr vect<int, dim> DI{1, 0, 0};
  constexpr vect<int, dim> DJ{0, 1, 0};
  constexpr vect<int, dim> DK{0, 0, 1};

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const auto derivx =
      [=] CCTK_DEVICE CCTK_HOST(const GF3D2<const CCTK_REAL> &u_, const PointDesc &p)
			    CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto I = p.I;
        if (p.NI[0] == 0)
          // interior
          // in direction parallel to face/edge, apply symmetric stencil
          // currently second-order accurate
          // TODO: apply same finite-difference order as interior
          return (u_(I + DI) - u_(I - DI)) / (2 * dx);
        if (p.NI[0] == +1)
          // upper boundary
          // in direction normal to face/edge/corner, apply asymmetric stencil
          // currently second-order accurate
          // TODO: apply same finite-difference order as interior
          return +(3 * u_(I) - 4 * u_(I - DI) + u_(I - 2 * DI)) / (2 * dx);
        if (p.NI[0] == -1)
          // lower boundary
          // in direction normal to face/edge/corner, apply asymmetric stencil
          // currently second-order accurate
          // TODO: apply same finite-difference order as interior
          return -(3 * u_(I) - 4 * u_(I + DI) + u_(I + 2 * DI)) / (2 * dx);
        assert(0);
      };

  const auto derivy =
      [=] CCTK_DEVICE CCTK_HOST(const GF3D2<const CCTK_REAL> &u_, const PointDesc &p)
			    CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto I = p.I;
        if (p.NI[1] == 0)
          // interior
          // in direction parallel to face/edge, apply symmetric stencil
          // currently second-order accurate
          // TODO: apply same finite-difference order as interior
          return (u_(I + DJ) - u_(I - DJ)) / (2 * dy);
        if (p.NI[1] == +1)
          // upper boundary
          // in direction normal to face/edge/corner, apply asymmetric stencil
          // currently second-order accurate
          // TODO: apply same finite-difference order as interior
          return +(3 * u_(I) - 4 * u_(I - DJ) + u_(I - 2 * DJ)) / (2 * dy);
        if (p.NI[1] == -1)
          // lower boundary
          // in direction normal to face/edge/corner, apply asymmetric stencil
          // currently second-order accurate
          // TODO: apply same finite-difference order as interior
          return -(3 * u_(I) - 4 * u_(I + DJ) + u_(I + 2 * DJ)) / (2 * dy);
        assert(0);
      };

  const auto derivz =
      [=] CCTK_DEVICE CCTK_HOST(const GF3D2<const CCTK_REAL> &u_, const PointDesc &p)
			    CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto I = p.I;
        if (p.NI[2] == 0)
          // interior
          // in direction parallel to face/edge, apply symmetric stencil
          // currently second-order accurate
          // TODO: apply same finite-difference order as interior
          return (u_(I + DK) - u_(I - DK)) / (2 * dz);
        if (p.NI[2] == +1)
          // upper boundary
          // in direction normal to face/edge/corner, apply asymmetric stencil
          // currently second-order accurate
          // TODO: apply same finite-difference order as interior
          return +(3 * u_(I) - 4 * u_(I - DK) + u_(I - 2 * DK)) / (2 * dz);
        if (p.NI[2] == -1)
          // lower boundary
          // in direction normal to face/edge/corner, apply asymmetric stencil
          // currently second-order accurate
          // TODO: apply same finite-difference order as interior
          return -(3 * u_(I) - 4 * u_(I + DK) + u_(I + 2 * DK)) / (2 * dz);
        assert(0);
      };

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_radbnd<0,0,0>(
    grid.nghostzones,
    [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
        CCTK_ATTRIBUTE_ALWAYS_INLINE {
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
    const CCTK_REAL r = sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z));

    // Find local wave speeds at radiative boundary point p.I
    const CCTK_REAL vx = v0 * p.x / r;
    const CCTK_REAL vy = v0 * p.y / r;
    const CCTK_REAL vz = v0 * p.z / r;
    //CCTK_REAL const vr = sqrt(pow2(vx) + pow2(vy) + pow2(vz));

    // Derivatives
    const CCTK_REAL varx = derivx(var, p);
    const CCTK_REAL vary = derivy(var, p);
    const CCTK_REAL varz = derivz(var, p);

    // radiative rhs
    rhs(p.I) = -vx * varx - vy * vary - vz * varz - v0 * (var(p.I) - var0) / r;

    if(radpower>0){
      // When solution is known to have a Coulomb-like component asymptotically,
      // this is estimated and extrapolated from physical interior points to
      // the radiative boundary.
      // I.e.:
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
      // the Coulomb term, d_t h / r ^n , is estimated by extrapolating from the
      // nearest interior points
      //
      // Displacement to get from p.I to interior point placed nghostpoints away
      vect<int, dim> displacement{grid.nghostzones[0]*p.NI[0],
                                  grid.nghostzones[1]*p.NI[1],
                                  grid.nghostzones[2]*p.NI[2]};
      vect<int, dim> intp = p.I - displacement;

      // coordinates at p.I-displacement
      const CCTK_REAL xint = p.x - displacement[0]*p.DX[0];
      const CCTK_REAL yint = p.y - displacement[1]*p.DX[1];
      const CCTK_REAL zint = p.z - displacement[2]*p.DX[2];
      const CCTK_REAL rint = sqrt(pow2(xint) + pow2(yint) + pow2(zint));

      // Find local wave speeds at physical point p.I-displacement
      const CCTK_REAL vxint = v0 * xint / rint;
      const CCTK_REAL vyint = v0 * yint / rint;
      const CCTK_REAL vzint = v0 * zint / rint;

      // Derivatives at physical point p.I-displacement
      const CCTK_REAL varxint = (var(intp + p.DI[0]) - var(intp - p.DI[0])) / (2 * p.DX[0]);
      const CCTK_REAL varyint = (var(intp + p.DI[1]) - var(intp - p.DI[1])) / (2 * p.DX[1]);
      const CCTK_REAL varzint = (var(intp + p.DI[2]) - var(intp - p.DI[2])) / (2 * p.DX[2]);

      // extrapolate Coulomb component, rescale to account for radial fall-off
      CCTK_REAL rad =
              - vxint * varxint
              - vyint * varyint
              - vzint * varzint
              - v0 * (var(intp) - var0) / rint;
      CCTK_REAL aux = (rhs(intp) - rad) * pow(rint/r,radpower);

      // radiative rhs with extrapolated Coulomb correction
      rhs(p.I) += aux;
    }

  });

}

} // namespace NewRadX
