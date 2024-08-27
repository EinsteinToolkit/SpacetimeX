#ifndef Z4CO_DERIVS_HXX
#define Z4CO_DERIVS_HXX

#include <derivs.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Z4co {
using namespace Arith;
using namespace Loop;

constexpr int deriv_o = 4;

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
apply_upwind_diss(const cGH *restrict const cctkGH, const GF3D2<const T> &gf_,
                  const vec<GF3D2<const T>, dim> &gf_betaG_,
                  const GF3D2<T> &gf_rhs_) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const vect<CCTK_REAL, dim> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  if (epsdiss == 0) {

    const Loop::GridDescBaseDevice grid(cctkGH);
    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          const vec<vreal, dim> betaG = gf_betaG_(mask, p.I);
          const vreal rhs_old = gf_rhs_(mask, p.I);
          const vreal rhs_new = rhs_old + Derivs::calc_deriv_upwind<deriv_o>(
                                              gf_, mask, p.I, dx, betaG);
          gf_rhs_.store(mask, p.I, rhs_new);
        });

  } else {

    const Loop::GridDescBaseDevice grid(cctkGH);
    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          const vec<vreal, dim> betaG = gf_betaG_(mask, p.I);
          const vreal rhs_old = gf_rhs_(mask, p.I);
          const vreal rhs_new =
              rhs_old +
              Derivs::calc_deriv_upwind<deriv_o>(gf_, mask, p.I, dx, betaG) +
              epsdiss * Derivs::calc_diss<deriv_o>(gf_, mask, p.I, dx);
          gf_rhs_.store(mask, p.I, rhs_new);
        });
  }
}

} // namespace Z4co

#endif // #ifndef Z4CO_DERIVS_HXX
