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

  simd<CCTK_REAL> (*calc_deriv_upwind)(
      const GF3D2<const CCTK_REAL> &, const simdl<CCTK_REAL> &,
      const vect<int, dim> &, const vect<CCTK_REAL, dim> &,
      const vec<simd<CCTK_REAL>, dim> &);

  simd<CCTK_REAL> (*calc_diss)(const GF3D2<const CCTK_REAL> &,
                               const simdl<CCTK_REAL> &, const vect<int, dim> &,
                               const vect<CCTK_REAL, dim> &);

  switch (deriv_order) {
  case 2: {
    calc_deriv_upwind = &Derivs::calc_deriv_upwind<2>;
    calc_diss = &Derivs::calc_diss<2>;
    break;
  }
  case 4: {
    calc_deriv_upwind = &Derivs::calc_deriv_upwind<4>;
    calc_diss = &Derivs::calc_diss<4>;
    break;
  }
  // case 6: {
  //   calc_deriv_upwind = &Derivs::calc_deriv_upwind<6>;
  //   calc_diss = &Derivs::calc_diss<6>;
  //   break;
  // }
  default:
    assert(0);
  }

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
          const vreal rhs_new =
              rhs_old + calc_deriv_upwind(gf_, mask, p.I, dx, betaG);
          gf_rhs_.store(mask, p.I, rhs_new);
        });

  } else {

    const Loop::GridDescBaseDevice grid(cctkGH);
    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          const vec<vreal, dim> betaG = gf_betaG_(mask, p.I);
          const vreal rhs_old = gf_rhs_(mask, p.I);
          const vreal rhs_new = rhs_old +
                                calc_deriv_upwind(gf_, mask, p.I, dx, betaG) +
                                epsdiss * calc_diss(gf_, mask, p.I, dx);
          gf_rhs_.store(mask, p.I, rhs_new);
        });
  }
}

} // namespace Z4co

#endif // #ifndef Z4CO_DERIVS_HXX
