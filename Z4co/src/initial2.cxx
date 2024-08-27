#include <derivs.hxx>
#include <loop_device.hxx>
#include <simd.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#ifdef __CUDACC__
#include <nvToolsExt.h>
#endif

namespace Z4co {
using namespace Arith;
using namespace Loop;
using namespace std;

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr vec<smat<T, D>, D>
calc_gammal(const smat<vec<T, D>, D> &dg) {
  // Gammal_abc
  return vec<smat<T, D>, D>([&](int a) ARITH_INLINE {
    return smat<T, D>([&](int b, int c) ARITH_INLINE {
      return (dg(a, b)(c) + dg(a, c)(b) - dg(b, c)(a)) / 2;
    });
  });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr vec<smat<T, D>, D>
calc_gamma(const smat<T, D> &gu, const vec<smat<T, D>, D> &Gammal) {
  // Gamma^a_bc
  return vec<smat<T, D>, D>([&](int a) ARITH_INLINE {
    return smat<T, D>([&](int b, int c) ARITH_INLINE {
      return sum<D>([&](int x)
                        ARITH_INLINE { return gu(a, x) * Gammal(x)(b, c); });
    });
  });
}

extern "C" void Z4co_Initial2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4co_Initial2;
  DECLARE_CCTK_PARAMETERS;

  const vect<CCTK_REAL, 3> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  vec<simd<CCTK_REAL>, dim> (*calc_deriv)(
      const GF3D2<const CCTK_REAL> &, const simdl<CCTK_REAL> &,
      const vect<int, dim> &, const vect<CCTK_REAL, dim> &);

  switch (deriv_order) {
  case 2:
    calc_deriv = &Derivs::calc_deriv<2>;
    break;
  case 4:
    calc_deriv = &Derivs::calc_deriv<4>;
    break;
  case 6:
    calc_deriv = &Derivs::calc_deriv<6>;
    break;
  default:
    assert(0);
  }

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_gammat1{
      GF3D2<const CCTK_REAL>(layout1, gammatxx),
      GF3D2<const CCTK_REAL>(layout1, gammatxy),
      GF3D2<const CCTK_REAL>(layout1, gammatxz),
      GF3D2<const CCTK_REAL>(layout1, gammatyy),
      GF3D2<const CCTK_REAL>(layout1, gammatyz),
      GF3D2<const CCTK_REAL>(layout1, gammatzz)};

  const vec<GF3D2<CCTK_REAL>, 3> gf_Gamt1{GF3D2<CCTK_REAL>(layout1, Gamtx),
                                          GF3D2<CCTK_REAL>(layout1, Gamty),
                                          GF3D2<CCTK_REAL>(layout1, Gamtz)};

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4co_Initial2::initial2");
#endif
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);

        // Load
        const smat<vreal, 3> gammat = gf_gammat1(mask, index1);

        // Calculate Z4c variables (only Gamt)
        const smat<vreal, 3> gammatu = calc_inv(gammat, vreal(1));

        const smat<vec<vreal, 3>, 3> dgammat([&](int a, int b) {
          return calc_deriv(gf_gammat1(a, b), mask, p.I, dx);
        });

        const vec<smat<vreal, 3>, 3> Gammatl = calc_gammal(dgammat);
        const vec<smat<vreal, 3>, 3> Gammat = calc_gamma(gammatu, Gammatl);
        const vec<vreal, 3> Gamt([&](int a) ARITH_INLINE {
          return sum_symm<3>([&](int x, int y) ARITH_INLINE {
            return gammatu(x, y) * Gammat(a)(x, y);
          });
        });

        // Store
        gf_Gamt1.store(mask, index1, Gamt);
      });
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace Z4co
