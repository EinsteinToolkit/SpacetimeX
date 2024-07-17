#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include <defs.hxx>
#include <mat.hxx>
#include <sum.hxx>
#include <vec.hxx>

namespace Z4c {
using namespace Arith;

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<vec<T, D>, D>
calc_dgu(const smat<T, D> &gu, const smat<vec<T, D>, D> &dg) {
  // g_xy = g_xa g_yb g^ab
  // g_xy,c = (g_xa g_yb g^ab),c
  //        = g_xa,c g_yb g^ab + g_xa g_yb,c g^ab + g_xa g_yb g^ab,c
  // g_xa g_yb g^ab,c = g_xy,c - g_xa,c g_yb g^ab - g_xa g_yb,c g^ab
  //                  = g_xy,c - g_xy,c - g_xy,c
  //                  = - g_xy,c
  // g^ab,c = - g^ax g^by g_xy,c
  return smat<vec<T, D>, D>([&](int a, int b) ARITH_INLINE {
    return vec<T, D>([&](int c) ARITH_INLINE {
      // return sum2sym([&](int x, int y) ARITH_INLINE {
      //   return -gu(a, x) * gu(b, y) * dg(x, y)(c);
      // });
      return sum<D>([&](int x) ARITH_INLINE {
        return -gu(a, x) * sum<D>([&](int y) ARITH_INLINE {
          return gu(b, y) * dg(x, y)(c);
        });
      });
    });
  });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<vec<T, D>, D>
calc_dAu(const smat<T, D> &gu, const smat<vec<T, D>, D> &dgu,
         const smat<T, D> &A, const smat<vec<T, D>, D> &dA) {
  // A^ab,c = (g^ax g^by A_xy),c
  //        = g^ax,c g^by A_xy + g^ax g^by,c A_xy + g^ax g^by A_xy,c
  return smat<vec<T, D>, D>([&](int a, int b) ARITH_INLINE {
    return vec<T, D>([&](int c) ARITH_INLINE {
      // return sum2sym([&](int x, int y) ARITH_INLINE {
      //   return dgu(a, x)(c) * gu(b, y) * A(x, y)   //
      //          + gu(a, x) * dgu(b, y)(c) * A(x, y) //
      //          + gu(a, x) * gu(b, y) * dA(x, y)(c);
      // });
      return sum<D>([&](int x, int y) ARITH_INLINE {
        return dgu(a, x)(c) * gu(b, y) * A(x, y)   //
               + gu(a, x) * dgu(b, y)(c) * A(x, y) //
               + gu(a, x) * gu(b, y) * dA(x, y)(c);
      });
    });
  });
}

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

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr vec<vec<vec<T, D>, D>, D>
calc_gammalu(const smat<T, D> &gu, const vec<smat<T, D>, D> &Gammal) {
  // Gamma_ab^c
  return vec<vec<vec<T, D>, D>, D>([&](int a) ARITH_INLINE {
    return vec<vec<T, D>, D>([&](int b) ARITH_INLINE {
      return vec<T, D>([&](int c) ARITH_INLINE {
        return sum<D>([&](int x)
                          ARITH_INLINE { return Gammal(a)(b, x) * gu(x, c); });
      });
    });
  });
}

} // namespace Z4c

#endif // #ifndef PHYSICS_HXX
