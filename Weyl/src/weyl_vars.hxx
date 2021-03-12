#ifndef WEYL_HXX
#define WEYL_HXX

#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"

#include <cplx.hxx>

#include <cmath>

namespace Weyl {

template <typename T> struct weyl_vars_noderivs {

  // Position
  const vec4<T, UP> coord;

  // ADM variables
  const mat3<T, DN, DN> gamma;
  const T alpha;
  const vec3<T, UP> beta;

  // Intermediate quantities
  const vec3<T, DN> betal;

  // 4-metric
  const mat4<T, DN, DN> g;

  // Inverse 4-metric
  const T detg;
  const mat4<T, UP, UP> gu;

  // Tetrad
  const vec4<T, UP> et, ephi, etheta, er;
  const vec4<T, UP> l, n;
  const vec4<cplx<T>, UP> m;

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST weyl_vars_noderivs(
      const T time, const vec3<T, UP> &coord3, const mat3<T, DN, DN> &gamma,
      const T &alpha, const vec3<T, UP> &beta)
      : coord(time, coord3(0), coord3(1), coord3(2)),
        //
        gamma(gamma), alpha(alpha), beta(beta),
        //
        betal([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
            return gamma(a, x) * beta(x);
          });
        }),
        //
        g([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          if (a == 0 && b == 0)
            return -pow(alpha, 2);
          if (a == 0)
            return betal(b - 1);
          if (b == 0)
            return betal(a - 1);
          return gamma(a - 1, b - 1);
        }),
        detg(g.det()),                       //
        gu(g.inv(detg)),                     //
        et(calc_et(gu)),                     //
        ephi(calc_ephi(coord, g)),           //
        etheta(calc_etheta(coord, g, ephi)), //
        er(calc_er(coord, g, etheta, ephi)), //
        l([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return (et(a) + er(a)) / sqrt(T(2));
        }),
        n([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return (et(a) - er(a)) / sqrt(T(2));
        }),
        m([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return cplx<T>(etheta(a), ephi(a)) / sqrt(T(2));
        })
  //
  {}
};

template <typename T> struct weyl_vars : weyl_vars_noderivs<T> {

  // C++ is tedious:

  // Position
  using weyl_vars_noderivs<T>::coord;

  // ADM variables
  using weyl_vars_noderivs<T>::alpha;
  using weyl_vars_noderivs<T>::beta;
  using weyl_vars_noderivs<T>::gamma;

  // Intermediate quantities
  using weyl_vars_noderivs<T>::betal;

  // 4-metric
  using weyl_vars_noderivs<T>::g;

  // Inverse 4-metric
  using weyl_vars_noderivs<T>::detg;
  using weyl_vars_noderivs<T>::gu;

  // Tetrad
  using weyl_vars_noderivs<T>::et;
  using weyl_vars_noderivs<T>::ephi;
  using weyl_vars_noderivs<T>::etheta;
  using weyl_vars_noderivs<T>::er;
  using weyl_vars_noderivs<T>::l;
  using weyl_vars_noderivs<T>::n;
  using weyl_vars_noderivs<T>::m;

  // Time derivatives of ADM variables
  const mat3<T, DN, DN> k;
  const T dtalpha;
  const vec3<T, UP> dtbeta;

  // Spatial derivatives of ADM variables
  const mat3<vec3<T, DN>, DN, DN> dgamma;
  const mat3<mat3<T, DN, DN>, DN, DN> ddgamma;
  const vec3<T, DN> dalpha;

  // Second time derivatives of ADM variables
  const mat3<T, DN, DN> dtk;
  const T dt2alpha;
  const vec3<T, UP> dt2beta;

  // Space-time derivatives of ADM variables
  const mat3<vec3<T, DN>, DN, DN> dk;
  const vec3<T, DN> ddtalpha;
  const vec3<vec3<T, DN>, UP> ddtbeta;

  // Second spatial derivatives of ADM variables
  const mat3<T, DN, DN> ddalpha;
  const vec3<vec3<T, DN>, UP> dbeta;
  const vec3<mat3<T, DN, DN>, UP> ddbeta;

  // Intermediate quantities
  const mat3<T, DN, DN> dtgamma;
  const vec3<T, DN> dtbetal;
  const vec3<vec3<T, DN>, DN> dbetal;

  const mat3<T, DN, DN> dt2gamma;
  const mat3<vec3<T, DN>, DN, DN> ddtgamma;
  const vec3<T, DN> dt2betal;
  const vec3<vec3<T, DN>, DN> ddtbetal;
  const vec3<mat3<T, DN, DN>, DN> ddbetal;

  // Derivatives of 4-metric
  const mat4<vec4<T, DN>, DN, DN> dg;
  const mat4<mat4<T, DN, DN>, DN, DN> ddg;
  const mat4<vec4<T, DN>, UP, UP> dgu;

  // Christoffel symbol
  const vec4<mat4<T, DN, DN>, DN> Gammal;
  const vec4<mat4<T, DN, DN>, UP> Gamma;

  const vec4<mat4<vec4<T, DN>, DN, DN>, DN> dGammal;
  const vec4<mat4<vec4<T, DN>, DN, DN>, UP> dGamma;

  // Riemann, Ricci, Weyl
  // TODO: Use Rm(a,b,c,d) == Rm(c,d,a,b)
  const amat4<amat4<T, DN, DN>, DN, DN> Rm;
  const mat4<T, DN, DN> R;
  const T Rsc;
  const amat4<amat4<T, DN, DN>, DN, DN> C;

  // Ricci and Weyl scalars
  const T Lambda;
  const T Phi00, Phi11, Phi22;
  const cplx<T> Phi10, Phi20, Phi21;
  const cplx<T> Psi0, Psi1, Psi2, Psi3, Psi4;

  // Gradient of tetrad
  const vec4<vec4<T, DN>, UP> det, dephi, detheta, der;
  const vec4<vec4<T, DN>, UP> dl, dn;
  const vec4<vec4<cplx<T>, DN>, UP> dm;

  // Newman-Penrose spin coefficients
  const cplx<T> npkappa, npsigma, nprho, nptau, npepsilon, npbeta, npalpha,
      npgamma, nppi, npmu, nplambda, npnu;

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST weyl_vars(
      const T time, const vec3<T, UP> &coord3,
      //
      const mat3<T, DN, DN> &gamma, const T &alpha, const vec3<T, UP> &beta,
      //
      const mat3<T, DN, DN> &k, const T &dtalpha, const vec3<T, UP> &dtbeta,
      //
      const mat3<vec3<T, DN>, DN, DN> &dgamma, const vec3<T, DN> &dalpha,
      const vec3<vec3<T, DN>, UP> &dbeta,
      //
      const mat3<T, DN, DN> &dtk, const T &dt2alpha, const vec3<T, UP> &dt2beta,
      //
      const mat3<vec3<T, DN>, DN, DN> &dk, const vec3<T, DN> &ddtalpha,
      const vec3<vec3<T, DN>, UP> &ddtbeta,
      //
      const mat3<mat3<T, DN, DN>, DN, DN> &ddgamma,
      const mat3<T, DN, DN> &ddalpha, const vec3<mat3<T, DN, DN>, UP> &ddbeta)
      : weyl_vars_noderivs<T>(time, coord3, gamma, alpha, beta),
        // Time derivatives of ADM variables
        k(k), dtalpha(dtalpha), dtbeta(dtbeta),
        // Spatial derivatives of ADM variables
        dgamma(dgamma), ddgamma(ddgamma), dalpha(dalpha),
        // Second time derivatives of ADM variables
        dtk(dtk), dt2alpha(dt2alpha), dt2beta(dt2beta),
        // Space-time derivatives of ADM variables
        dk(dk), ddtalpha(ddtalpha), ddtbeta(ddtbeta),
        // Second spatial derivatives of ADM variables
        ddalpha(ddalpha), dbeta(dbeta), ddbeta(ddbeta),
        //
        // dt gamma_ij = -2 alpha K_ij
        //               + gamma_kj beta^k,i + gamma_ik beta^k,j
        //               + beta^k gamma_ij,k
        dtgamma([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return -2 * alpha * k(a, b) //
                 + sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                     return gamma(x, b) * dbeta(x)(a);
                   }) //
                 + sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                     return gamma(a, x) * dbeta(x)(b);
                   }) //
                 + sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                     return beta(x) * dgamma(a, b)(x);
                   });
        }),
        dtbetal([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
            return dtgamma(a, x) * beta(x) + gamma(a, x) * dtbeta(x);
          });
        }),
        dbetal([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return vec3<T, DN>([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int b) {
            return sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
              return dgamma(a, x)(b) * beta(x) + gamma(a, x) * dbeta(x)(b);
            });
          });
        }),
        //
        dt2gamma([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return -2 * dtalpha * k(a, b)  //
                 - 2 * alpha * dtk(a, b) //
                 + sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                     return dtgamma(x, b) * dbeta(x)(a) +
                            gamma(x, b) * ddtbeta(x)(a);
                   }) //
                 + sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                     return dtgamma(a, x) * dbeta(x)(b) +
                            gamma(a, x) * ddtbeta(x)(b);
                   }) //
                 + sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                     return dtbeta(x) * dgamma(a, b)(x) +
                            beta(x) * ddtgamma(a, b)(x);
                   });
        }),
        ddtgamma([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return vec3<T, DN>([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int c) {
            return -2 * dalpha(c) * k(a, b)  //
                   - 2 * alpha * dk(a, b)(c) //
                   + sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                       return dgamma(x, b)(c) * dbeta(x)(a) +
                              gamma(x, b) * ddbeta(x)(a, c);
                     }) //
                   + sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                       return dgamma(a, x)(c) * dbeta(x)(b) +
                              gamma(a, x) * ddbeta(x)(b, c);
                     }) //
                   + sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                       return dbeta(x)(c) * dgamma(a, b)(x) +
                              beta(x) * ddgamma(a, b)(x, c);
                     });
          });
        }),
        dt2betal([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
            return dt2gamma(a, x) * beta(x)        //
                   + 2 * dtgamma(a, x) * dtbeta(x) //
                   + dtgamma(a, x) * dt2beta(x);
          });
        }),
        ddtbetal([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return vec3<T, DN>([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int b) {
            return sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
              return ddtgamma(a, x)(b) * beta(x)   //
                     + dgamma(a, x)(b) * dtbeta(x) //
                     + dtgamma(a, x) * dbeta(x)(b) //
                     + gamma(a, x) * ddtbeta(x)(b);
            });
          });
        }),
        ddbetal([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return mat3<T, DN, DN>(
              [&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int b, int c) {
                return sum1([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int x) {
                  return ddgamma(a, x)(b, c) * beta(x)   //
                         + dgamma(a, x)(b) * dbeta(x)(c) //
                         + dgamma(a, x)(c) * dbeta(x)(b) //
                         + gamma(a, x) * ddbeta(x)(b, c);
                });
              });
        }),
        //
        dg([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return vec4<T, DN>([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int c) {
            if (c == 0) {
              if (a == 0 && b == 0)
                return dtalpha;
              if (a == 0)
                return dtbetal(b - 1);
              if (b == 0)
                return dtbetal(a - 1);
              return dtgamma(a - 1, b - 1);
            }
            if (a == 0 && b == 0)
              return dalpha(c - 1);
            if (a == 0)
              return dbetal(b - 1)(c - 1);
            if (b == 0)
              return dbetal(a - 1)(c - 1);
            return dgamma(a - 1, b - 1)(c - 1);
          });
        }),
        //
        ddg([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return mat4<T, DN, DN>(
              [&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int c, int d) {
                if (c == 0 && d == 0) {
                  if (a == 0 && b == 0)
                    return dt2alpha;
                  if (a == 0)
                    return dt2betal(b - 1);
                  if (b == 0)
                    return dt2betal(a - 1);
                  return dt2gamma(a - 1, b - 1);
                }
                if (c == 0) {
                  if (a == 0 && b == 0)
                    return ddtalpha(d - 1);
                  if (a == 0)
                    return ddtbetal(b - 1)(d - 1);
                  if (b == 0)
                    return ddtbetal(a - 1)(d - 1);
                  return ddtgamma(a - 1, b - 1)(d - 1);
                }
                if (d == 0) {
                  if (a == 0 && b == 0)
                    return ddtalpha(c - 1);
                  if (a == 0)
                    return ddtbetal(b - 1)(c - 1);
                  if (b == 0)
                    return ddtbetal(a - 1)(c - 1);
                  return ddtgamma(a - 1, b - 1)(c - 1);
                }
                if (a == 0 && b == 0)
                  return ddalpha(c - 1, d - 1);
                if (a == 0)
                  return ddbetal(b - 1)(c - 1, d - 1);
                if (b == 0)
                  return ddbetal(a - 1)(c - 1, d - 1);
                return ddgamma(a - 1, b - 1)(c - 1, d - 1);
              });
        }),
        //
        dgu(calc_dgu(gu, dg)),
        //
        Gammal(calc_gammal(dg)),       //
        Gamma(calc_gamma(gu, Gammal)), //
        //
        dGammal(calc_dgammal(ddg)),                    //
        dGamma(calc_dgamma(gu, dgu, Gammal, dGammal)), //
        //
        Rm(calc_riemann(g, Gamma, dGamma)), //
        R(calc_ricci(gu, Rm)),              //
        Rsc(R.trace(gu)),                   //
        C(calc_weyl(g, Rm, R, Rsc)),
        //
        // Badri Krishnan's PhD thesis, appendix A
        Lambda(Rsc / 24), //
        Phi00(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return R(a, b) * l(a) * l(b) / 2;
        })),
        Phi11(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return R(a, b) * (l(a) * n(b) + real(m(a) * conj(m(b)))) / 4;
        })),
        Phi22(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return R(a, b) * n(a) * n(b) / 2;
        })),
        Phi10(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return R(a, b) * l(a) * conj(m(b)) / T(2);
        })),
        Phi20(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return R(a, b) * conj(m(a)) * conj(m(b)) / T(2);
        })),
        Phi21(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return R(a, b) * conj(m(a)) * n(b) / T(2);
        })),
        Psi0(
            sum44([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b, int c, int d) {
              return C(a, b)(c, d) * l(a) * m(b) * l(c) * m(d);
            })),
        Psi1(
            sum44([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b, int c, int d) {
              return C(a, b)(c, d) * l(a) * m(b) * l(c) * n(d);
            })),
        Psi2(
            sum44([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b, int c, int d) {
              return C(a, b)(c, d) * l(a) * m(b) * conj(m(c)) * n(d);
            })),
        Psi3(
            sum44([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b, int c, int d) {
              return C(a, b)(c, d) * l(a) * n(b) * conj(m(c)) * n(d);
            })),
        Psi4(
            sum44([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b, int c, int d) {
              return C(a, b)(c, d) * conj(m(a)) * n(b) * conj(m(c)) * n(d);
            })),
        det(calc_det(gu, dgu, et, Gamma)),                                    //
        dephi(calc_dephi(coord, g, dg, ephi, Gamma)),                         //
        detheta(calc_detheta(coord, g, dg, ephi, dephi, etheta, Gamma)),      //
        der(calc_der(coord, g, dg, ephi, dephi, etheta, detheta, er, Gamma)), //
        dl([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return (det(a) + der(a)) / sqrt(T(2));
        }),
        dn([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return (det(a) - der(a)) / sqrt(T(2));
        }),
        dm([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a) {
          return vec4<cplx<T>, DN>([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int b) {
            return cplx<T>(detheta(a)(b), dephi(a)(b)) / sqrt(T(2));
          });
        }),
        npkappa(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return -m(a) * l(b) * dl(a)(b);
        })),
        npsigma(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return -m(a) * m(b) * dl(a)(b);
        })),
        nprho(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return -m(a) * conj(m(b)) * dl(a)(b);
        })),
        nptau(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return -m(a) * n(b) * dl(a)(b);
        })),
        npepsilon(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return (conj(m(a)) * l(b) * dm(a)(b) - n(a) * l(b) * dl(a)(b)) / T(2);
        })),
        npbeta(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return (conj(m(a)) * m(b) * dm(a)(b) - n(a) * m(b) * dl(a)(b)) / T(2);
        })),
        npalpha(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return (conj(m(a)) * conj(m(b)) * dm(a)(b) -
                  n(a) * conj(m(b)) * dl(a)(b)) /
                 T(2);
        })),
        npgamma(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return (conj(m(a)) * n(b) * dm(a)(b) - n(a) * n(b) * dl(a)(b)) / T(2);
        })),
        nppi(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return conj(m(a)) * l(b) * dn(a)(b);
        })),
        npmu(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return conj(m(a)) * m(b) * dn(a)(b);
        })),
        nplambda(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return conj(m(a)) * conj(m(b)) * dn(a)(b);
        })),
        npnu(sum42([&] CCTK_ATTRIBUTE_ALWAYS_INLINE(int a, int b) {
          return conj(m(a)) * n(b) * dn(a)(b);
        }))
  //
  {}
};

} // namespace Weyl

#endif // #ifndef WEYL_HXX
