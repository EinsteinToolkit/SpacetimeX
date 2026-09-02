#include "physics.hxx"

#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef __CUDACC__
#include <nvtx3/nvToolsExt.h>
#endif

#include <cmath>
#include <sstream>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_Enforce(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Z4c_Enforce;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<CCTK_REAL> &gf_chi = chi;

  const GF3D2<CCTK_REAL> &gf_alphaG = alphaG;

  const smat<GF3D2<CCTK_REAL>, 3> gf_gammat{
      gammatxx, gammatxy, gammatxz, gammatyy, gammatyz, gammatzz,
  };

  const smat<GF3D2<CCTK_REAL>, 3> gf_At{
      Atxx, Atxy, Atxz, Atyy, Atyz, Atzz,
  };

  const GF3D2<CCTK_REAL> &gf_Kh = Kh;

  const vec<GF3D2<CCTK_REAL>, 3> gf_Gamt{Gamtx, Gamty, Gamtz};

  const GF3D2<CCTK_REAL> &gf_Theta = Theta;

  const GF3D2<CCTK_REAL> &gf_A = A;

  const vec<GF3D2<CCTK_REAL>, 3> gf_betaG{betaGx, betaGy, betaGz};

  const vec<GF3D2<CCTK_REAL>, 3> gf_B{Bx, By, Bz};

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const auto delta3 = one<smat<vreal, 3> >()();

#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4c_Enforce::enforce");
#endif
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);

        // Load
        const vreal chi_old = gf_chi(mask, index1);
        const vreal alphaG_old = gf_alphaG(mask, index1);

        const smat<vreal, 3> gammat_unclamped = gf_gammat(mask, index1);
        const smat<vreal, 3> At_unclamped = gf_At(mask, index1);
        const vec<vreal, 3> betaG_old = gf_betaG(mask, index1);

        const vreal Kh_old = gf_Kh(mask, index1);
        const vec<vreal, 3> Gamt_old = gf_Gamt(mask, index1);
        const vreal Theta_old = gf_Theta(mask, index1);
        const vreal A_old = gf_A(mask, index1);
        const vec<vreal, 3> B_old = gf_B(mask, index1);

        // Enforce floors and ceilings

        const vreal chi_floored = fmax(vreal(chi_floor - 1), chi_old);
        const vreal chi =
            chi_max > 0 ? fmin(vreal(chi_max - 1), chi_floored) : chi_floored;
        const vreal alphaG_floored = fmax(vreal(alphaG_floor - 1), alphaG_old);
        const vreal alphaG = alpha_max > 0
                                 ? fmin(vreal(alpha_max - 1), alphaG_floored)
                                 : alphaG_floored;

        // Clamps (see param.ccl). They act only where values are extreme.
        const auto clampto = [=](const vreal &x, const CCTK_REAL bound)
                                 ARITH_INLINE {
          return bound > 0 ? fmax(vreal(-bound), fmin(vreal(bound), x)) : x;
        };
        const auto clamp = [=](const vreal &x) ARITH_INLINE {
          return clampto(x, clamp_max);
        };

        const vreal Kh = clamp(Kh_old);
        const vec<vreal, 3> Gamt(
            [&](int a) ARITH_INLINE { return clamp(Gamt_old(a)); });
        const vreal Theta = clamp(Theta_old);
        const vreal A_ = clamp(A_old);
        const vec<vreal, 3> B_([&](int a) ARITH_INLINE { return clamp(B_old(a)); });
        const vec<vreal, 3> betaG([&](int a) ARITH_INLINE {
          return clampto(betaG_old(a), beta_max);
        });

        // Conformal metric. Where delta + gammat has stopped being positive
        // definite (Sylvester's criterion, with a small margin) it is reset to
        // flat: the det normalisation below would otherwise rescale an
        // indefinite matrix by a negative factor. Then the components are
        // clamped, det gamma_tilde = 1 is enforced, and the components are
        // clamped again, because normalising a metric with a tiny determinant
        // can undo the first clamp. Where the second clamp acts, det
        // gamma_tilde deviates from 1; the inverse used below is computed from
        // the actual determinant. None of this changes anything where no bound
        // is reached.
        const auto posdef = [=](const smat<vreal, 3> &g) ARITH_INLINE {
          const vreal m1 = g(0, 0);
          const vreal m2 = g(0, 0) * g(1, 1) - g(0, 1) * g(0, 1);
          const vreal m3 = calc_det(g);
          const CCTK_REAL eps = 1.0e-6;
          return (m1 > vreal(eps)) && (m2 > vreal(eps)) && (m3 > vreal(eps));
        };
        const auto reset_where = [=](const vbool &bad,
                                     const smat<vreal, 3> &g) ARITH_INLINE {
          return smat<vreal, 3>([&](int a, int b) ARITH_INLINE {
            return if_else(bad, vreal(0), g(a, b));
          });
        };
        const smat<vreal, 3> gammat_c1_raw([&](int a, int b) ARITH_INLINE {
          return clampto(gammat_unclamped(a, b), gammat_max);
        });
        const smat<vreal, 3> gammat_c1 =
            reset_where(!posdef(delta3 + gammat_c1_raw), gammat_c1_raw);

        // Enforce algebraic constraints
        // See arXiv:1212.2901 [gr-qc].

        const vreal detgammat_c1 = calc_det(delta3 + gammat_c1);
        const vreal chi1_c1 = 1 / cbrt(detgammat_c1) - 1;
        const smat<vreal, 3> gammat_n([&](int a, int b) ARITH_INLINE {
          return (1 + chi1_c1) * (delta3(a, b) + gammat_c1(a, b)) -
                 delta3(a, b);
        });
        const smat<vreal, 3> gammat_c2([&](int a, int b) ARITH_INLINE {
          return clampto(gammat_n(a, b), gammat_max);
        });
        const smat<vreal, 3> gammat =
            reset_where(!posdef(delta3 + gammat_c2), gammat_c2);
        const vreal detgammat = calc_det(delta3 + gammat);
#ifdef CCTK_DEBUG
        if (gammat_max == 0) {
          const vreal gammat_norm = maxabs(delta3 + gammat);
          const vreal gammat_scale = gammat_norm;
#if !defined __CUDACC__ && !defined __HIPCC__
          if (!(all(fabs(detgammat - 1) <= 1.0e-12 * gammat_scale))) {
            ostringstream buf;
            buf << "det gammat is not one: gammat=" << gammat
                << " det(gammat)=" << detgammat;
            CCTK_VERROR("%s", buf.str().c_str());
          }
#endif
          assert(all(fabs(detgammat - 1) <= 1.0e-12 * gammat_scale));
        }
#endif

        const smat<vreal, 3> gammatu =
            calc_inv(delta3 + gammat, detgammat) - delta3;

        // Trace-free part of A_tilde, from the clamped A_tilde with the inverse
        // of the final conformal metric; clamped once more afterwards, since
        // the trace removal can exceed the bound where the metric is extreme.
        // The stored A_tilde is trace-free wherever no bound is reached.
        const smat<vreal, 3> At_c1([&](int a, int b) ARITH_INLINE {
          return clamp(At_unclamped(a, b));
        });
        const vreal traceAt_c1 = sum_symm<3>([&](int x, int y) ARITH_INLINE {
          return (delta3(x, y) + gammatu(x, y)) * At_c1(x, y);
        });
        const smat<vreal, 3> At_tf([&](int a, int b) ARITH_INLINE {
          return At_c1(a, b) - traceAt_c1 / 3 * (delta3(a, b) + gammat(a, b));
        });
        const smat<vreal, 3> At([&](int a, int b) ARITH_INLINE {
          return clamp(At_tf(a, b));
        });
#ifdef CCTK_DEBUG
        if (clamp_max == 0 && gammat_max == 0) {
          const vreal traceAt = sum_symm<3>([&](int x, int y) ARITH_INLINE {
            return (delta3(x, y) + gammatu(x, y)) * At(x, y);
          });
          const vreal gammatu_norm = maxabs(delta3 + gammatu);
          const vreal At_norm = maxabs(At);
          const vreal At_scale =
              fmax(fmax(maxabs(delta3 + gammat), gammatu_norm), At_norm);
#if !defined __CUDACC__ && !defined __HIPCC__
          if (!(all(fabs(traceAt) <= 1.0e-12 * At_scale))) {
            ostringstream buf;
            buf << "tr At: At=" << At << " tr(At)=" << traceAt;
            CCTK_VERROR("%s", buf.str().c_str());
          }
#endif
          assert(all(fabs(traceAt) <= 1.0e-12 * At_scale));
        }
#endif

        // Store
        gf_chi.store(mask, index1, chi);
        gf_gammat.store(mask, index1, gammat);
        gf_Kh.store(mask, index1, Kh);
        gf_At.store(mask, index1, At);
        gf_Gamt.store(mask, index1, Gamt);
        gf_Theta.store(mask, index1, Theta);
        gf_alphaG.store(mask, index1, alphaG);
        gf_betaG.store(mask, index1, betaG);
        gf_A.store(mask, index1, A_);
        gf_B.store(mask, index1, B_);
      });
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace Z4c
