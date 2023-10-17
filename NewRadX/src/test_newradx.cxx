#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>
#include <newradx.hxx>

namespace NewRadX {

#define ijk p.I
#define im3jk p.I - 3 * p.DI[0]
#define im2jk p.I - 2 * p.DI[0]
#define im1jk p.I - 1 * p.DI[0]
#define ip3jk p.I + 3 * p.DI[0]
#define ip2jk p.I + 2 * p.DI[0]
#define ip1jk p.I + 1 * p.DI[0]

#define ijm3k p.I - 3 * p.DI[1]
#define ijm2k p.I - 2 * p.DI[1]
#define ijm1k p.I - 1 * p.DI[1]
#define ijp3k p.I + 3 * p.DI[1]
#define ijp2k p.I + 2 * p.DI[1]
#define ijp1k p.I + 1 * p.DI[1]

#define ijkm3 p.I - 3 * p.DI[2]
#define ijkm2 p.I - 2 * p.DI[2]
#define ijkm1 p.I - 1 * p.DI[2]
#define ijkp3 p.I + 3 * p.DI[2]
#define ijkp2 p.I + 2 * p.DI[2]
#define ijkp1 p.I + 1 * p.DI[2]

// Flat Laplacian operator at 4th order accuracy
template <typename T>
constexpr T Laplacian_4thorder(const Loop::GF3D2<T> gf,
                               const Loop::PointDesc &p) {
  return 1.0 / 12.0 *
         ((-gf(ip2jk) + 16 * gf(ip1jk) - 30 * gf(ijk)
           -gf(im2jk) + 16 * gf(im1jk)) /
              (p.DX[0] * p.DX[0]) +
          (-gf(ijp2k) + 16 * gf(ijp1k) - 30 * gf(ijk)
           -gf(ijm2k) + 16 * gf(ijm1k)) /
              (p.DX[1] * p.DX[1]) +
          (-gf(ijkp2) + 16 * gf(ijkp1) - 30 * gf(ijk)
           -gf(ijkm2) + 16 * gf(ijkm1)) /
              (p.DX[2] * p.DX[2]));
}

// 5th order Kreiss-Oliger dissipation
template <typename T>
constexpr T KO_diss_5thorder(const Loop::GF3D2<T> gf,
                             const Loop::PointDesc &p) {
  return 1.0 / 64.0 *
         ((gf(im3jk) - 6.0 * gf(im2jk) + 15.0 * gf(im1jk) - 20.0 * gf(ijk) +
           gf(ip3jk) - 6.0 * gf(ip2jk) + 15.0 * gf(ip1jk)) /
              p.DX[0] +
          (gf(ijm3k) - 6.0 * gf(ijm2k) + 15.0 * gf(ijm1k) - 20.0 * gf(ijk) +
           gf(ijp3k) - 6.0 * gf(ijp2k) + 15.0 * gf(ijp1k)) /
              p.DX[1] +
          (gf(ijkm3) - 6.0 * gf(ijkm2) + 15.0 * gf(ijkm1) - 20.0 * gf(ijk) +
           gf(ijkp3) - 6.0 * gf(ijkp2) + 15.0 * gf(ijkp1)) /
              p.DX[2]);
}

// Test scalar field initial data
template <typename T>
constexpr void gaussian(const T t, const T A, const T W,
                        const T x, const T y, const T z,
                        T &uu, T &vv) {
  using std::exp, std::pow, std::sqrt;

  const T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

  const auto f = [&](const T v) {
    return A * exp(-pow(v, 2) / (2 * pow(W, 2)));
  };

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    // L'Hôpital
    uu = f(-t) * (1 - pow(t/W, 2));
    vv = f(-t) * (3 - pow(t/W, 2)) * (-t/pow(W,2));
  } else {
    uu = 0.5 * (f(r - t) * (r - t) +
                f(r + t) * (r + t)) / r;
    vv = 0.5 * (f(r - t) * (pow((r - t) / W, 2) - 1) -
                f(r + t) * (pow((r + t) / W, 2) - 1)) / r;
  }
}

// Scheduled functions
extern "C" void NewRadX_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_NewRadX_Init;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            gaussian(cctk_time,
                     gaussian_a0, gaussian_w0,
                     p.x - gaussian_x0, p.y - gaussian_y0,
                     p.z - gaussian_z0, uu(p.I), vv(p.I));
      });
}

extern "C" void NewRadX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_NewRadX_RHS;
  DECLARE_CCTK_PARAMETERS;

  if (grid.nghostzones[0] != 3 || grid.nghostzones[1] != 3 ||
      grid.nghostzones[2] != 3) {
    CCTK_ERROR("Invalid number of ghost zones for 4th order FD");
  }

  grid.loop_int_device<0, 0, 0>(
    grid.nghostzones,
    [=] CCTK_DEVICE(const Loop::PointDesc &p)
        CCTK_ATTRIBUTE_ALWAYS_INLINE {
          uu_rhs(p.I) = vv(p.I);
          vv_rhs(p.I) = Laplacian_4thorder(uu, p);
        });

  if (test_use_newradx) {
    NewRadX_Apply(cctkGH, uu, uu_rhs, 0, 1, n_falloff);
    NewRadX_Apply(cctkGH, vv, vv_rhs, 0, 1, n_falloff + 1);
  }

  if (test_add_dissipation) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              uu_rhs(p.I) += eps_dissipation * KO_diss_5thorder(uu, p);
              vv_rhs(p.I) += eps_dissipation * KO_diss_5thorder(vv, p);
        });
  }
}

extern "C" void NewRadX_CompareSolution(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_NewRadX_CompareSolution;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            CCTK_REAL u_analytic, v_analytic;
            gaussian(cctk_time,
                     gaussian_a0, gaussian_w0,
                     p.x - gaussian_x0, p.y - gaussian_y0,
                     p.z - gaussian_z0, u_analytic, v_analytic);
            uu_err(p.I) = uu(p.I) - u_analytic;
            vv_err(p.I) = vv(p.I) - v_analytic;
      });
}

extern "C" void NewRadX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_NewRadX_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            using std::abs, std::max;
            CCTK_REAL maxabs_duu = 0;
            for (int d = 0; d < Loop::dim; ++d) {
              const int ni = d == 0 ? 1 : 2;
              const int nj = d == 1 ? 1 : 2;
              const int nk = d == 2 ? 1 : 2;
              for (int dk = 0; dk < nk; ++dk) {
                for (int dj = 0; dj < nj; ++dj) {
                  for (int di = 0; di < ni; ++di) {
                    const auto I = p.I + di * p.DI[0] + dj * p.DI[1] + dk * p.DI[2];
                    maxabs_duu = max(maxabs_duu, abs(uu(I + p.DI[d]) - uu(I)));
                  }
                }
              }
            }
            regrid_error(p.I) = maxabs_duu;
      });
}

} // namespace NewRadX
