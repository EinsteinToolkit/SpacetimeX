#include "weyl.hxx"

#include <cctk.h>

#ifdef __CUDACC__
// Disable CCTK_DEBUG since the debug information takes too much
// parameter space to launch the kernels
#ifdef CCTK_DEBUG
#undef CCTK_DEBUG
#endif
#endif

#include "derivs_spacetimex.hxx"
#include "physics.hxx"
#include "weyl_vars.hxx"

#include <defs.hxx>
#include <loop_device.hxx>
#include <mat.hxx>
#include <rten.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Weyl {
using namespace Arith;
using namespace Loop;
using namespace std;

namespace {
inline vect<int, dim> box_int_imin(const GridDescBase &grid,
                                   const std::array<int, dim> &nghostzones) {
  vect<int, dim> imin, imax;
  grid.box_int<0, 0, 0>(nghostzones, imin, imax);
  return imin;
}

inline vect<int, dim> box_int_imax(const GridDescBase &grid,
                                   const std::array<int, dim> &nghostzones) {
  vect<int, dim> imin, imax;
  grid.box_int<0, 0, 0>(nghostzones, imin, imax);
  return imax;
}
} // namespace

#define GETVAR2(TYPE, NAME)                                                    \
  ([&]() {                                                                     \
    DECLARE_CCTK_ARGUMENTS_Weyl_Weyl;                                          \
    return GF3D2<TYPE>(layout1, NAME);                                         \
  }())

#define GETVAR5(TYPE, NAME)                                                    \
  ([&]() {                                                                     \
    DECLARE_CCTK_ARGUMENTS_Weyl_Weyl;                                          \
    return GF3D5<TYPE>(layout5, NAME);                                         \
  }())

gfs_t::gfs_t(const cGH *const cctkGH)
    : cctkGH(cctkGH),
      //
      indextype{0, 0, 0}, nghostzones{cctkGH->cctk_nghostzones[0],
                                      cctkGH->cctk_nghostzones[1],
                                      cctkGH->cctk_nghostzones[2]},
      grid(cctkGH), imin(box_int_imin(grid, nghostzones)),
      imax(box_int_imax(grid, nghostzones)), layout1(cctkGH, indextype),
      layout5(cctkGH, indextype), layout0(imin, imax),
      //
      gf_alpha1(GETVAR2(const CCTK_REAL, alp)),
      gf_beta1{GETVAR2(const CCTK_REAL, betax), GETVAR2(const CCTK_REAL, betay),
               GETVAR2(const CCTK_REAL, betaz)},
      gf_gamma1{GETVAR2(const CCTK_REAL, gxx), GETVAR2(const CCTK_REAL, gxy),
                GETVAR2(const CCTK_REAL, gxz), GETVAR2(const CCTK_REAL, gyy),
                GETVAR2(const CCTK_REAL, gyz), GETVAR2(const CCTK_REAL, gzz)},
      gf_K1{GETVAR2(const CCTK_REAL, kxx), GETVAR2(const CCTK_REAL, kxy),
            GETVAR2(const CCTK_REAL, kxz), GETVAR2(const CCTK_REAL, kyy),
            GETVAR2(const CCTK_REAL, kyz), GETVAR2(const CCTK_REAL, kzz)},
      gf_dtalpha1(GETVAR2(const CCTK_REAL, dtalp)),
      gf_dtbeta1{GETVAR2(const CCTK_REAL, dtbetax),
                 GETVAR2(const CCTK_REAL, dtbetay),
                 GETVAR2(const CCTK_REAL, dtbetaz)},
      gf_dtK1{GETVAR2(const CCTK_REAL, kxx), GETVAR2(const CCTK_REAL, kxy),
              GETVAR2(const CCTK_REAL, kxz), GETVAR2(const CCTK_REAL, kyy),
              GETVAR2(const CCTK_REAL, kyz), GETVAR2(const CCTK_REAL, kzz)},
      gf_dt2alpha1(GETVAR2(const CCTK_REAL, dt2alp)),
      gf_dt2beta1{GETVAR2(const CCTK_REAL, dt2betax),
                  GETVAR2(const CCTK_REAL, dt2betay),
                  GETVAR2(const CCTK_REAL, dt2betaz)},
      //
      gf_Psi0re5(GETVAR5(CCTK_REAL, Psi0re)),
      gf_Psi0im5(GETVAR5(CCTK_REAL, Psi0im)),
      gf_Psi1re5(GETVAR5(CCTK_REAL, Psi1re)),
      gf_Psi1im5(GETVAR5(CCTK_REAL, Psi1im)),
      gf_Psi2re5(GETVAR5(CCTK_REAL, Psi2re)),
      gf_Psi2im5(GETVAR5(CCTK_REAL, Psi2im)),
      gf_Psi3re5(GETVAR5(CCTK_REAL, Psi3re)),
      gf_Psi3im5(GETVAR5(CCTK_REAL, Psi3im)),
      gf_Psi4re5(GETVAR5(CCTK_REAL, Psi4re)),
      gf_Psi4im5(GETVAR5(CCTK_REAL, Psi4im)),
      //
      nvars(371), ivar(0), vars(layout0, nvars),
      //
      gf_alpha0(make_gf()), gf_dalpha0(make_vec_gf()),
      gf_ddalpha0(make_mat_gf()), gf_beta0(make_vec_gf()),
      gf_dbeta0(make_vec_vec_gf()), gf_ddbeta0(make_vec_mat_gf()),
      gf_gamma0(make_mat_gf()), gf_dgamma0(make_mat_vec_gf()),
      gf_ddgamma0(make_mat_mat_gf()), gf_K0(make_mat_gf()),
      gf_dK0(make_mat_vec_gf()),
      //
      gf_dtalpha0(make_gf()), gf_ddtalpha0(make_vec_gf()),
      gf_dtbeta0(make_vec_gf()), gf_ddtbeta0(make_vec_vec_gf()),
      gf_dtK0(make_mat_gf()),
      //
      gf_dt2alpha0(make_gf()), gf_dt2beta0(make_vec_gf()),
      // Intermediate variables: 4-metric
      tile_g4(make_mat4_gf()), tile_dg4(make_mat4_vec4_gf()),
      tile_ddg4(make_mat4_mat4_gf()),
      // Intermediate variables: 4-curvature
      tile_Gamma4(make_vec4_mat4_gf()), tile_R4(make_mat4_gf()),
      tile_C4(make_rten4_gf())
//
{
  if (ivar != nvars)
    CCTK_VERROR("Allocated nvars=%d temporary variables, but used only %d of "
                "these. Update the definition of `nvars`.",
                nvars, ivar);

  calc_derivs2(cctkGH, gf_alpha1, gf_alpha0, gf_dalpha0, gf_ddalpha0, layout0);
  calc_derivs2(cctkGH, gf_beta1, gf_beta0, gf_dbeta0, gf_ddbeta0, layout0);
  calc_derivs2(cctkGH, gf_gamma1, gf_gamma0, gf_dgamma0, gf_ddgamma0, layout0);
  calc_derivs(cctkGH, gf_K1, gf_K0, gf_dK0, layout0);

  calc_derivs(cctkGH, gf_dtalpha1, gf_dtalpha0, gf_ddtalpha0, layout0);
  calc_derivs(cctkGH, gf_dtbeta1, gf_dtbeta0, gf_ddtbeta0, layout0);
  calc_copy(cctkGH, gf_dtK1, gf_dtK0, layout0);

  calc_copy(cctkGH, gf_dt2alpha1, gf_dt2alpha0, layout0);
  calc_copy(cctkGH, gf_dt2beta1, gf_dt2beta0, layout0);
}

#undef GETVAR2
#undef GETVAR5

extern "C" void Weyl_Weyl(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Weyl_Weyl;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2 + 1)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2 + 1);

  {
    const std::array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                           cctkGH->cctk_nghostzones[1],
                                           cctkGH->cctk_nghostzones[2]};
    const GridDescBaseDevice grid(cctkGH);
    vect<int, dim> imin, imax;
    grid.box_int<0, 0, 0>(nghostzones, imin, imax);
    const auto isize = imax - imin;
    const auto np = prod(isize);
    assert(np > 0);
  }

  const gfs_t gfs(cctkGH);
  gfs.calc_metric();
  gfs.calc_curvature();
  gfs.calc_scalars();
}

} // namespace Weyl
