#include "derivs.hxx"
#include "z4c_vars.hxx"

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

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_ADM(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM;
  DECLARE_CCTK_PARAMETERS;

  // The apply_upwind calls below use a deriv_order/2 + 1 wide stencil. This
  // routine runs at initial, before Z4c_RHS can perform the same check.
  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2 + 1)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2 + 1);

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_chi1(layout1, chi);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_gammat1{
      GF3D2<const CCTK_REAL>(layout1, gammatxx),
      GF3D2<const CCTK_REAL>(layout1, gammatxy),
      GF3D2<const CCTK_REAL>(layout1, gammatxz),
      GF3D2<const CCTK_REAL>(layout1, gammatyy),
      GF3D2<const CCTK_REAL>(layout1, gammatyz),
      GF3D2<const CCTK_REAL>(layout1, gammatzz)};

  const GF3D2<const CCTK_REAL> gf_Kh1(layout1, Kh);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_At1{
      GF3D2<const CCTK_REAL>(layout1, Atxx),
      GF3D2<const CCTK_REAL>(layout1, Atxy),
      GF3D2<const CCTK_REAL>(layout1, Atxz),
      GF3D2<const CCTK_REAL>(layout1, Atyy),
      GF3D2<const CCTK_REAL>(layout1, Atyz),
      GF3D2<const CCTK_REAL>(layout1, Atzz)};

  const vec<GF3D2<const CCTK_REAL>, 3> gf_Gamt1{
      GF3D2<const CCTK_REAL>(layout1, Gamtx),
      GF3D2<const CCTK_REAL>(layout1, Gamty),
      GF3D2<const CCTK_REAL>(layout1, Gamtz)};

  const GF3D2<const CCTK_REAL> gf_Theta1(layout1, Theta);

  const GF3D2<const CCTK_REAL> gf_alphaG1(layout1, alphaG);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_betaG1{
      GF3D2<const CCTK_REAL>(layout1, betaGx),
      GF3D2<const CCTK_REAL>(layout1, betaGy),
      GF3D2<const CCTK_REAL>(layout1, betaGz)};

  const GF3D2<const CCTK_REAL> gf_A1(layout1, A);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_B1{
      GF3D2<const CCTK_REAL>(layout1, Bx), GF3D2<const CCTK_REAL>(layout1, By),
      GF3D2<const CCTK_REAL>(layout1, Bz)};

  const smat<GF3D2<CCTK_REAL>, 3> gf_g1{
      GF3D2<CCTK_REAL>(layout1, gxx), GF3D2<CCTK_REAL>(layout1, gxy),
      GF3D2<CCTK_REAL>(layout1, gxz), GF3D2<CCTK_REAL>(layout1, gyy),
      GF3D2<CCTK_REAL>(layout1, gyz), GF3D2<CCTK_REAL>(layout1, gzz)};

  const smat<GF3D2<CCTK_REAL>, 3> gf_K1{
      GF3D2<CCTK_REAL>(layout1, kxx), GF3D2<CCTK_REAL>(layout1, kxy),
      GF3D2<CCTK_REAL>(layout1, kxz), GF3D2<CCTK_REAL>(layout1, kyy),
      GF3D2<CCTK_REAL>(layout1, kyz), GF3D2<CCTK_REAL>(layout1, kzz)};

  const GF3D2<CCTK_REAL> gf_alp1(layout1, alp);

  const GF3D2<CCTK_REAL> gf_dtalp1(layout1, dtalp);

  const vec<GF3D2<CCTK_REAL>, 3> gf_beta1{GF3D2<CCTK_REAL>(layout1, betax),
                                          GF3D2<CCTK_REAL>(layout1, betay),
                                          GF3D2<CCTK_REAL>(layout1, betaz)};

  const vec<GF3D2<CCTK_REAL>, 3> gf_dtbeta1{GF3D2<CCTK_REAL>(layout1, dtbetax),
                                            GF3D2<CCTK_REAL>(layout1, dtbetay),
                                            GF3D2<CCTK_REAL>(layout1, dtbetaz)};

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4c_ADM::adm");
#endif
  grid.loop_all_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);

        // Load and calculate
        const z4c_vars_noderivs<vreal> vars(
            set_Theta_zero, kappa1, kappa2, f_mu_L, f_mu_S, eta, evolveA,
            evolveB, alphaDriver, betaDriver, //
            gf_chi1(mask, index1), gf_gammat1(mask, index1),
            gf_Kh1(mask, index1), gf_At1(mask, index1), gf_Gamt1(mask, index1),
            gf_Theta1(mask, index1), gf_alphaG1(mask, index1),
            gf_betaG1(mask, index1), gf_A1(mask, index1),
            gf_B1(mask, index1), //
            Arith::nan<vreal>()(), Arith::nan<vec<vreal, 3> >()(),
            Arith::nan<smat<vreal, 3> >()());

        // Store
        gf_g1.store(mask, index1, vars.g);
        gf_K1.store(mask, index1, vars.K);
        gf_alp1.store(mask, index1, vars.alpha);
        gf_dtalp1.store(mask, index1, vars.dtalpha);
        gf_beta1.store(mask, index1, vars.beta);
        gf_dtbeta1.store(mask, index1, vars.dtbeta);
      });
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif

  // vars.dtalpha and vars.dtbeta are only the source terms of the gauge
  // conditions. The conditions themselves are advected along the shift,
  //
  //     d/dt alpha  = -alpha f_mu_L Khat + beta^i d_i alpha
  //
  // and correspondingly for the shift, so ADMBaseX::dtlapse and dtshift are
  // wrong by beta^i d_i alpha (a few percent near a black hole horizon)
  // unless that term is added here as well. It is added for both gauges: with
  // evolveA the stored value is A, and the lapse is still advected on top of
  // it.
  //
  // Kreiss-Oliger dissipation is deliberately NOT added. It is applied to the
  // evolved right hand sides in rhs.cxx because it stabilises the discrete
  // update, but it is not part of d/dt alpha in any sense that a thorn
  // reading ADMBaseX::dtlapse would want.
  //
  // A derivative needs neighbours, so unlike the loop above this covers the
  // interior only. schedule.ccl synchronises dtlapse and dtshift afterwards,
  // which carries the corrected interior values into the ghost zones and
  // applies the configured CarpetX outer boundary condition; only with
  // CarpetX::boundary_* = "none" does the outer boundary layer keep the
  // source term alone.
  apply_upwind(cctkGH, gf_alphaG1, gf_betaG1, gf_dtalp1);

  for (int a = 0; a < 3; ++a)
    apply_upwind(cctkGH, gf_betaG1(a), gf_betaG1, gf_dtbeta1(a));
}

} // namespace Z4c
