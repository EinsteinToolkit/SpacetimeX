#include <defs.hxx>
#include <dual.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace KerrSchildX {

// How to regularize the Kerr-Schild auxiliary radius `r` near the ring
// singularity. The two smooth variants follow `Exact`'s
// `src/metrics/Kerr_KerrSchild.F`.
enum class regularization_t { hard_max, power_law, parabolic };

// Decode the `r_regularization` parameter. Called on the host, once per
// scheduled routine, so that the loop body sees a plain enum.
static regularization_t decode_regularization(const char *const keyword) {
  if (CCTK_Equals(keyword, "max"))
    return regularization_t::hard_max;
  if (CCTK_Equals(keyword, "power-law"))
    return regularization_t::power_law;
  if (CCTK_Equals(keyword, "parabolic"))
    return regularization_t::parabolic;
  CCTK_VERROR("Unknown r_regularization \"%s\"", keyword);
  return regularization_t::hard_max; // unreachable
}

// Regularize the Kerr-Schild auxiliary radius. `r2` is the square of the
// unregularized radius, and is non-negative up to round-off.
//
// The ring singularity sits at r = 0, inside the horizon, but it is on the
// grid because Kerr-Schild coordinates are horizon penetrating. All three
// variants keep `r` bounded away from zero; only the latter two are
// differentiable, which matters for an evolution.
template <typename T>
static CCTK_DEVICE T
regularize_r(const T r2, const CCTK_REAL epsilon, const int power,
             const regularization_t mode) {
  using std::max, std::sqrt;

  switch (mode) {

  case regularization_t::hard_max:
    // Legacy behaviour: floor r^2 directly. Continuous, but has a kink at
    // r^2 = epsilon. Note that here `epsilon` bounds r^2, not r.
    return sqrt(max(T(epsilon), r2));

  case regularization_t::power_law: {
    // r <- (r^power + epsilon^power)^(1/power), smooth everywhere.
    // `power` is restricted to a power of two (see ParamCheck) so that the
    // root is a chain of square roots; `Arith::dual` implements `sqrt`
    // exactly but only supports integer exponents in `pow`.
    const int nroots = power == 2 ? 1 : power == 4 ? 2 : 3;
    T u = max(T(0), r2);
    for (int i = 1; i < nroots; ++i)
      u = u * u; // r^2 -> r^power
    CCTK_REAL epsp = 1;
    for (int i = 0; i < power; ++i)
      epsp *= epsilon;
    u = u + T(epsp);
    for (int i = 0; i < nroots; ++i)
      u = sqrt(u);
    return u;
  }

  case regularization_t::parabolic: {
    // For r < epsilon replace r by an even polynomial in r, matched
    // smoothly at r = epsilon and positive at r = 0 -- a "stuffed" black
    // hole. Written in terms of u = r^2 so that no square root is taken
    // inside the singular region.
    const CCTK_REAL e = epsilon;
    const T u = max(T(0), r2);
    if (u >= T(e * e))
      return sqrt(u);
    switch (power) {
    case 2: {
      const CCTK_REAL c0 = e / 2, c1 = 1 / (2 * e);
      return T(c0) + u * T(c1);
    }
    case 4: {
      const CCTK_REAL c0 = 3 * e / 8, c1 = 3 / (4 * e), c2 = -1 / (8 * e * e * e);
      return T(c0) + u * (T(c1) + u * T(c2));
    }
    case 6: {
      const CCTK_REAL e3 = e * e * e, e5 = e3 * e * e;
      const CCTK_REAL c0 = 5 * e / 16, c1 = 15 / (16 * e), c2 = -5 / (16 * e3),
                      c3 = 1 / (16 * e5);
      return T(c0) + u * (T(c1) + u * (T(c2) + u * T(c3)));
    }
    default: { // power == 8
      const CCTK_REAL e3 = e * e * e, e5 = e3 * e * e, e7 = e5 * e * e;
      const CCTK_REAL c0 = 35 * e / 128, c1 = 35 / (32 * e), c2 = -35 / (64 * e3),
                      c3 = 7 / (32 * e5), c4 = -5 / (128 * e7);
      return T(c0) + u * (T(c1) + u * (T(c2) + u * (T(c3) + u * T(c4))));
    }
    }
  }
  }

  // Unreachable
  return sqrt(max(T(epsilon), r2));
}

template <typename T>
static CCTK_DEVICE void kerr_schild(
    // Parameters
    const T m, const T a, const CCTK_REAL epsilon, const int power,
    const regularization_t mode,
    // Current point
    const T t, const T x, const T y, const T z,
    // Downstairs metric
    T &gtt, T &gtx, T &gty, T &gtz, T &gxx, T &gxy, T &gxz, T &gyy, T &gyz,
    T &gzz,
    // Upstairs metric
    T &gutt, T &gutx, T &guty, T &gutz, T &guxx, T &guxy, T &guxz, T &guyy,
    T &guyz, T &guzz) {
  using Arith::pow2;
  using std::sqrt;

  // Coordinate distance to the centre of the black hole, which is at the
  // origin.
  auto rho2 = pow2(x) + pow2(y) + pow2(z);

  // Spherical auxiliary coordinate r and angle theta in the BH rest frame.
  auto r2 = T(0.5) * (rho2 - pow2(a)) +
            sqrt(T(0.25) * pow2(rho2 - pow2(a)) + pow2(a * z));
  auto r = regularize_r(r2, epsilon, power, mode);

  auto costheta = z / r;

  // Coefficient H. Note this transforms as a scalar.
  auto H = m * r / (pow2(r) + pow2(a * costheta));

  // Components of l_a
  auto lt = T(1);
  auto lx = (r * x + a * y) / (pow2(r) + pow2(a));
  auto ly = (r * y - a * x) / (pow2(r) + pow2(a));
  auto lz = z / r;

  // Downstairs metric. g_ab = flat_ab + H l_a l_b
  gtt = -1 + 2 * H * lt * lt;
  gtx = 2 * H * lt * lx;
  gty = 2 * H * lt * ly;
  gtz = 2 * H * lt * lz;
  gxx = 1 + 2 * H * lx * lx;
  gyy = 1 + 2 * H * ly * ly;
  gzz = 1 + 2 * H * lz * lz;
  gxy = 2 * H * lx * ly;
  gyz = 2 * H * ly * lz;
  gxz = 2 * H * lx * lz;

  // Upstairs metric. g^ab = flat^ab - H l^a l^b
  // Notice that g^ab = g_ab and l^i = l_i and l^0 = - l_0 in flat
  // spacetime.
  gutt = -1 - 2 * H * lt * lt;
  gutx = 2 * H * lt * lx;
  guty = 2 * H * lt * ly;
  gutz = 2 * H * lt * lz;
  guxx = 1 - 2 * H * lx * lx;
  guyy = 1 - 2 * H * ly * ly;
  guzz = 1 - 2 * H * lz * lz;
  guxy = -2 * H * lx * ly;
  guyz = -2 * H * ly * lz;
  guxz = -2 * H * lx * lz;
}

template <typename T>
static CCTK_DEVICE void kerr_schild_derivs(
    // Parameters
    const T m, const T a, const CCTK_REAL epsilon, const int power,
    const regularization_t mode,
    // Current point
    const T t, const T x, const T y, const T z,
    // Derivative direction, set one of these to one, the others to zero
    const int dt, const int dx, const int dy, const int dz,
    // Downstairs metric
    T &gtt, T &gtx, T &gty, T &gtz, T &gxx, T &gxy, T &gxz, T &gyy, T &gyz,
    T &gzz,
    // Upstairs metric
    T &gutt, T &gutx, T &guty, T &gutz, T &guxx, T &guxy, T &guxz, T &guyy,
    T &guyz, T &guzz) {
  using Arith::dual;

  // Dual numbers for metric derivatives. A dual number holds the function
  // value (`.val`) and its derivative (`.eps`). similar to a complex
  // number.
  using DUAL_REAL = dual<CCTK_REAL>;
  DUAL_REAL dgtt, dgtx, dgty, dgtz, dgxx, dgxy, dgxz, dgyy, dgyz, dgzz;
  DUAL_REAL dgutt, dgutx, dguty, dgutz, dguxx, dguxy, dguxz, dguyy, dguyz,
      dguzz;

  // Derivative of metric
  kerr_schild<DUAL_REAL>(
      // Parameters
      m, a, epsilon, power, mode,
      // Current point
      DUAL_REAL(t, dt), DUAL_REAL(x, dx), DUAL_REAL(y, dy), DUAL_REAL(z, dz),
      // Downstairs metric
      dgtt, dgtx, dgty, dgtz, dgxx, dgxy, dgxz, dgyy, dgyz, dgzz,
      // Upstairs metric
      dgutt, dgutx, dguty, dgutz, dguxx, dguxy, dguxz, dguyy, dguyz, dguzz);
  gtt = dgtt.eps;
  gtx = dgtx.eps;
  gty = dgty.eps;
  gtz = dgtz.eps;
  gxx = dgxx.eps;
  gxy = dgxy.eps;
  gxz = dgxz.eps;
  gyy = dgyy.eps;
  gyz = dgyz.eps;
  gzz = dgzz.eps;
  gutt = dgutt.eps;
  gutx = dgutx.eps;
  guty = dguty.eps;
  gutz = dgutz.eps;
  guxx = dguxx.eps;
  guxy = dguxy.eps;
  guxz = dguxz.eps;
  guyy = dguyy.eps;
  guyz = dguyz.eps;
  guzz = dguzz.eps;
}

// The 3+1 ADM quantities at one point.
template <typename T> struct adm_vars_t {
  T gxx, gxy, gxz, gyy, gyz, gzz;
  T kxx, kxy, kxz, kyy, kyz, kzz;
  T alp;
  T betax, betay, betaz;
  T dtalp;
  T dtbetax, dtbetay, dtbetaz;
};

// Evaluate the Kerr-Schild ADM data at one point.
//
// `with_curv` and `with_dtgauge` select how much work is done: the lapse and
// shift need the metric only, the gauge time derivatives additionally need
// the time derivative of the metric, and the extrinsic curvature needs all
// four derivative directions. Fields that are not requested are left
// untouched.
template <typename T>
static CCTK_DEVICE adm_vars_t<T>
kerr_schild_adm(const T m, const T a, const CCTK_REAL epsilon, const int power,
                const regularization_t mode, const T t, const T x, const T y,
                const T z, const bool with_curv, const bool with_dtgauge) {
  using Arith::pow2, Arith::pown;
  using std::sqrt;

  adm_vars_t<T> v{};

  // Downstairs and upstairs metric
  T gtt, gtx, gty, gtz, gxx, gxy, gxz, gyy, gyz, gzz;
  T gutt, gutx, guty, gutz, guxx, guxy, guxz, guyy, guyz, guzz;
  kerr_schild(m, a, epsilon, power, mode, t, x, y, z, //
              gtt, gtx, gty, gtz, gxx, gxy, gxz, gyy, gyz, gzz,
              gutt, gutx, guty, gutz, guxx, guxy, guxz, guyy, guyz, guzz);

  v.gxx = gxx;
  v.gxy = gxy;
  v.gxz = gxz;
  v.gyy = gyy;
  v.gyz = gyz;
  v.gzz = gzz;

  // Calculate lapse and shift from the upper metric
  const T alp = 1 / sqrt(-gutt);
  const T betax = -gutx / gutt;
  const T betay = -guty / gutt;
  const T betaz = -gutz / gutt;

  v.alp = alp;
  v.betax = betax;
  v.betay = betay;
  v.betaz = betaz;

  if (!(with_curv || with_dtgauge))
    return v;

  // Time derivative of metric
  T dtgtt, dtgtx, dtgty, dtgtz, dtgxx, dtgxy, dtgxz, dtgyy, dtgyz, dtgzz;
  T dtgutt, dtgutx, dtguty, dtgutz, dtguxx, dtguxy, dtguxz, dtguyy, dtguyz,
      dtguzz;
  kerr_schild_derivs(m, a, epsilon, power, mode, t, x, y, z, 1, 0, 0, 0, //
                     dtgtt, dtgtx, dtgty, dtgtz, dtgxx, dtgxy, dtgxz, dtgyy,
                     dtgyz, dtgzz, dtgutt, dtgutx, dtguty, dtgutz, dtguxx,
                     dtguxy, dtguxz, dtguyy, dtguyz, dtguzz);

  // Calculate time derivatives of lapse and shift. The metric is
  // stationary, so these vanish analytically; they are computed rather than
  // hard-wired to zero so that a future boosted or moving black hole works
  // without further changes.
  v.dtalp = 0.5 / pown(sqrt(-gutt), 3) * dtgutt;
  v.dtbetax = (-dtgutx * gutt + gutx * dtgutt) / pow2(gutt);
  v.dtbetay = (-dtguty * gutt + guty * dtgutt) / pow2(gutt);
  v.dtbetaz = (-dtgutz * gutt + gutz * dtgutt) / pow2(gutt);

  if (!with_curv)
    return v;

  // x derivative of metric
  T dxgtt, dxgtx, dxgty, dxgtz, dxgxx, dxgxy, dxgxz, dxgyy, dxgyz, dxgzz;
  T dxgutt, dxgutx, dxguty, dxgutz, dxguxx, dxguxy, dxguxz, dxguyy, dxguyz,
      dxguzz;
  kerr_schild_derivs(m, a, epsilon, power, mode, t, x, y, z, 0, 1, 0, 0, //
                     dxgtt, dxgtx, dxgty, dxgtz, dxgxx, dxgxy, dxgxz, dxgyy,
                     dxgyz, dxgzz, dxgutt, dxgutx, dxguty, dxgutz, dxguxx,
                     dxguxy, dxguxz, dxguyy, dxguyz, dxguzz);

  // y derivative of metric
  T dygtt, dygtx, dygty, dygtz, dygxx, dygxy, dygxz, dygyy, dygyz, dygzz;
  T dygutt, dygutx, dyguty, dygutz, dyguxx, dyguxy, dyguxz, dyguyy, dyguyz,
      dyguzz;
  kerr_schild_derivs(m, a, epsilon, power, mode, t, x, y, z, 0, 0, 1, 0, //
                     dygtt, dygtx, dygty, dygtz, dygxx, dygxy, dygxz, dygyy,
                     dygyz, dygzz, dygutt, dygutx, dyguty, dygutz, dyguxx,
                     dyguxy, dyguxz, dyguyy, dyguyz, dyguzz);

  // z derivative of metric
  T dzgtt, dzgtx, dzgty, dzgtz, dzgxx, dzgxy, dzgxz, dzgyy, dzgyz, dzgzz;
  T dzgutt, dzgutx, dzguty, dzgutz, dzguxx, dzguxy, dzguxz, dzguyy, dzguyz,
      dzguzz;
  kerr_schild_derivs(m, a, epsilon, power, mode, t, x, y, z, 0, 0, 0, 1, //
                     dzgtt, dzgtx, dzgty, dzgtz, dzgxx, dzgxy, dzgxz, dzgyy,
                     dzgyz, dzgzz, dzgutt, dzgutx, dzguty, dzgutz, dzguxx,
                     dzguxy, dzguxz, dzguyy, dzguyz, dzguzz);

  // Calculate space derivatives of shift
  const T dxbetax = (-dxgutx * gutt + gutx * dxgutt) / pow2(gutt);
  const T dxbetay = (-dxguty * gutt + guty * dxgutt) / pow2(gutt);
  const T dxbetaz = (-dxgutz * gutt + gutz * dxgutt) / pow2(gutt);

  const T dybetax = (-dygutx * gutt + gutx * dygutt) / pow2(gutt);
  const T dybetay = (-dyguty * gutt + guty * dygutt) / pow2(gutt);
  const T dybetaz = (-dygutz * gutt + gutz * dygutt) / pow2(gutt);

  const T dzbetax = (-dzgutx * gutt + gutx * dzgutt) / pow2(gutt);
  const T dzbetay = (-dzguty * gutt + guty * dzgutt) / pow2(gutt);
  const T dzbetaz = (-dzgutz * gutt + gutz * dzgutt) / pow2(gutt);

  // Extrinsic curvature
  // d_t g_ij = -2 \alpha K_ij + \beta^k d_k g_ij
  //            + g_kj d_i \beta^k + g_ik d_j \beta^k
  v.kxx = (-dtgxx + (dxgxx * betax + dygxx * betay + dzgxx * betaz +
                     2 * (dxbetax * gxx + dxbetay * gxy + dxbetaz * gxz))) /
          (2 * alp);
  v.kyy = (-dtgyy + (dxgyy * betax + dygyy * betay + dzgyy * betaz +
                     2 * (dybetax * gxy + dybetay * gyy + dybetaz * gyz))) /
          (2 * alp);
  v.kzz = (-dtgzz + (dxgzz * betax + dygzz * betay + dzgzz * betaz +
                     2 * (dzbetax * gxz + dzbetay * gyz + dzbetaz * gzz))) /
          (2 * alp);
  v.kxy = (-dtgxy + (dxgxy * betax + dygxy * betay + dzgxy * betaz +
                     dxbetax * gxy + dxbetay * gyy + dxbetaz * gyz +
                     dybetax * gxx + dybetay * gxy + dybetaz * gxz)) /
          (2 * alp);
  v.kyz = (-dtgyz + (dxgyz * betax + dygyz * betay + dzgyz * betaz +
                     dybetax * gxz + dybetay * gyz + dybetaz * gzz +
                     dzbetax * gxy + dzbetay * gyy + dzbetaz * gyz)) /
          (2 * alp);
  v.kxz = (-dtgxz + (dxgxz * betax + dygxz * betay + dzgxz * betaz +
                     dxbetax * gxz + dxbetay * gyz + dxbetaz * gzz +
                     dzbetax * gxx + dzbetay * gxy + dzbetaz * gxz)) /
          (2 * alp);

  return v;
}

extern "C" void KerrSchildX_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_KerrSchildX_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  using std::abs;
  if (abs(spin) >= mass)
    CCTK_VPARAMWARN(
        "Spin parameter %g must have absolute value less than mass %g",
        double(spin), double(mass));

  if (CCTK_Equals(r_regularization, "power-law"))
    if (!(power == 2 || power == 4 || power == 8))
      CCTK_VPARAMWARN("r_regularization=\"power-law\" requires power to be 2, "
                      "4, or 8, but power=%d",
                      int(power));

  if (CCTK_Equals(r_regularization, "parabolic"))
    if (!(power == 2 || power == 4 || power == 6 || power == 8))
      CCTK_VPARAMWARN("r_regularization=\"parabolic\" requires power to be 2, "
                      "4, 6, or 8, but power=%d",
                      int(power));
}

extern "C" void KerrSchildX_InitialData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_KerrSchildX_InitialData;
  DECLARE_CCTK_PARAMETERS;

  const regularization_t mode = decode_regularization(r_regularization);

  // Rename grid functions out of the way. This way we can have local
  // variables with nice, short names.
  const auto &gxx_ = gxx;
  const auto &gxy_ = gxy;
  const auto &gxz_ = gxz;
  const auto &gyy_ = gyy;
  const auto &gyz_ = gyz;
  const auto &gzz_ = gzz;
  const auto &kxx_ = kxx;
  const auto &kxy_ = kxy;
  const auto &kxz_ = kxz;
  const auto &kyy_ = kyy;
  const auto &kyz_ = kyz;
  const auto &kzz_ = kzz;

  cctk_grid.loop_all_device<0, 0, 0>(
      cctk_grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto v = kerr_schild_adm<CCTK_REAL>(
            mass, spin, epsilon, power, mode, cctk_time, vcoordx(p.I),
            vcoordy(p.I), vcoordz(p.I), true, false);
        gxx_(p.I) = v.gxx;
        gxy_(p.I) = v.gxy;
        gxz_(p.I) = v.gxz;
        gyy_(p.I) = v.gyy;
        gyz_(p.I) = v.gyz;
        gzz_(p.I) = v.gzz;
        kxx_(p.I) = v.kxx;
        kxy_(p.I) = v.kxy;
        kxz_(p.I) = v.kxz;
        kyy_(p.I) = v.kyy;
        kyz_(p.I) = v.kyz;
        kzz_(p.I) = v.kzz;
      });
}

extern "C" void KerrSchildX_InitialLapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_KerrSchildX_InitialLapse;
  DECLARE_CCTK_PARAMETERS;

  const regularization_t mode = decode_regularization(r_regularization);

  const auto &alp_ = alp;

  cctk_grid.loop_all_device<0, 0, 0>(
      cctk_grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto v = kerr_schild_adm<CCTK_REAL>(
            mass, spin, epsilon, power, mode, cctk_time, vcoordx(p.I),
            vcoordy(p.I), vcoordz(p.I), false, false);
        alp_(p.I) = v.alp;
      });
}

extern "C" void KerrSchildX_InitialShift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_KerrSchildX_InitialShift;
  DECLARE_CCTK_PARAMETERS;

  const regularization_t mode = decode_regularization(r_regularization);

  const auto &betax_ = betax;
  const auto &betay_ = betay;
  const auto &betaz_ = betaz;

  cctk_grid.loop_all_device<0, 0, 0>(
      cctk_grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto v = kerr_schild_adm<CCTK_REAL>(
            mass, spin, epsilon, power, mode, cctk_time, vcoordx(p.I),
            vcoordy(p.I), vcoordz(p.I), false, false);
        betax_(p.I) = v.betax;
        betay_(p.I) = v.betay;
        betaz_(p.I) = v.betaz;
      });
}

extern "C" void KerrSchildX_InitialDtLapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_KerrSchildX_InitialDtLapse;
  DECLARE_CCTK_PARAMETERS;

  const regularization_t mode = decode_regularization(r_regularization);

  const auto &dtalp_ = dtalp;

  cctk_grid.loop_all_device<0, 0, 0>(
      cctk_grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto v = kerr_schild_adm<CCTK_REAL>(
            mass, spin, epsilon, power, mode, cctk_time, vcoordx(p.I),
            vcoordy(p.I), vcoordz(p.I), false, true);
        dtalp_(p.I) = v.dtalp;
      });
}

extern "C" void KerrSchildX_InitialDtShift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_KerrSchildX_InitialDtShift;
  DECLARE_CCTK_PARAMETERS;

  const regularization_t mode = decode_regularization(r_regularization);

  const auto &dtbetax_ = dtbetax;
  const auto &dtbetay_ = dtbetay;
  const auto &dtbetaz_ = dtbetaz;

  cctk_grid.loop_all_device<0, 0, 0>(
      cctk_grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto v = kerr_schild_adm<CCTK_REAL>(
            mass, spin, epsilon, power, mode, cctk_time, vcoordx(p.I),
            vcoordy(p.I), vcoordz(p.I), false, true);
        dtbetax_(p.I) = v.dtbetax;
        dtbetay_(p.I) = v.dtbetay;
        dtbetaz_(p.I) = v.dtbetaz;
      });
}

} // namespace KerrSchildX
