#include <defs.hxx>
#include <dual.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace KerrSchild {

template <typename T>
static CCTK_DEVICE void kerr_schild(
    // Parameters
    const T m, const T a, const T epsilon,
    // Current point
    const T t, const T x, const T y, const T z,
    // Downstairs metric
    T &gtt, T &gtx, T &gty, T &gtz, T &gxx, T &gxy, T &gxz, T &gyy, T &gyz,
    T &gzz,
    // Upstairs metric
    T &gutt, T &gutx, T &guty, T &gutz, T &guxx, T &guxy, T &guxz, T &guyy,
    T &guyz, T &guzz) {
  using Arith::pow2;
  using std::max, std::sqrt;

  // Coordinate distance to center of
  // black hole. Note it moves!
  auto rho2 = pow2(x) + pow2(y) + pow2(z);

  // Spherical auxiliary coordinate r and angle theta in BH rest frame.
  auto r2 = T(0.5) * (rho2 - pow2(a)) +
            sqrt(T(0.25) * pow2(rho2 - pow2(a)) + pow2(a * z));
  auto r = sqrt(max(epsilon, r2));

  auto costheta = z / r;

  // Coefficient H. Note this transforms as a scalar.
  auto H = m * r / (pow2(r) + pow2(a * costheta));

  // Components of l_a
  auto lt = 1;
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
    const T m, const T a, const T epsilon,
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
      m, a, epsilon,
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

extern "C" void KerrSchildX_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_KerrSchildX_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  using std::abs;
  if (abs(spin) >= mass)
    CCTK_VPARAMWARN(
        "Spin parameter %g must have absolute value less than mass %g",
        double(spin), double(mass));
}

extern "C" void KerrSchildX_InitialData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_KerrSchildX_InitialData;
  DECLARE_CCTK_PARAMETERS;

  using Arith::pow2, Arith::pown;
  using std::sqrt;

  const bool write_metric = CCTK_Equals(initial_data, "kerrschildx");
  const bool write_lapse = CCTK_Equals(initial_lapse, "kerrschildx");
  const bool write_shift = CCTK_Equals(initial_shift, "kerrschildx");
  const bool write_dtlapse = CCTK_Equals(initial_dtlapse, "kerrschildx");
  const bool write_dtshift = CCTK_Equals(initial_dtshift, "kerrschildx");

  assert(write_metric || write_lapse || write_shift || write_dtlapse ||
         write_dtshift);

  // Rename grid functions out of the way. This way we can have local variables
  // with nice, short names.
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
  const auto &alp_ = alp;
  const auto &betax_ = betax;
  const auto &betay_ = betay;
  const auto &betaz_ = betaz;
  const auto &dtalp_ = dtalp;
  const auto &dtbetax_ = dtbetax;
  const auto &dtbetay_ = dtbetay;
  const auto &dtbetaz_ = dtbetaz;

  cctk_grid.loop_all_device<0, 0, 0>(
      cctk_grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Current point
        auto t = cctk_time;
        auto x = vcoordx(p.I);
        auto y = vcoordy(p.I);
        auto z = vcoordz(p.I);

        // Downstairs and upstairs metric
        CCTK_REAL gtt, gtx, gty, gtz, gxx, gxy, gxz, gyy, gyz, gzz;
        CCTK_REAL gutt, gutx, guty, gutz, guxx, guxy, guxz, guyy, guyz, guzz;
        kerr_schild(
            // Parameters
            mass, spin, epsilon,
            // Current point
            t, x, y, z,
            // Downstairs metric
            gtt, gtx, gty, gtz, gxx, gxy, gxz, gyy, gyz, gzz,
            // Upstairs metric
            gutt, gutx, guty, gutz, guxx, guxy, guxz, guyy, guyz, guzz);

        // Time derivative of metric
        CCTK_REAL dtgtt, dtgtx, dtgty, dtgtz, dtgxx, dtgxy, dtgxz, dtgyy, dtgyz,
            dtgzz;
        CCTK_REAL dtgutt, dtgutx, dtguty, dtgutz, dtguxx, dtguxy, dtguxz,
            dtguyy, dtguyz, dtguzz;
        kerr_schild_derivs(
            // Parameters
            mass, spin, epsilon,
            // Current point
            t, x, y, z,
            // Derivative direction
            1, 0, 0, 0,
            // Downstairs metric
            dtgtt, dtgtx, dtgty, dtgtz, dtgxx, dtgxy, dtgxz, dtgyy, dtgyz,
            dtgzz,
            // Upstairs metric
            dtgutt, dtgutx, dtguty, dtgutz, dtguxx, dtguxy, dtguxz, dtguyy,
            dtguyz, dtguzz);

        // x derivative of metric
        CCTK_REAL dxgtt, dxgtx, dxgty, dxgtz, dxgxx, dxgxy, dxgxz, dxgyy, dxgyz,
            dxgzz;
        CCTK_REAL dxgutt, dxgutx, dxguty, dxgutz, dxguxx, dxguxy, dxguxz,
            dxguyy, dxguyz, dxguzz;
        kerr_schild_derivs(
            // Parameters
            mass, spin, epsilon,
            // Current point
            t, x, y, z,
            // Derivative direction
            0, 1, 0, 0,
            // Downstairs metric
            dxgtt, dxgtx, dxgty, dxgtz, dxgxx, dxgxy, dxgxz, dxgyy, dxgyz,
            dxgzz,
            // Upstairs metric
            dxgutt, dxgutx, dxguty, dxgutz, dxguxx, dxguxy, dxguxz, dxguyy,
            dxguyz, dxguzz);

        // y derivative of metric
        CCTK_REAL dygtt, dygtx, dygty, dygtz, dygxx, dygxy, dygxz, dygyy, dygyz,
            dygzz;
        CCTK_REAL dygutt, dygutx, dyguty, dygutz, dyguxx, dyguxy, dyguxz,
            dyguyy, dyguyz, dyguzz;
        kerr_schild_derivs(
            // Parameters
            mass, spin, epsilon,
            // Current point
            t, x, y, z,
            // Derivative direction
            0, 0, 1, 0,
            // Downstairs metric
            dygtt, dygtx, dygty, dygtz, dygxx, dygxy, dygxz, dygyy, dygyz,
            dygzz,
            // Upstairs metric
            dygutt, dygutx, dyguty, dygutz, dyguxx, dyguxy, dyguxz, dyguyy,
            dyguyz, dyguzz);

        // z derivative of metric
        CCTK_REAL dzgtt, dzgtx, dzgty, dzgtz, dzgxx, dzgxy, dzgxz, dzgyy, dzgyz,
            dzgzz;
        CCTK_REAL dzgutt, dzgutx, dzguty, dzgutz, dzguxx, dzguxy, dzguxz,
            dzguyy, dzguyz, dzguzz;
        kerr_schild_derivs(
            // Parameters
            mass, spin, epsilon,
            // Current point
            t, x, y, z,
            // Derivative direction
            0, 0, 0, 1,
            // Downstairs metric
            dzgtt, dzgtx, dzgty, dzgtz, dzgxx, dzgxy, dzgxz, dzgyy, dzgyz,
            dzgzz,
            // Upstairs metric
            dzgutt, dzgutx, dzguty, dzgutz, dzguxx, dzguxy, dzguxz, dzguyy,
            dzguyz, dzguzz);

        // Calculate lapse and shift from the upper metric
        auto alp = 1 / sqrt(-gutt);

        auto betax = -gutx / gutt;
        auto betay = -guty / gutt;
        auto betaz = -gutz / gutt;

        // Calculate space derivatives of shift
        auto dxbetax = (-dxgutx * gutt + gutx * dxgutt) / pow2(gutt);
        auto dxbetay = (-dxguty * gutt + guty * dxgutt) / pow2(gutt);
        auto dxbetaz = (-dxgutz * gutt + gutz * dxgutt) / pow2(gutt);

        auto dybetax = (-dygutx * gutt + gutx * dygutt) / pow2(gutt);
        auto dybetay = (-dyguty * gutt + guty * dygutt) / pow2(gutt);
        auto dybetaz = (-dygutz * gutt + gutz * dygutt) / pow2(gutt);

        auto dzbetax = (-dzgutx * gutt + gutx * dzgutt) / pow2(gutt);
        auto dzbetay = (-dzguty * gutt + guty * dzgutt) / pow2(gutt);
        auto dzbetaz = (-dzgutz * gutt + gutz * dzgutt) / pow2(gutt);

        // Calculate time derivatives of lapse and shift
        auto dtalp = 0.5 / pown(sqrt(-gutt), 3) * dtgutt;

        auto dtbetax = (-dtgutx * gutt + gutx * dtgutt) / pow2(gutt);
        auto dtbetay = (-dtguty * gutt + guty * dtgutt) / pow2(gutt);
        auto dtbetaz = (-dtgutz * gutt + gutz * dtgutt) / pow2(gutt);

        // Extrinsic curvature
        // d_t g_ij = -2 \alpha K_ij + d_i \beta_j + d_j \beta_i
        auto kxx =
            (-dtgxx + (dxgxx * betax + dygxx * betay + dzgxx * betaz +
                       2 * (dxbetax * gxx + dxbetay * gxy + dxbetaz * gxz))) /
            (2 * alp);
        auto kyy =
            (-dtgyy + (dxgyy * betax + dygyy * betay + dzgyy * betaz +
                       2 * (dybetax * gxy + dybetay * gyy + dybetaz * gyz))) /
            (2 * alp);
        auto kzz =
            (-dtgzz + (dxgzz * betax + dygzz * betay + dzgzz * betaz +
                       2 * (dzbetax * gxz + dzbetay * gyz + dzbetaz * gzz))) /
            (2 * alp);
        auto kxy = (-dtgxy + (dxgxy * betax + dygxy * betay + dzgxy * betaz +
                              dxbetax * gxy + dxbetay * gyy + dxbetaz * gyz +
                              dybetax * gxx + dybetay * gxy + dybetaz * gxz)) /
                   (2 * alp);
        auto kyz = (-dtgyz + (dxgyz * betax + dygyz * betay + dzgyz * betaz +
                              dybetax * gxz + dybetay * gyz + dybetaz * gzz +
                              dzbetax * gxy + dzbetay * gyy + dzbetaz * gyz)) /
                   (2 * alp);
        auto kxz = (-dtgxz + (dxgxz * betax + dygxz * betay + dzgxz * betaz +
                              dxbetax * gxz + dxbetay * gyz + dxbetaz * gzz +
                              dzbetax * gxx + dzbetay * gxy + dzbetaz * gxz)) /
                   (2 * alp);

        if (write_metric) {
          gxx_(p.I) = gxx;
          gxy_(p.I) = gxy;
          gxz_(p.I) = gxz;
          gyy_(p.I) = gyy;
          gyz_(p.I) = gyz;
          gzz_(p.I) = gzz;
          kxx_(p.I) = kxx;
          kxy_(p.I) = kxy;
          kxz_(p.I) = kxz;
          kyy_(p.I) = kyy;
          kyz_(p.I) = kyz;
          kzz_(p.I) = kzz;
        }
        if (write_lapse) {
          // Calculate lapse from the upper metric
          alp_(p.I) = 1 / sqrt(-gutt);
        }
        if (write_shift) {
          // Calculate shift from the upper metric
          betax_(p.I) = -gutx / gutt;
          betay_(p.I) = -guty / gutt;
          betaz_(p.I) = -gutz / gutt;
        }
        if (write_dtlapse) {
          // Metric is stationary
          dtalp_(p.I) = dtalp;
        }
        if (write_dtshift) {
          // Metric is stationary
          dtbetax_(p.I) = dtbetax;
          dtbetay_(p.I) = dtbetay;
          dtbetaz_(p.I) = dtbetaz;
        }
      });
}
} // namespace KerrSchild
