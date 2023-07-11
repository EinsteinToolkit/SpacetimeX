#include "discretization.hxx"
#include "sYlm.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <limits>
#include <random>

namespace AHFinder {

template <typename T, typename U>
constexpr bool isapprox(const T x, const U y) {
  using std::abs, std::pow;
  using R = decltype(abs(declval<T>() - declval<U>()));
  const R eps = pow(std::numeric_limits<R>::epsilon(), R(3) / R(4));
  return abs(x - y) <= eps;
}

template <typename T, typename U>
constexpr bool isapproxv(const T x, const U y) {
  using std::abs, std::pow;
  using R = decltype(abs(declval<typename T::value_type>() -
                         declval<typename U::value_type>()));
  const R eps = pow(std::numeric_limits<R>::epsilon(), R(3) / R(4));
  return maxabs(x - y) <= eps;
}

template <typename T> void test_aij_alm() {
  CCTK_VINFO("test_aij_alm...");

  const int nmodes = 5;
  const geom_t geom(nmodes);

  // Loop over all modes
  int count = 0;
  for (int ll = 0; ll < nmodes; ++ll) {
    for (int mm = 0; mm <= ll; ++mm) {
      for (int cc = 0; cc <= min(mm, 1); ++cc) { // re / im
        ++count;

        // Define coefficients (for a real function)
        alm_t<std::complex<T> > alm(geom, 0);
        for (int l = 0; l <= geom.lmax; ++l)
          for (int m = 0; m <= l; ++m)
            alm(l, m) = 0;
        if (mm == 0) {
          alm(ll, mm) = 1;
        } else {
          alm(ll, mm) = cc == 0 ? 1.0 : 1.0i;
          alm(ll, -mm) = T(bitsign(mm)) * conj(alm(ll, mm));
        }

        // Evaluate
        const aij_t<T> aij = evaluate(alm);
        // const aij_t<std::complex<T> > aij = evaluate(alm, 0);

        // Check
        if (ll <= sYlm_lmax) {
          for (int i = 0; i < geom.ntheta; ++i) {
            for (int j = 0; j < geom.nphi; ++j) {
              const T theta = geom.coord_theta(i, j);
              const T phi = geom.coord_phi(i, j);
              const T a = aij(i, j);
              std::complex<T> sc;
              if (mm == 0) {
                sc = sYlm(0, ll, mm, theta, phi);
              } else {
                const std::complex<T> c = cc == 0 ? 1.0 : 1.0i;
                sc = c * sYlm(0, ll, mm, theta, phi) +
                     T(bitsign(mm)) * conj(c) * sYlm(0, ll, -mm, theta, phi);
              }
              assert(isapprox(imag(sc), T(0)));
              const T s = real(sc);
              assert(isapprox(a, s));
            }
          }
        }

        // Expand
        const alm_t<std::complex<T> > alm1 = expand(aij, 0);
        // const alm_t<std::complex<T>> alm1 = expand(aij, 0);

        // Check
        for (int l = 0; l <= geom.lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            if (!(isapprox(alm1(l, m), alm(l, m))))
              CCTK_VINFO("ll=%d,mm=%d,cc=%d l=%d,m=%d "
                         "alm=(%.17g,%.17g) alm1=(%.17g,%.17g)",
                         ll, mm, cc, l, m, real(alm(l, m)), imag(alm(l, m)),
                         real(alm1(l, m)), imag(alm1(l, m)));
            assert(isapprox(alm1(l, m), alm(l, m)));
          }
        }

#if 0

        // Gradient
        const alm_t<std::complex<T> > dalm = grad(alm);

        // Evaluate
        const aij_t<std::complex<T> > daij = evaluate_grad(dalm);

        // Check
        if (ll <= sYlm_lmax) {
          for (int i = 0; i < geom.ntheta; ++i) {
            for (int j = 0; j < geom.nphi; ++j) {
              const T theta = geom.coord_theta(i, j);
              const T phi = geom.coord_phi(i, j);
              const std::complex<T> da = daij(i, j);
              array<std::complex<T>, 2> dsc;
              if (mm == 0) {
                dsc = dsYlm(0, ll, mm, theta, phi);
              } else {
                const std::complex<T> c = cc == 0 ? 1.0 : 1.0i;
                const array<std::complex<T>, 2> dsc_p =
                    dsYlm(0, ll, mm, theta, phi);
                const array<std::complex<T>, 2> dsc_m =
                    dsYlm(0, ll, -mm, theta, phi);
                for (int n = 0; n < 2; ++n) {
                  dsc[n] = c * dsc_p[n] + T(bitsign(mm)) * conj(c) * dsc_m[n];
                  assert(isapprox(imag(dsc[n]), T(0)));
                }
              }
              std::complex<T> ds{real(dsc[0]), real(dsc[1])};
              if (!(isapprox(da, ds)))
                CCTK_VINFO("ll=%d,mm=%d,cc=%d i=%d,j=%d,theta=%f,phi=%f "
                           "da=(%f,%f) ds=(%f,%f)",
                           ll, mm, cc, i, j, theta, phi, real(da), imag(da),
                           real(ds), imag(ds));
              assert(isapprox(da, ds));
            }
          }
        }

        // Expand
        const alm_t<std::complex<T> > dalm1 = expand_grad(daij);

        // Check
        for (int l = 0; l <= geom.lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            if (!(isapprox(dalm1(l, m), dalm(l, m))))
              CCTK_VINFO("ll=%d,mm=%d,cc=%d l=%d,m=%d "
                         "dalm=(%.17g,%.17g) dalm1=(%.17g,%.17g)",
                         ll, mm, cc, l, m, real(dalm(l, m)), imag(dalm(l, m)),
                         real(dalm1(l, m)), imag(dalm1(l, m)));
            assert(isapprox(dalm1(l, m), dalm(l, m)));
          }
        }

        // Divergence
        const alm_t<std::complex<T> > lalm = div(dalm);

        // Check
        for (int l = 0; l <= geom.lmax; ++l)
          for (int m = -l; m <= l; ++m)
            assert(isapprox(lalm(l, m), T(-l * (l + 1)) * alm(l, m)));

#endif
      }
    }
  }
  assert(count == nmodes * nmodes);
}

template <typename T> void test_scalar_aij_alm() {
  CCTK_VINFO("test_scalar_aij_alm...");

  const int nmodes = 5;
  const geom_t geom(nmodes);

  std::random_device rdev;
  std::default_random_engine reng(rdev());
  std::uniform_real_distribution<T> rdist(-1, 1);
  const auto rand = [&]() { return rdist(reng); };
  const auto randc = [&]() {
    return std::complex<T>(rdist(reng), rdist(reng));
  };

  for (int iter = 0; iter < 100; ++iter) {
    scalar_aij_t<T> n(geom), x(geom), y(geom), z(geom);
    n = 0;
    fmap_(x, [&](auto &a) { a = rand(); });
    fmap_(y, [&](auto &a) { a = rand(); });
    fmap_(z, [&](auto &a) { a = rand(); });
    const T a = rand();
    const T b = rand();

    assert(isapproxv(n, n));
    assert(!isapproxv(n, x));
    assert(isapproxv(n + x, x));
    assert(isapproxv(0 * x, n));
    assert(isapproxv(x + y, y + x));
    assert(isapproxv((x + y) + z, x + (y + z)));
    assert(isapproxv(+x, x));
    assert(isapproxv(-x, -1 * x));
    assert(isapproxv(x - y, x + (-y)));
    assert(isapproxv(x - x, n));
    assert(isapproxv(x * a, a * x));
    assert(isapproxv(a * (x + y), a * x + a * y));
    assert(isapproxv((a + b) * x, a * x + b * x));

    scalar_alm_t<std::complex<T> > xlm = expand(x);
    scalar_aij_t<T> x2 = evaluate(xlm);
    scalar_alm_t<std::complex<T> > xlm2 = expand(x2);
    scalar_aij_t<T> x3 = evaluate(xlm2);
    assert(isapproxv(x3, x2));
  }

  for (int iter = 0; iter < 100; ++iter) {
    scalar_alm_t<std::complex<T> > n(geom), x(geom), y(geom), z(geom);
    n = 0;
    fmap_(x, [&](auto &a) { a = randc(); });
    fmap_(y, [&](auto &a) { a = randc(); });
    fmap_(z, [&](auto &a) { a = randc(); });
    const std::complex<T> a = randc();
    const std::complex<T> b = randc();

    assert(isapproxv(n, n));
    assert(!isapproxv(n, x));
    assert(isapproxv(n + x, x));
    assert(isapproxv(0 * x, n));
    assert(isapproxv(x + y, y + x));
    assert(isapproxv((x + y) + z, x + (y + z)));
    assert(isapproxv(+x, x));
    assert(isapproxv(-x, -1 * x));
    assert(isapproxv(x - y, x + (-y)));
    assert(isapproxv(x - x, n));
    assert(isapproxv(x * a, a * x));
    assert(isapproxv(a * (x + y), a * x + a * y));
    assert(isapproxv((a + b) * x, a * x + b * x));

    scalar_aij_t<T> xij = evaluate(x);
    scalar_alm_t<std::complex<T> > x2 = expand(xij);
    scalar_aij_t<T> xij2 = evaluate(x2);
    scalar_alm_t<std::complex<T> > x3 = expand(xij2);
    assert(isapproxv(x3, x2));
  }
}

template <typename T> void test_vector_aij_alm() {
  CCTK_VINFO("test_vector_aij_alm...");

  const int nmodes = 5;
  const geom_t geom(nmodes);

  std::random_device rdev;
  std::default_random_engine reng(rdev());
  std::uniform_real_distribution<T> rdist(-1, 1);
  const auto rand = [&]() { return rdist(reng); };
  const auto randc = [&]() {
    return std::complex<T>(rdist(reng), rdist(reng));
  };

  for (int iter = 0; iter < 100; ++iter) {
    vector_aij_t<T> n(geom), x(geom), y(geom), z(geom);
    n = 0;
    fmap_(x, [&](auto &a) { a = rand(); });
    fmap_(y, [&](auto &a) { a = rand(); });
    fmap_(z, [&](auto &a) { a = rand(); });
    const T a = rand();
    const T b = rand();

    assert(isapproxv(n, n));
    assert(!isapproxv(n, x));
    assert(isapproxv(n + x, x));
    assert(isapproxv(0 * x, n));
    assert(isapproxv(x + y, y + x));
    assert(isapproxv((x + y) + z, x + (y + z)));
    assert(isapproxv(+x, x));
    assert(isapproxv(-x, -1 * x));
    assert(isapproxv(x - y, x + (-y)));
    assert(isapproxv(x - x, n));
    assert(isapproxv(x * a, a * x));
    assert(isapproxv(a * (x + y), a * x + a * y));
    assert(isapproxv((a + b) * x, a * x + b * x));

    vector_alm_t<std::complex<T> > xlm = expand(x);
    vector_aij_t<T> x2 = evaluate(xlm);
    vector_alm_t<std::complex<T> > xlm2 = expand(x2);
    vector_aij_t<T> x3 = evaluate(xlm2);
    assert(isapproxv(x3, x2));
  }

  for (int iter = 0; iter < 100; ++iter) {
    vector_alm_t<std::complex<T> > n(geom), x(geom), y(geom), z(geom);
    n = 0;
    fmap_(x, [&](auto &a) { a = randc(); });
    fmap_(y, [&](auto &a) { a = randc(); });
    fmap_(z, [&](auto &a) { a = randc(); });
    const std::complex<T> a = randc();
    const std::complex<T> b = randc();

    assert(isapproxv(n, n));
    assert(!isapproxv(n, x));
    assert(isapproxv(n + x, x));
    assert(isapproxv(0 * x, n));
    assert(isapproxv(x + y, y + x));
    assert(isapproxv((x + y) + z, x + (y + z)));
    assert(isapproxv(+x, x));
    assert(isapproxv(-x, -1 * x));
    assert(isapproxv(x - y, x + (-y)));
    assert(isapproxv(x - x, n));
    assert(isapproxv(x * a, a * x));
    assert(isapproxv(a * (x + y), a * x + a * y));
    assert(isapproxv((a + b) * x, a * x + b * x));

    vector_aij_t<T> xij = evaluate(x);
    vector_alm_t<std::complex<T> > x2 = expand(xij);
    vector_aij_t<T> xij2 = evaluate(x2);
    vector_alm_t<std::complex<T> > x3 = expand(xij2);
    assert(isapproxv(x3, x2));
  }
}

template <typename T> void test_tensor_aij_alm() {
  CCTK_VINFO("test_tensor_aij_alm...");

  const int nmodes = 5;
  const geom_t geom(nmodes);

  std::random_device rdev;
  std::default_random_engine reng(rdev());
  std::uniform_real_distribution<T> rdist(-1, 1);
  const auto rand = [&]() { return rdist(reng); };
  const auto randc = [&]() {
    return std::complex<T>(rdist(reng), rdist(reng));
  };

  for (int iter = 0; iter < 100; ++iter) {
    tensor_aij_t<T> n(geom), x(geom), y(geom), z(geom);
    n = 0;
    fmap_(x, [&](auto &a) { a = rand(); });
    fmap_(y, [&](auto &a) { a = rand(); });
    fmap_(z, [&](auto &a) { a = rand(); });
    const T a = rand();
    const T b = rand();

    assert(isapproxv(n, n));
    assert(!isapproxv(n, x));
    assert(isapproxv(n + x, x));
    assert(isapproxv(0 * x, n));
    assert(isapproxv(x + y, y + x));
    assert(isapproxv((x + y) + z, x + (y + z)));
    assert(isapproxv(+x, x));
    assert(isapproxv(-x, -1 * x));
    assert(isapproxv(x - y, x + (-y)));
    assert(isapproxv(x - x, n));
    assert(isapproxv(x * a, a * x));
    assert(isapproxv(a * (x + y), a * x + a * y));
    assert(isapproxv((a + b) * x, a * x + b * x));

    tensor_alm_t<std::complex<T> > xlm = expand(x);
    tensor_aij_t<T> x2 = evaluate(xlm);
    tensor_alm_t<std::complex<T> > xlm2 = expand(x2);
    tensor_aij_t<T> x3 = evaluate(xlm2);
    assert(isapproxv(x3, x2));
  }

  for (int iter = 0; iter < 100; ++iter) {
    tensor_alm_t<std::complex<T> > n(geom), x(geom), y(geom), z(geom);
    n = 0;
    fmap_(x, [&](auto &a) { a = randc(); });
    fmap_(y, [&](auto &a) { a = randc(); });
    fmap_(z, [&](auto &a) { a = randc(); });
    const std::complex<T> a = randc();
    const std::complex<T> b = randc();

    assert(isapproxv(n, n));
    assert(!isapproxv(n, x));
    assert(isapproxv(n + x, x));
    assert(isapproxv(0 * x, n));
    assert(isapproxv(x + y, y + x));
    assert(isapproxv((x + y) + z, x + (y + z)));
    assert(isapproxv(+x, x));
    assert(isapproxv(-x, -1 * x));
    assert(isapproxv(x - y, x + (-y)));
    assert(isapproxv(x - x, n));
    assert(isapproxv(x * a, a * x));
    assert(isapproxv(a * (x + y), a * x + a * y));
    assert(isapproxv((a + b) * x, a * x + b * x));

    tensor_aij_t<T> xij = evaluate(x);
    tensor_alm_t<std::complex<T> > x2 = expand(xij);
    tensor_aij_t<T> xij2 = evaluate(x2);
    tensor_alm_t<std::complex<T> > x3 = expand(xij2);
    assert(isapproxv(x3, x2));
  }
}

template <typename T> void test_tensor3_aij_alm() {
  CCTK_VINFO("test_tensor3_aij_alm...");

  const int nmodes = 5;
  const geom_t geom(nmodes);

  std::random_device rdev;
  std::default_random_engine reng(rdev());
  std::uniform_real_distribution<T> rdist(-1, 1);
  const auto rand = [&]() { return rdist(reng); };
  const auto randc = [&]() {
    return std::complex<T>(rdist(reng), rdist(reng));
  };

  for (int iter = 0; iter < 100; ++iter) {
    tensor3_aij_t<T> n(geom), x(geom), y(geom), z(geom);
    n = 0;
    fmap_(x, [&](auto &a) { a = rand(); });
    fmap_(y, [&](auto &a) { a = rand(); });
    fmap_(z, [&](auto &a) { a = rand(); });
    const T a = rand();
    const T b = rand();

    assert(isapproxv(n, n));
    assert(!isapproxv(n, x));
    assert(isapproxv(n + x, x));
    assert(isapproxv(0 * x, n));
    assert(isapproxv(x + y, y + x));
    assert(isapproxv((x + y) + z, x + (y + z)));
    assert(isapproxv(+x, x));
    assert(isapproxv(-x, -1 * x));
    assert(isapproxv(x - y, x + (-y)));
    assert(isapproxv(x - x, n));
    assert(isapproxv(x * a, a * x));
    assert(isapproxv(a * (x + y), a * x + a * y));
    assert(isapproxv((a + b) * x, a * x + b * x));

    tensor3_alm_t<std::complex<T> > xlm = expand(x);
    tensor3_aij_t<T> x2 = evaluate(xlm);
    tensor3_alm_t<std::complex<T> > xlm2 = expand(x2);
    tensor3_aij_t<T> x3 = evaluate(xlm2);
    assert(isapproxv(x3, x2));
  }

  for (int iter = 0; iter < 100; ++iter) {
    tensor3_alm_t<std::complex<T> > n(geom), x(geom), y(geom), z(geom);
    n = 0;
    fmap_(x, [&](auto &a) { a = randc(); });
    fmap_(y, [&](auto &a) { a = randc(); });
    fmap_(z, [&](auto &a) { a = randc(); });
    const std::complex<T> a = randc();
    const std::complex<T> b = randc();

    assert(isapproxv(n, n));
    assert(!isapproxv(n, x));
    assert(isapproxv(n + x, x));
    assert(isapproxv(0 * x, n));
    assert(isapproxv(x + y, y + x));
    assert(isapproxv((x + y) + z, x + (y + z)));
    assert(isapproxv(+x, x));
    assert(isapproxv(-x, -1 * x));
    assert(isapproxv(x - y, x + (-y)));
    assert(isapproxv(x - x, n));
    assert(isapproxv(x * a, a * x));
    assert(isapproxv(a * (x + y), a * x + a * y));
    assert(isapproxv((a + b) * x, a * x + b * x));

    tensor3_aij_t<T> xij = evaluate(x);
    tensor3_alm_t<std::complex<T> > x2 = expand(xij);
    tensor3_aij_t<T> xij2 = evaluate(x2);
    tensor3_alm_t<std::complex<T> > x3 = expand(xij2);
    assert(isapproxv(x3, x2));
  }
}

template <typename T> void test_scalar_gradient() {
  CCTK_VINFO("test_scalar_gradient...");

  const int nmodes = 5;
  const geom_t geom(nmodes);

  scalar_aij_t<T> Y0p0ij(geom), Y1m1ij(geom), Y1p0ij(geom), Y1p1ij(geom);
  for (int j = 0; j < geom.nphi; ++j) {
#pragma omp simd
    for (int i = 0; i < geom.ntheta; ++i) {
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);
      using std::cos, std::sin;
      Y0p0ij()(i, j) = 1;
      Y1m1ij()(i, j) = sin(theta) * sin(phi);
      Y1p0ij()(i, j) = cos(theta);
      Y1p1ij()(i, j) = sin(theta) * cos(phi);
    }
  }

  const auto Y0p0lm = expand(Y0p0ij);
  const auto Y1m1lm = expand(Y1m1ij);
  const auto Y1p0lm = expand(Y1p0ij);
  const auto Y1p1lm = expand(Y1p1ij);

  const auto dY0p0lm = gradient(Y0p0lm);
  const auto dY1m1lm = gradient(Y1m1lm);
  const auto dY1p0lm = gradient(Y1p0lm);
  const auto dY1p1lm = gradient(Y1p1lm);

  const auto dY0p0ij = evaluate(dY0p0lm);
  const auto dY1m1ij = evaluate(dY1m1lm);
  const auto dY1p0ij = evaluate(dY1p0lm);
  const auto dY1p1ij = evaluate(dY1p1lm);

  for (int j = 0; j < geom.nphi; ++j) {
    for (int i = 0; i < geom.ntheta; ++i) {
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);
      using std::cos, std::sin;
      // Y0p0ij()(i, j) = 1;
      assert(isapprox(dY0p0ij(0)(i, j), 0));
      assert(isapprox(dY0p0ij(1)(i, j), 0));
      // Y1m1ij()(i, j) = sin(theta) * sin(phi);
      assert(isapprox(dY1m1ij(0)(i, j), cos(theta) * sin(phi)));
      assert(isapprox(dY1m1ij(1)(i, j), cos(phi)));
      // Y1p0ij()(i, j) = cos(theta);
      assert(isapprox(dY1p0ij(0)(i, j), -sin(theta)));
      assert(isapprox(dY1p0ij(1)(i, j), 0));
      // Y1p1ij()(i, j) = sin(theta) * cos(phi);
      assert(isapprox(dY1p1ij(0)(i, j), cos(theta) * cos(phi)));
      assert(isapprox(dY1p1ij(1)(i, j), -sin(phi)));
    }
  }
}

extern "C" void AHFinder_test_discretization(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_test_discretization;

  CCTK_VINFO("Testing discretization methods...");
  test_aij_alm<CCTK_REAL>();
  test_scalar_aij_alm<CCTK_REAL>();
  test_vector_aij_alm<CCTK_REAL>();
  test_tensor_aij_alm<CCTK_REAL>();
  test_tensor3_aij_alm<CCTK_REAL>();
  test_scalar_gradient<CCTK_REAL>();
  CCTK_VINFO("Done.");
}

} // namespace AHFinder
