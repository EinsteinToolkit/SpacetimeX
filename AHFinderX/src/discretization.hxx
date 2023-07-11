#ifndef DISCRETIZATION_HXX
#define DISCRETIZATION_HXX

#include <cctk.h>

#include <ssht/ssht.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace AHFinder {

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct abs_result { using value_type = T; };
template <typename T> struct abs_result<std::complex<T> > {
  using value_type = T;
};
template <typename T> using abs_result_t = typename abs_result<T>::value_type;

////////////////////////////////////////////////////////////////////////////////

struct geom_t {
  const int nmodes;

  // Grid
  const int ntheta;
  const int nphi;
  const int npoints;

  // Coefficients
  const int lmax;
  const int ncoeffs;

  geom_t() = delete;
  geom_t(const int nmodes)
      : nmodes(nmodes), ntheta(nmodes), nphi(2 * nmodes - 1),
        npoints(ntheta * nphi), lmax(nmodes - 1), ncoeffs(nmodes * nmodes) {}

  double coord_theta(const int i, const int j) const {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const double theta = ssht_sampling_mw_t2theta(i, nmodes);
    assert(theta >= 0 && theta <= M_PI);
    return theta;
  }

  double coord_phi(const int i, const int j) const {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const double phi = ssht_sampling_mw_p2phi(j, nmodes);
    assert(phi >= 0 && phi < 2 * M_PI);
    return phi;
  }

  double coord_dtheta(const int i, const int j) const {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const double theta_m =
        i == 0 ? 0.0 : (coord_theta(i - 1, j) + coord_theta(i, j)) / 2;
    const double theta_p =
        i == ntheta - 1 ? M_PI
                        : (coord_theta(i, j) + coord_theta(i + 1, j)) / 2;
    return theta_p - theta_m;
  }

  double coord_dphi(const int i, const int j) const { return 2 * M_PI / nphi; }

  int gind(const int i, const int j) const {
    // 0 <= i < ntheta
    // 0 <= j < nphi
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const int ind = j + nphi * i;
    assert(ind >= 0 && ind < npoints);
    return ind;
  }

  int cind(const int l, const int m) const {
    // 0 <= l <= lmax
    // -l <= m <= l
    assert(l >= 0 && l <= lmax);
    assert(m >= -l && m <= l);
    int ind;
    ssht_sampling_elm2ind(&ind, l, m);
    assert(ind >= 0 && ind < ncoeffs);
    return ind;
  }

  friend std::ostream &operator<<(std::ostream &os, const geom_t &geom) {
    return os << "geom_t{ntheta:" << geom.ntheta << ",nphi:" << geom.nphi
              << ",npoints:" << geom.npoints << ",lmax:" << geom.lmax
              << ",ncoeffs:" << geom.ncoeffs << "}";
  }
};

////////////////////////////////////////////////////////////////////////////////

// Provide std::vector space operations for type `V<T>`
template <template <typename> typename V, typename T> struct vectorspace_mixin {
  using value_type = T;

  template <typename U> V<T> &operator=(const U &a) {
    fmap_(*(V<T> *)this, [&](auto &x) { x = a; });
    return *(V<T> *)this;
  }
  template <typename U> V<T> &operator+=(const U &a) {
    fmap_(*(V<T> *)this, [&](auto &x) { x += a; });
    return *(V<T> *)this;
  }
  template <typename U> V<T> &operator-=(const U &a) {
    fmap_(*(V<T> *)this, [&](auto &x) { x -= a; });
    return *(V<T> *)this;
  }
  template <typename U> V<T> &operator*=(const U &a) {
    fmap_(*(V<T> *)this, [&](auto &x) { x *= a; });
    return *(V<T> *)this;
  }
  template <typename U> V<T> &operator/=(const U &a) {
    fmap_(*(V<T> *)this, [&](auto &x) { x /= a; });
    return *(V<T> *)this;
  }

  friend V<T> operator+(const V<T> &xs) {
    return fmap([](auto x) { return +x; }, xs);
  }
  friend V<T> operator-(const V<T> &xs) {
    return fmap([](auto x) { return -x; }, xs);
  }
  friend V<T> abs(const V<T> &xs) {
    return fmap(
        [](auto x) {
          using std::abs;
          return abs(x);
        },
        xs);
  }
  friend V<T> conj(const V<T> &xs) {
    return fmap(
        [](auto x) {
          using std::conj;
          return conj(x);
        },
        xs);
  }
  friend V<T> operator+(const V<T> &xs, const V<T> &ys) {
    return fmap([](auto x, auto y) { return x + y; }, xs, ys);
  }
  friend V<T> operator-(const V<T> &xs, const V<T> &ys) {
    return fmap([](auto x, auto y) { return x - y; }, xs, ys);
  }
  friend V<T> operator*(const T &x, const V<T> &ys) {
    return fmap([&](auto y) { return x * y; }, ys);
  }
  friend V<T> operator*(const V<T> &xs, const T &y) {
    return fmap([&](auto x) { return x * y; }, xs);
  }
  friend V<T> operator/(const V<T> &xs, const T &y) {
    return fmap([&](auto x) { return x / y; }, xs);
  }
  friend V<T> max(const V<T> &xs, const V<T> &ys) {
    return fmap(
        [](auto x, auto y) {
          using std::max;
          return max(x, y);
        },
        xs, ys);
  }
  friend V<T> min(const V<T> &xs, const V<T> &ys) {
    return fmap(
        [](auto x, auto y) {
          using std::min;
          return min(x, y);
        },
        xs, ys);
  }

  using maxabs_t = abs_result_t<T>;
  friend maxabs_t maxabs(const V<T> &xs) {
    return fmapreduce(
        [](auto x) {
          using std::abs;
          return abs(x);
        },
        [](auto x, auto y) {
          using std::max;
          return max(x, y);
        },
        maxabs_t(0), xs);
  }
  friend T maximum(const V<T> &xs) {
    return fmapreduce([](auto x) { return x; },
                      [](auto x, auto y) {
                        using std::max;
                        return max(x, y);
                      },
                      std::numeric_limits<T>::lowest(), xs);
  }
  friend T minimum(const V<T> &xs) {
    return fmapreduce([](auto x) { return x; },
                      [](auto x, auto y) {
                        using std::min;
                        return min(x, y);
                      },
                      std::numeric_limits<T>::max(), xs);
  }

  friend bool operator==(const V<T> &xs, const V<T> &ys) {
    return fmapreduce([&](auto x, auto y) { return x == y; },
                      [&](auto x, auto y) { return x && y; }, true, xs, ys);
  }
  friend bool operator!=(const V<T> &xs, const V<T> &ys) { return !(xs == ys); }
};

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct alm_t : vectorspace_mixin<alm_t, T> {
  const geom_t &geom;

  int spin;
  std::vector<T> alm;

  alm_t() = delete;
  alm_t(const geom_t &geom, int spin)
      : geom(geom), spin(spin), alm(geom.ncoeffs) {}

  alm_t(const alm_t &) = default;
  alm_t(alm_t &&) = default;
  alm_t &operator=(const alm_t &other) {
    assert(&geom == &other.geom);
    spin = other.spin;
    alm = other.alm;
    return *this;
  }
  alm_t &operator=(alm_t &&other) {
    assert(&geom == &other.geom);
    spin = other.spin;
    alm = std::move(other.alm);
    return *this;
  }

  const T *data() const { return alm.data(); }
  T *data() { return alm.data(); }
  const T &operator()(const int l, const int m) const {
    return alm.at(geom.cind(l, m));
  }
  T &operator()(const int l, const int m) { return alm.at(geom.cind(l, m)); }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const alm_t &xs,
                      const alm_t<Args> &...ys) {
    using std::all_of;
    const std::array<bool, sizeof...(Args)> check_ys{&ys.geom == &xs.geom &&
                                                     ys.spin == xs.spin...};
    assert(all_of(check_ys.begin(), check_ys.end(), [](auto x) { return x; }));
#pragma omp simd
    for (std::size_t n = 0; n < xs.alm.size(); ++n)
      r = op(r, f(xs.alm[n], ys.alm[n]...));
    return r;
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend alm_t<R> fmap(const F &f, const alm_t &xs, const alm_t<Args> &...ys) {
    using std::all_of;
    const std::array<bool, sizeof...(Args)> check_ys{&ys.geom == &xs.geom &&
                                                     ys.spin == xs.spin...};
    assert(all_of(check_ys.begin(), check_ys.end(), [](auto x) { return x; }));
    alm_t<R> r(xs.geom, xs.spin);
#pragma omp simd
    for (std::size_t n = 0; n < xs.alm.size(); ++n)
      r.alm[n] = f(xs.alm[n], ys.alm[n]...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(alm_t &xs, const F &f, const alm_t<Args> &...ys) {
    using std::all_of;
    const std::array<bool, sizeof...(Args)> check_ys{&ys.geom == &xs.geom &&
                                                     ys.spin == xs.spin...};
    assert(all_of(check_ys.begin(), check_ys.end(), [](auto x) { return x; }));
#pragma omp simd
    for (std::size_t n = 0; n < xs.alm.size(); ++n)
      f(xs.alm[n], ys.alm[n]...);
  }

  // template <typename U> alm_t &operator=(const U &x) {
  //   return (vectorspace_mixin<alm_t, T> &)(*this) = x;
  // }
  using vectorspace_mixin<alm_t, T>::operator=;
  using vectorspace_mixin<alm_t, T>::operator+=;
  using vectorspace_mixin<alm_t, T>::operator-=;
  using vectorspace_mixin<alm_t, T>::operator*=;
  using vectorspace_mixin<alm_t, T>::operator/=;

  friend std::ostream &operator<<(std::ostream &os, const alm_t &alm) {
    os << "alm_t{geom:" << alm.geom << ",spin:" << alm.spin << ",coeffs:[";
    for (auto x : alm.alm)
      os << x << ",";
    os << "]}";
    return os;
  }
};

template <typename T> struct aij_t : vectorspace_mixin<aij_t, T> {
  const geom_t &geom;

  std::vector<T> aij;

  aij_t() = delete;
  aij_t(const geom_t &geom) : geom(geom), aij(geom.npoints, NAN) {}

  aij_t(const aij_t &) = default;
  aij_t(aij_t &&) = default;
  aij_t &operator=(const aij_t &other) {
    assert(&geom == &other.geom);
    aij = other.aij;
    return *this;
  }
  aij_t &operator=(aij_t &&other) {
    assert(&geom == &other.geom);
    aij = std::move(other.aij);
    return *this;
  }

  const T *data() const { return aij.data(); }
  T *data() { return aij.data(); }
  const T &operator()(const int i, const int j) const {
    return aij.at(geom.gind(i, j));
  }
  T &operator()(const int i, const int j) { return aij.at(geom.gind(i, j)); }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const aij_t &xs,
                      const aij_t<Args> &...ys) {
    using std::all_of;
    const std::array<bool, sizeof...(Args)> check_ys{&ys.geom == &xs.geom...};
    assert(all_of(check_ys.begin(), check_ys.end(), [](auto x) { return x; }));
#pragma omp simd
    for (std::size_t n = 0; n < xs.aij.size(); ++n)
      r = op(r, f(xs.aij[n], ys.aij[n]...));
    return r;
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend aij_t<R> fmap(const F &f, const aij_t &xs, const aij_t<Args> &...ys) {
    using std::all_of;
    const std::array<bool, sizeof...(Args)> check_ys{&ys.geom == &xs.geom...};
    assert(all_of(check_ys.begin(), check_ys.end(), [](auto x) { return x; }));
    aij_t<R> r(xs.geom);
#pragma omp simd
    for (std::size_t n = 0; n < xs.aij.size(); ++n)
      r.aij[n] = f(xs.aij[n], ys.aij[n]...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(aij_t &xs, const F &f, const aij_t<Args> &...ys) {
    using std::all_of;
    const std::array<bool, sizeof...(Args)> check_ys{&ys.geom == &xs.geom...};
    assert(all_of(check_ys.begin(), check_ys.end(), [](auto x) { return x; }));
#pragma omp simd
    for (std::size_t n = 0; n < xs.aij.size(); ++n)
      f(xs.aij[n], ys.aij[n]...);
  }

  // template <typename U> aij_t &operator=(const U &x) {
  //   return (vectorspace_mixin<aij_t, T> &)(*this) = x;
  // }
  using vectorspace_mixin<aij_t, T>::operator=;
  using vectorspace_mixin<aij_t, T>::operator+=;
  using vectorspace_mixin<aij_t, T>::operator-=;
  using vectorspace_mixin<aij_t, T>::operator*=;
  using vectorspace_mixin<aij_t, T>::operator/=;

  friend std::ostream &operator<<(std::ostream &os, const aij_t &aij) {
    os << "aij_t{geom:" << aij.geom << ",values:[";
    for (auto x : aij.aij)
      os << x << ",";
    os << "]}";
    return os;
  }
};

template <typename T>
alm_t<T> coefficients_from_const(const geom_t &geom, const T r0,
                                 const T r1z = 0) {
  alm_t<T> alm(geom, 0);
  using std::sqrt;
  for (int l = 0; l <= geom.lmax; ++l)
#pragma omp simd
    for (int m = -l; m <= l; ++m)
      alm(l, m) = l == 0 && m == 0 ? sqrt(T(4 * M_PI)) * r0 : 0;
  alm(1, 0) = T(2) * r1z;
  return alm;
}

template <typename T> T average(const alm_t<std::complex<T> > &alm) {
  using std::sqrt;
  return real(alm(0, 0)) / sqrt(4 * T(M_PI));
}

template <typename T>
alm_t<std::complex<T> > expand(const aij_t<T> &aij, const int spin) {
  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]
  alm_t<std::complex<T> > alm(aij.geom, spin);
  ssht_core_mw_forward_sov_conv_sym_real(alm.data(), aij.data(),
                                         alm.geom.nmodes, method, verbosity);
  return alm;
}

template <typename T>
alm_t<std::complex<T> > expand(const aij_t<std::complex<T> > &aij,
                               const int spin) {
  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]
  alm_t<std::complex<T> > alm(aij.geom, spin);
  ssht_core_mw_forward_sov_conv_sym(alm.data(), aij.data(), alm.geom.nmodes,
                                    spin, method, verbosity);
  return alm;
}

template <typename T> aij_t<T> evaluate(const alm_t<std::complex<T> > &alm) {
  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]
  aij_t<T> aij(alm.geom);
  ssht_core_mw_inverse_sov_sym_real(aij.data(), alm.data(), aij.geom.nmodes,
                                    method, verbosity);
  return aij;
}

template <typename T>
aij_t<std::complex<T> > evaluate(const alm_t<std::complex<T> > &alm,
                                 const int spin) {
  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]
  const geom_t &geom = alm.geom;
  aij_t<std::complex<T> > aij(geom);
  ssht_core_mw_inverse_sov_sym(aij.data(), alm.data(), aij.geom.nmodes, spin,
                               method, verbosity);
  return aij;
}

// \dh
template <typename T>
alm_t<std::complex<T> > eth(const alm_t<std::complex<T> > &alm) {
  const geom_t &geom = alm.geom;
  const int s = alm.spin;
  alm_t<std::complex<T> > blm(geom, s + 1);
  using std::sqrt;
  for (int l = 0; l <= geom.lmax; ++l)
#pragma omp simd
    for (int m = -l; m <= l; ++m)
      blm(l, m) = +sqrt(T((l - s) * (l + s + 1))) * alm(l, m);
  return blm;
}

// \bar\dh
template <typename T>
alm_t<std::complex<T> > eth_bar(const alm_t<std::complex<T> > &alm) {
  const geom_t &geom = alm.geom;
  const int s = alm.spin;
  alm_t<std::complex<T> > blm(geom, s - 1);
  using std::sqrt;
  for (int l = 0; l <= geom.lmax; ++l)
#pragma omp simd
    for (int m = -l; m <= l; ++m)
      blm(l, m) = -sqrt(T((l + s) * (l - s + 1))) * alm(l, m);
  return blm;
}

////////////////////////////////////////////////////////////////////////////////

// e_\theta, e_\phi are orthonormal dyad std::vectors
//
//      m = e_\theta + im e_\phi
// \bar m = e_\theta - im e_\phi
//
//      m^a      m_a = 0
// \bar m^a \bar m_a = 0
//      m^a \bar m_a = 2
// \bar m^a      m_a = 2
//
// e_\theta = 1/2   (m + \bar m)
// e_\phi   = 1/2im (m - \bar m)

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct scalar_aij_t : vectorspace_mixin<scalar_aij_t, T> {
  const geom_t &geom;

  aij_t<T> elts;

  scalar_aij_t() = delete;
  scalar_aij_t(const geom_t &geom) : geom(geom), elts(geom) {}
  scalar_aij_t(const scalar_aij_t &) = default;
  scalar_aij_t(scalar_aij_t &&) = default;
  scalar_aij_t &operator=(const scalar_aij_t &other) {
    assert(&geom == &other.geom);
    elts = other.elts;
    return *this;
  }
  scalar_aij_t &operator=(scalar_aij_t &&other) {
    assert(&geom == &other.geom);
    elts = std::move(other.elts);
    return *this;
  }

  const aij_t<T> &operator()() const { return elts; }
  aij_t<T> &operator()() { return elts; }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const scalar_aij_t &xs,
                      const scalar_aij_t<Args> &...ys) {
    return fmapreduce(f, op, r, xs(), ys()...);
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend scalar_aij_t<R> fmap(const F &f, const scalar_aij_t &xs,
                              const scalar_aij_t<Args> &...ys) {
    scalar_aij_t<R> r(xs.geom);
    fmap_(
        r(), [&](R &r, T x, Args... y) { r = f(x, y...); }, xs(), ys()...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(scalar_aij_t &xs, const F &f,
                    const scalar_aij_t<Args> &...ys) {
    fmap_(xs(), f, ys()...);
  }

  using vectorspace_mixin<scalar_aij_t, T>::operator=;
  using vectorspace_mixin<scalar_aij_t, T>::operator+=;
  using vectorspace_mixin<scalar_aij_t, T>::operator-=;
  using vectorspace_mixin<scalar_aij_t, T>::operator*=;
  using vectorspace_mixin<scalar_aij_t, T>::operator/=;

  friend std::ostream &operator<<(std::ostream &os, const scalar_aij_t &saij) {
    return os << "scalar_aij_t{geom:" << saij.geom << ",elts:" << saij.elts
              << "}";
  }
};

template <typename T> struct scalar_alm_t : vectorspace_mixin<scalar_alm_t, T> {
  const geom_t &geom;

  alm_t<T> elts;

  scalar_alm_t() = delete;
  scalar_alm_t(const geom_t &geom) : geom(geom), elts(geom, 0) {}
  scalar_alm_t(const scalar_alm_t &) = default;
  scalar_alm_t(scalar_alm_t &&) = default;
  scalar_alm_t &operator=(const scalar_alm_t &other) {
    assert(&geom == &other.geom);
    elts = other.elts;
    return *this;
  }
  scalar_alm_t &operator=(scalar_alm_t &&other) {
    assert(&geom == &other.geom);
    elts = std::move(other.elts);
    return *this;
  }

  const alm_t<T> &operator()() const { return elts; }
  alm_t<T> &operator()() { return elts; }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const scalar_alm_t &xs,
                      const scalar_alm_t<Args> &...ys) {
    return fmapreduce(f, op, r, xs(), ys()...);
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend scalar_alm_t<R> fmap(const F &f, const scalar_alm_t &xs,
                              const scalar_alm_t<Args> &...ys) {
    scalar_alm_t<R> r(xs.geom);
    fmap_(
        r(), [&](R &r, T x, Args... y) { r = f(x, y...); }, xs(), ys()...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(scalar_alm_t &xs, const F &f,
                    const scalar_alm_t<Args> &...ys) {
    fmap_(xs(), f, ys()...);
  }

  using vectorspace_mixin<scalar_alm_t, T>::operator=;
  using vectorspace_mixin<scalar_alm_t, T>::operator+=;
  using vectorspace_mixin<scalar_alm_t, T>::operator-=;
  using vectorspace_mixin<scalar_alm_t, T>::operator*=;
  using vectorspace_mixin<scalar_alm_t, T>::operator/=;

  friend std::ostream &operator<<(std::ostream &os, const scalar_alm_t &salm) {
    return os << "scalar_alm_t{geom:" << salm.geom << ",elts:" << salm.elts
              << "}";
  }
};

template <typename T>
scalar_alm_t<T> scalar_from_const(const geom_t &geom, const T r0,
                                  const T r1z = 0) {
  scalar_alm_t<T> salm(geom);
  salm() = coefficients_from_const(geom, r0, r1z);
  return salm;
}

template <typename T> T average(const scalar_alm_t<std::complex<T> > &alm) {
  return average(alm());
}

template <typename T>
scalar_alm_t<std::complex<T> > expand(const scalar_aij_t<T> &saij) {
  const geom_t &geom = saij.geom;

  scalar_aij_t<std::complex<T> > sbij(geom);
  for (int j = 0; j < geom.nphi; ++j)
#pragma omp simd
    for (int i = 0; i < geom.ntheta; ++i)
      sbij()(i, j) = saij()(i, j);

  scalar_alm_t<std::complex<T> > salm(geom);
  salm() = expand(sbij(), 0);

  return salm;
}

template <typename T>
scalar_aij_t<T> evaluate(const scalar_alm_t<std::complex<T> > &salm) {
  const geom_t &geom = salm.geom;

  scalar_aij_t<std::complex<T> > sbij(geom);
  sbij() = evaluate(salm(), 0);

  scalar_aij_t<T> saij(geom);
  for (int j = 0; j < geom.nphi; ++j)
#pragma omp simd
    for (int i = 0; i < geom.ntheta; ++i)
      saij()(i, j) = real(sbij()(i, j));

  return saij;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct vector_aij_t : vectorspace_mixin<vector_aij_t, T> {
  const geom_t &geom;

  std::array<aij_t<T>, 2> elts;

  vector_aij_t() = delete;
  vector_aij_t(const geom_t &geom)
      : geom(geom), elts{aij_t<T>(geom), aij_t<T>(geom)} {}
  vector_aij_t(const vector_aij_t &) = default;
  vector_aij_t(vector_aij_t &&) = default;
  vector_aij_t &operator=(const vector_aij_t &other) {
    assert(&geom == &other.geom);
    elts = other.elts;
    return *this;
  }
  vector_aij_t &operator=(vector_aij_t &&other) {
    assert(&geom == &other.geom);
    elts = std::move(other.elts);
    return *this;
  }

  const aij_t<T> &operator()(int i) const {
    assert(i >= 0 && i < 2);
    return elts.at(i);
  }
  aij_t<T> &operator()(int i) {
    assert(i >= 0 && i < 2);
    return elts.at(i);
  }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const vector_aij_t &xs,
                      const vector_aij_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      r = fmapreduce(f, op, r, xs(i), ys(i)...);
    return r;
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend vector_aij_t<R> fmap(const F &f, const vector_aij_t &xs,
                              const vector_aij_t<Args> &...ys) {
    vector_aij_t<R> r(xs.geom);
    for (int i = 0; i < 2; ++i)
      fmap_(
          r(i), [&](R &r, T x, Args... y) { r = f(x, y...); }, xs(i), ys(i)...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(vector_aij_t &xs, const F &f,
                    const vector_aij_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      fmap_(xs(i), f, ys(i)...);
  }

  using vectorspace_mixin<vector_aij_t, T>::operator=;
  using vectorspace_mixin<vector_aij_t, T>::operator+=;
  using vectorspace_mixin<vector_aij_t, T>::operator-=;
  using vectorspace_mixin<vector_aij_t, T>::operator*=;
  using vectorspace_mixin<vector_aij_t, T>::operator/=;

  friend std::ostream &operator<<(std::ostream &os, const vector_aij_t &vaij) {
    return os << "vector_aij_t{geom:" << vaij.geom << ",elts:["
              << vaij.elts.at(0) << "," << vaij.elts.at(1) << "]}";
  }
};

template <typename T> struct vector_alm_t : vectorspace_mixin<vector_alm_t, T> {
  const geom_t &geom;

  std::array<alm_t<T>, 2> elts;

  vector_alm_t() = delete;
  vector_alm_t(const geom_t &geom)
      : geom(geom), elts{alm_t<T>(geom, +1), alm_t<T>(geom, -1)} {}
  vector_alm_t(const vector_alm_t &) = default;
  vector_alm_t(vector_alm_t &&) = default;
  vector_alm_t &operator=(const vector_alm_t &other) {
    assert(&geom == &other.geom);
    elts = other.elts;
    return *this;
  }
  vector_alm_t &operator=(vector_alm_t &&other) {
    assert(&geom == &other.geom);
    elts = std::move(other.elts);
    return *this;
  }

  const alm_t<T> &operator()(int i) const { return elts.at(i); }
  alm_t<T> &operator()(int i) { return elts.at(i); }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const vector_alm_t &xs,
                      const vector_alm_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      r = fmapreduce(f, op, r, xs(i), ys(i)...);
    return r;
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend vector_alm_t<R> fmap(const F &f, const vector_alm_t &xs,
                              const vector_alm_t<Args> &...ys) {
    vector_alm_t<R> r(xs.geom);
    for (int i = 0; i < 2; ++i)
      fmap_(
          r(i), [&](R &r, T x, Args... y) { r = f(x, y...); }, xs(i), ys(i)...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(vector_alm_t &xs, const F &f,
                    const vector_alm_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      fmap_(xs(i), f, ys(i)...);
  }

  using vectorspace_mixin<vector_alm_t, T>::operator=;
  using vectorspace_mixin<vector_alm_t, T>::operator+=;
  using vectorspace_mixin<vector_alm_t, T>::operator-=;
  using vectorspace_mixin<vector_alm_t, T>::operator*=;
  using vectorspace_mixin<vector_alm_t, T>::operator/=;

  friend std::ostream &operator<<(std::ostream &os, const vector_alm_t &valm) {
    return os << "vector_alm_t{geom:" << valm.geom << ",elts:["
              << valm.elts.at(0) << "," << valm.elts.at(1) << "]}";
  }
};

template <typename T>
vector_alm_t<std::complex<T> > expand(const vector_aij_t<T> &vaij) {
  using namespace std::literals::complex_literals;

  const geom_t &geom = vaij.geom;

  vector_aij_t<std::complex<T> > vbij(geom);
  for (int j = 0; j < geom.nphi; ++j) {
#pragma omp simd
    for (int i = 0; i < geom.ntheta; ++i) {
      // b0 =      m_x a^x
      // b1 = \bar m_x a^x

      // m[p] = {m, \bar m}
      const std::array<std::array<std::complex<T>, 2>, 2> m{
          {{{1, +1i}}, {{1, -1i}}}};
      std::array<std::array<std::complex<T>, 2>, 2> b;
      for (int p = 0; p < 2; ++p) {
        std::complex<T> s = 0;
        for (int x = 0; x < 2; ++x)
          s += m[p][x] * vaij(x)(i, j);
        vbij(p)(i, j) = s;
      }
    }
  }

  vector_alm_t<std::complex<T> > valm(geom);
  valm(0) = expand(vbij(0), +1);
  valm(1) = expand(vbij(1), -1);

  return valm;
}

template <typename T>
vector_aij_t<T> evaluate(const vector_alm_t<std::complex<T> > &valm) {
  using namespace std::literals::complex_literals;

  const geom_t &geom = valm.geom;

  vector_aij_t<std::complex<T> > vbij(geom);
  vbij(0) = evaluate(valm(0), +1);
  vbij(1) = evaluate(valm(1), -1);

  vector_aij_t<T> vaij(geom);
  for (int j = 0; j < geom.nphi; ++j) {
#pragma omp simd
    for (int i = 0; i < geom.ntheta; ++i) {
      // b0 =      m_x a^x
      // b1 = \bar m_x a^x

      // m[p] = {m, \bar m}
      const std::array<std::array<std::complex<T>, 2>, 2> minv2{
          {{{+1, +1}}, {{-1i, +1i}}}}; // minv2[x][p] = 2 * inv(m)
      for (int x = 0; x < 2; ++x) {
        std::complex<T> s = 0;
        for (int p = 0; p < 2; ++p)
          s += minv2[x][p] * vbij(p)(i, j);
        vaij(x)(i, j) = real(s) / 2;
      }
    }
  }

  return vaij;
}

template <typename T>
vector_alm_t<std::complex<T> >
gradient(const scalar_alm_t<std::complex<T> > &salm) {
  const geom_t &geom = salm.geom;

  vector_alm_t<std::complex<T> > valm(geom);

  valm(0) = -eth(salm());
  valm(1) = -eth_bar(salm());

  return valm;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct tensor_aij_t : vectorspace_mixin<tensor_aij_t, T> {
  const geom_t &geom;

  std::array<std::array<aij_t<T>, 2>, 2> elts;

  tensor_aij_t() = delete;
  tensor_aij_t(const geom_t &geom)
      : geom(geom), elts{{{{aij_t<T>(geom), aij_t<T>(geom)}},
                          {{aij_t<T>(geom), aij_t<T>(geom)}}}} {}
  tensor_aij_t(const tensor_aij_t &) = default;
  tensor_aij_t(tensor_aij_t &&) = default;
  tensor_aij_t &operator=(const tensor_aij_t &other) {
    assert(&geom == &other.geom);
    elts = other.elts;
    return *this;
  }
  tensor_aij_t &operator=(tensor_aij_t &&other) {
    assert(&geom == &other.geom);
    elts = std::move(other.elts);
    return *this;
  }

  const aij_t<T> &operator()(int i, int j) const {
    assert(i >= 0 && i < 2);
    assert(j >= 0 && j < 2);
    return elts.at(i).at(j);
  }
  aij_t<T> &operator()(int i, int j) {
    assert(i >= 0 && i < 2);
    assert(j >= 0 && j < 2);
    return elts.at(i).at(j);
  }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const tensor_aij_t &xs,
                      const tensor_aij_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        r = fmapreduce(f, op, r, xs(i, j), ys(i, j)...);
    return r;
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend tensor_aij_t<R> fmap(const F &f, const tensor_aij_t &xs,
                              const tensor_aij_t<Args> &...ys) {
    tensor_aij_t<R> r(xs.geom);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        fmap_(
            r(i, j), [&](R &r, T x, Args... y) { r = f(x, y...); }, xs(i, j),
            ys(i, j)...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(tensor_aij_t &xs, const F &f,
                    const tensor_aij_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        fmap_(xs(i, j), f, ys(i, j)...);
  }

  using vectorspace_mixin<tensor_aij_t, T>::operator=;
  using vectorspace_mixin<tensor_aij_t, T>::operator+=;
  using vectorspace_mixin<tensor_aij_t, T>::operator-=;
  using vectorspace_mixin<tensor_aij_t, T>::operator*=;
  using vectorspace_mixin<tensor_aij_t, T>::operator/=;
};

template <typename T> struct tensor_alm_t : vectorspace_mixin<tensor_alm_t, T> {
  const geom_t &geom;

  std::array<std::array<alm_t<T>, 2>, 2> elts;

  tensor_alm_t() = delete;
  tensor_alm_t(const geom_t &geom)
      : geom(geom), elts{{{{alm_t<T>(geom, +2), alm_t<T>(geom, 0)}},
                          {{alm_t<T>(geom, 0), alm_t<T>(geom, -2)}}}} {}
  tensor_alm_t(const tensor_alm_t &) = default;
  tensor_alm_t(tensor_alm_t &&) = default;
  tensor_alm_t &operator=(const tensor_alm_t &other) {
    assert(&geom == &other.geom);
    elts = other.elts;
    return *this;
  }
  tensor_alm_t &operator=(tensor_alm_t &&other) {
    assert(&geom == &other.geom);
    elts = std::move(other.elts);
    return *this;
  }

  const alm_t<T> &operator()(int i, int j) const { return elts.at(i).at(j); }
  alm_t<T> &operator()(int i, int j) { return elts.at(i).at(j); }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const tensor_alm_t &xs,
                      const tensor_alm_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        r = fmapreduce(f, op, r, xs(i, j), ys(i, j)...);
    return r;
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend tensor_alm_t<R> fmap(const F &f, const tensor_alm_t &xs,
                              const tensor_alm_t<Args> &...ys) {
    tensor_alm_t<R> r(xs.geom);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        fmap_(
            r(i, j), [&](R &r, T x, Args... y) { r = f(x, y...); }, xs(i, j),
            ys(i, j)...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(tensor_alm_t &xs, const F &f,
                    const tensor_alm_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        fmap_(xs(i, j), f, ys(i, j)...);
  }

  using vectorspace_mixin<tensor_alm_t, T>::operator=;
  using vectorspace_mixin<tensor_alm_t, T>::operator+=;
  using vectorspace_mixin<tensor_alm_t, T>::operator-=;
  using vectorspace_mixin<tensor_alm_t, T>::operator*=;
  using vectorspace_mixin<tensor_alm_t, T>::operator/=;
};

template <typename T>
tensor_alm_t<std::complex<T> > expand(const tensor_aij_t<T> &taij) {
  using namespace std::literals::complex_literals;

  const geom_t &geom = taij.geom;

  tensor_aij_t<std::complex<T> > tbij(geom);
  for (int j = 0; j < geom.nphi; ++j) {
#pragma omp simd
    for (int i = 0; i < geom.ntheta; ++i) {
      // b00 =      m_x      m_y a^xy
      // b01 =      m_x \bar m_y a^xy
      // b10 = \bar m_x      m_y a^xy
      // b11 = \bar m_x \bar m_y a^xy

      // m[p] = {m, \bar m}
      const std::array<std::array<std::complex<T>, 2>, 2> m{
          {{{1, +1i}}, {{1, -1i}}}};
      for (int q = 0; q < 2; ++q) {
        for (int p = 0; p < 2; ++p) {
          std::complex<T> s = 0;
          for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
              s += m[p][x] * m[q][y] * taij(x, y)(i, j);
          tbij(p, q)(i, j) = s;
        }
      }
    }
  }

  tensor_alm_t<std::complex<T> > talm(geom);
  talm(0, 0) = expand(tbij(0, 0), +2);
  talm(0, 1) = expand(tbij(0, 1), 0);
  talm(1, 0) = expand(tbij(1, 0), 0);
  talm(1, 1) = expand(tbij(1, 1), -2);

  return talm;
}

template <typename T>
tensor_aij_t<T> evaluate(const tensor_alm_t<std::complex<T> > &talm) {
  using namespace std::literals::complex_literals;

  const geom_t &geom = talm.geom;

  tensor_aij_t<std::complex<T> > tbij(geom);
  tbij(0, 0) = evaluate(talm(0, 0), +2);
  tbij(0, 1) = evaluate(talm(0, 1), 0);
  tbij(1, 0) = evaluate(talm(1, 0), 0);
  tbij(1, 1) = evaluate(talm(1, 1), -2);

  tensor_aij_t<T> taij(geom);
  for (int j = 0; j < geom.nphi; ++j) {
#pragma omp simd
    for (int i = 0; i < geom.ntheta; ++i) {
      // b00 =      m_x      m_y a^xy
      // b01 =      m_x \bar m_y a^xy
      // b10 = \bar m_x      m_y a^xy
      // b11 = \bar m_x \bar m_y a^xy

      // m[p] = {m, \bar m}
      const std::array<std::array<std::complex<T>, 2>, 2> minv2{
          {{{+1, +1}}, {{-1i, +1i}}}}; // minv2[x][p] = 2 * inv(m)
      for (int y = 0; y < 2; ++y) {
        for (int x = 0; x < 2; ++x) {
          std::complex<T> s = 0;
          for (int q = 0; q < 2; ++q)
            for (int p = 0; p < 2; ++p)
              s += minv2[x][p] * minv2[y][q] * tbij(p, q)(i, j);
          taij(x, y)(i, j) = real(s) / 4;
        }
      }
    }
  }

  return taij;
}

template <typename T>
tensor_alm_t<std::complex<T> >
gradient(const vector_alm_t<std::complex<T> > &valm) {
  const geom_t &geom = valm.geom;

  tensor_alm_t<std::complex<T> > talm(geom);

  talm(0, 0) = -eth(valm(0));
  talm(0, 1) = -ethbar(valm(0));
  talm(1, 0) = -eth(valm(1));
  talm(1, 1) = -ethbar(valm(1));

  return valm;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct tensor3_aij_t : vectorspace_mixin<tensor3_aij_t, T> {
  const geom_t &geom;

  std::array<std::array<std::array<aij_t<T>, 2>, 2>, 2> elts;

  tensor3_aij_t() = delete;
  tensor3_aij_t(const geom_t &geom)
      : geom(geom), elts{{{{{{aij_t<T>(geom), aij_t<T>(geom)}},
                            {{aij_t<T>(geom), aij_t<T>(geom)}}}},
                          {{{{aij_t<T>(geom), aij_t<T>(geom)}},
                            {{aij_t<T>(geom), aij_t<T>(geom)}}}}}} {}
  tensor3_aij_t(const tensor3_aij_t &) = default;
  tensor3_aij_t(tensor3_aij_t &&) = default;
  tensor3_aij_t &operator=(const tensor3_aij_t &other) {
    assert(&geom == &other.geom);
    elts = other.elts;
    return *this;
  }
  tensor3_aij_t &operator=(tensor3_aij_t &&other) {
    assert(&geom == &other.geom);
    elts = std::move(other.elts);
    return *this;
  }

  const aij_t<T> &operator()(int i, int j, int k) const {
    assert(i >= 0 && i < 2);
    assert(j >= 0 && j < 2);
    assert(k >= 0 && k < 2);
    return elts.at(i).at(j).at(k);
  }
  aij_t<T> &operator()(int i, int j, int k) {
    assert(i >= 0 && i < 2);
    assert(j >= 0 && j < 2);
    assert(k >= 0 && k < 2);
    return elts.at(i).at(j).at(k);
  }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const tensor3_aij_t &xs,
                      const tensor3_aij_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          r = fmapreduce(f, op, r, xs(i, j, k), ys(i, j, k)...);
    return r;
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend tensor3_aij_t<R> fmap(const F &f, const tensor3_aij_t &xs,
                               const tensor3_aij_t<Args> &...ys) {
    tensor3_aij_t<R> r(xs.geom);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          fmap_(
              r(i, j, k), [&](R &r, T x, Args... y) { r = f(x, y...); },
              xs(i, j, k), ys(i, j, k)...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(tensor3_aij_t &xs, const F &f,
                    const tensor3_aij_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          fmap_(xs(i, j, k), f, ys(i, j, k)...);
  }

  using vectorspace_mixin<tensor3_aij_t, T>::operator=;
  using vectorspace_mixin<tensor3_aij_t, T>::operator+=;
  using vectorspace_mixin<tensor3_aij_t, T>::operator-=;
  using vectorspace_mixin<tensor3_aij_t, T>::operator*=;
  using vectorspace_mixin<tensor3_aij_t, T>::operator/=;
};

template <typename T>
struct tensor3_alm_t : vectorspace_mixin<tensor3_alm_t, T> {
  const geom_t &geom;

  std::array<std::array<std::array<alm_t<T>, 2>, 2>, 2> elts;

  tensor3_alm_t() = delete;
  tensor3_alm_t(const geom_t &geom)
      : geom(geom), elts{{{{{{alm_t<T>(geom, +3), alm_t<T>(geom, +1)}},
                            {{alm_t<T>(geom, +1), alm_t<T>(geom, -1)}}}},
                          {{{{alm_t<T>(geom, +1), alm_t<T>(geom, -1)}},
                            {{alm_t<T>(geom, -1), alm_t<T>(geom, -3)}}}}}} {}
  tensor3_alm_t(const tensor3_alm_t &) = default;
  tensor3_alm_t(tensor3_alm_t &&) = default;
  tensor3_alm_t &operator=(const tensor3_alm_t &other) {
    assert(&geom == &other.geom);
    elts = other.elts;
    return *this;
  }
  tensor3_alm_t &operator=(tensor3_alm_t &&other) {
    assert(&geom == &other.geom);
    elts = std::move(other.elts);
    return *this;
  }

  const alm_t<T> &operator()(int i, int j, int k) const {
    return elts.at(i).at(j).at(k);
  }
  alm_t<T> &operator()(int i, int j, int k) { return elts.at(i).at(j).at(k); }

  template <typename F, typename Op, typename R, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > >,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<Op(R, R)> > > >
  friend R fmapreduce(const F &f, const Op &op, R r, const tensor3_alm_t &xs,
                      const tensor3_alm_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          r = fmapreduce(f, op, r, xs(i, j, k), ys(i, j, k)...);
    return r;
  }
  template <typename F, typename... Args,
            typename R = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T, Args...)> > > >
  friend tensor3_alm_t<R> fmap(const F &f, const tensor3_alm_t &xs,
                               const tensor3_alm_t<Args> &...ys) {
    tensor3_alm_t<R> r(xs.geom);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          fmap_(
              r(i, j, k), [&](R &r, T x, Args... y) { r = f(x, y...); },
              xs(i, j, k), ys(i, j, k)...);
    return r;
  }
  template <typename F, typename... Args,
            typename = std::remove_cv_t<
                std::remove_reference_t<std::result_of_t<F(T &, Args...)> > > >
  friend void fmap_(tensor3_alm_t &xs, const F &f,
                    const tensor3_alm_t<Args> &...ys) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          fmap_(xs(i, j, k), f, ys(i, j, k)...);
  }

  using vectorspace_mixin<tensor3_alm_t, T>::operator=;
  using vectorspace_mixin<tensor3_alm_t, T>::operator+=;
  using vectorspace_mixin<tensor3_alm_t, T>::operator-=;
  using vectorspace_mixin<tensor3_alm_t, T>::operator*=;
  using vectorspace_mixin<tensor3_alm_t, T>::operator/=;
};

template <typename T>
tensor3_alm_t<std::complex<T> > expand(const tensor3_aij_t<T> &taij) {
  using namespace std::literals::complex_literals;

  const geom_t &geom = taij.geom;

  tensor3_aij_t<std::complex<T> > tbij(geom);
  for (int j = 0; j < geom.nphi; ++j) {
#pragma omp simd
    for (int i = 0; i < geom.ntheta; ++i) {
      // b000 =      m_x      m_y      m_z a^xyz
      // b001 =      m_x \bar m_y      m_z a^xyz
      // b010 = \bar m_x      m_y      m_z a^xyz
      // b011 = \bar m_x \bar m_y      m_z a^xyz
      // b100 =      m_x      m_y \bar m_z a^xyz
      // b101 =      m_x \bar m_y \bar m_z a^xyz
      // b110 = \bar m_x      m_y \bar m_z a^xyz
      // b111 = \bar m_x \bar m_y \bar m_z a^xyz

      // m[p] = {m, \bar m}
      const std::array<std::array<std::complex<T>, 2>, 2> m{
          {{{1, +1i}}, {{1, -1i}}}};
      for (int r = 0; r < 2; ++r) {
        for (int q = 0; q < 2; ++q) {
          for (int p = 0; p < 2; ++p) {
            std::complex<T> s = 0;
            for (int z = 0; z < 2; ++z)
              for (int y = 0; y < 2; ++y)
                for (int x = 0; x < 2; ++x)
                  s += m[p][x] * m[q][y] * m[r][z] * taij(x, y, z)(i, j);
            tbij(p, q, r)(i, j) = s;
          }
        }
      }
    }
  }

  tensor3_alm_t<std::complex<T> > talm(geom);
  talm(0, 0, 0) = expand(tbij(0, 0, 0), +3);
  talm(0, 0, 1) = expand(tbij(0, 0, 1), +1);
  talm(0, 1, 0) = expand(tbij(0, 1, 0), +1);
  talm(0, 1, 1) = expand(tbij(0, 1, 1), -1);
  talm(1, 0, 0) = expand(tbij(1, 0, 0), +1);
  talm(1, 0, 1) = expand(tbij(1, 0, 1), -1);
  talm(1, 1, 0) = expand(tbij(1, 1, 0), -1);
  talm(1, 1, 1) = expand(tbij(1, 1, 1), -3);

  return talm;
}

template <typename T>
tensor3_aij_t<T> evaluate(const tensor3_alm_t<std::complex<T> > &talm) {
  using namespace std::literals::complex_literals;

  const geom_t &geom = talm.geom;

  tensor3_aij_t<std::complex<T> > tbij(geom);
  tbij(0, 0, 0) = evaluate(talm(0, 0, 0), +3);
  tbij(0, 0, 1) = evaluate(talm(0, 0, 1), +1);
  tbij(0, 1, 0) = evaluate(talm(0, 1, 0), +1);
  tbij(0, 1, 1) = evaluate(talm(0, 1, 1), -1);
  tbij(1, 0, 0) = evaluate(talm(1, 0, 0), +1);
  tbij(1, 0, 1) = evaluate(talm(1, 0, 1), -1);
  tbij(1, 1, 0) = evaluate(talm(1, 1, 0), -1);
  tbij(1, 1, 1) = evaluate(talm(1, 1, 1), -3);

  tensor3_aij_t<T> taij(geom);
  for (int j = 0; j < geom.nphi; ++j) {
#pragma omp simd
    for (int i = 0; i < geom.ntheta; ++i) {
      // b000 =      m_x      m_y      m_z a^xyz
      // b001 =      m_x \bar m_y      m_z a^xyz
      // b010 = \bar m_x      m_y      m_z a^xyz
      // b011 = \bar m_x \bar m_y      m_z a^xyz
      // b100 =      m_x      m_y \bar m_z a^xyz
      // b101 =      m_x \bar m_y \bar m_z a^xyz
      // b110 = \bar m_x      m_y \bar m_z a^xyz
      // b111 = \bar m_x \bar m_y \bar m_z a^xyz

      // m[p] = {m, \bar m}
      const std::array<std::array<std::complex<T>, 2>, 2> minv2{
          {{{+1, +1}}, {{-1i, +1i}}}}; // minv2[x][p] = 2 * inv(m)
      for (int z = 0; z < 2; ++z) {
        for (int y = 0; y < 2; ++y) {
          for (int x = 0; x < 2; ++x) {
            std::complex<T> s = 0;
            for (int r = 0; r < 2; ++r)
              for (int q = 0; q < 2; ++q)
                for (int p = 0; p < 2; ++p)
                  s += minv2[x][p] * minv2[y][q] * minv2[z][r] *
                       tbij(p, q, r)(i, j);
            taij(x, y, z)(i, j) = real(s) / 8;
          }
        }
      }
    }
  }

  return taij;
}

////////////////////////////////////////////////////////////////////////////////

#if 0

template <typename T>
alm_t<std::complex<T> > filter(const alm_t<std::complex<T> > &alm,
                               const int lmax) {
  const geom_t &geom = alm.geom;
  alm_t<std::complex<T> > blm(geom, alm.spin);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = l <= lmax ? alm(l, m) : 0;
  return blm;
}

template <typename T>
alm_t<std::complex<T> > grad(const alm_t<std::complex<T> > &alm) {
  const geom_t &geom = alm.geom;
  alm_t<std::complex<T> > blm(geom, 1);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = -sqrt(T(l * (l + 1))) * alm(l, m);
  return blm;
}

template <typename T>
alm_t<std::complex<T> > div(const alm_t<std::complex<T> > &alm) {
  const geom_t &geom = alm.geom;
  alm_t<std::complex<T> > blm(geom, 0);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = sqrt(T(l * (l + 1))) * alm(l, m);
  return blm;
}

// laplace = div grad
template <typename T>
alm_t<std::complex<T> > laplace(const alm_t<std::complex<T> > &alm) {
  const geom_t &geom = alm.geom;
  alm_t<std::complex<T> > blm(geom, 0);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = -T(l * (l + 1)) * alm(l, m);
  return blm;
}

template <typename T>
alm_t<std::complex<T> > expand_grad(const aij_t<std::complex<T> > &aij) {
  const geom_t &geom = aij.geom;
  aij_t<std::complex<T> > aij_p(geom);
  // aij_t<std::complex<T> > aij_m(geom);
  for (int i = 0; i < geom.ntheta; ++i) {
    for (int j = 0; j < geom.nphi; ++j) {
      aij_p(i, j) = aij(i, j);
      // aij_m(i, j) = conj(aij(i, j));
    }
  }
  const alm_t<std::complex<T> > alm_p = expand(aij_p, +1);
  // const alm_t<std::complex<T>> alm_m = expand(aij_m, -1);
  alm_t<std::complex<T> > alm(geom, 0);
  for (int l = 0; l <= geom.lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      // alm(l, m) = 1 / T(2) * (alm_p(l, m) - alm_m(l, m));
      // const std::complex<T> blm = -1i / T(2) * (alm_p(l, m) + alm_m(l, m));
      // assert(abs(blm) <= 1.0e-12);
      alm(l, m) = alm_p(l, m);
    }
  }
  return alm;
}

template <typename T>
aij_t<std::complex<T> > evaluate_grad(const alm_t<std::complex<T> > &alm) {
  const geom_t &geom = alm.geom;
  alm_t<std::complex<T> > alm_p(geom, 0);
  // alm_t<std::complex<T>> alm_m(geom);
  for (int l = 0; l <= geom.lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      alm_p(l, m) = alm(l, m);
      // alm_m(l, m) = -alm(l, m);
    }
  }
  const aij_t<std::complex<T> > aij_p = evaluate(alm_p, +1);
  // const aij_t<std::complex<T> > aij_m = evaluate(alm_m, -1);
  aij_t<std::complex<T> > aij(geom);
  for (int i = 0; i < geom.ntheta; ++i) {
    for (int j = 0; j < geom.nphi; ++j) {
      // aij(i, j) = aij_p(i, j) + conj(aij_m(i, j));
      aij(i, j) = aij_p(i, j);
    }
  }
  return aij;
}

#endif

} // namespace AHFinder

#endif // #ifndef DISCRETIZATION_HXX
