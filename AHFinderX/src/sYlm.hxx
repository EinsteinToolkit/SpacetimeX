#ifndef SYLM_HXX
#define SYLM_HXX

// Implement a few spin-weighted spherical harmonics for testing

// #include <iostream>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <type_traits>

namespace AHFinder {
using namespace std;
using namespace std::complex_literals;

// These spherical harmonics are implemented as defined in the
// "Spin-Weighted Spherical Harmonics Calculator" at
// <https://www.black-holes.org/code/SpinWeightedSphericalHarmonics.nb>,
// which is referenced by Wikipedia on
// <https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics>.

constexpr int sYlm_smin = -1;
constexpr int sYlm_smax = 1;
constexpr int sYlm_lmax = 2;

constexpr int bitsign(const int i) { return i % 2 == 0 ? 1 : -1; }

// s = -1

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == -1 && l == 1 && m == -1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), -phi) * sqrt(3 / M_PI) *
         pow(sin(theta / 2), 2);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == -1 && l == 1 && m == 0), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * sqrt(3 / (2 * M_PI)) * sin(theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == -1 && l == 1 && m == 1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), phi) * sqrt(3 / M_PI) * pow(cos(theta / 2), 2);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == -1 && l == 2 && m == -2), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), -2 * phi) * sqrt(5 / M_PI) *
         pow(sin(theta / 2), 2) * sin(theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == -1 && l == 2 && m == -1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), -phi) * sqrt(5 / M_PI) * (1 + 2 * cos(theta)) *
         pow(sin(theta / 2), 2);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == -1 && l == 2 && m == 0), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * sqrt(15 / (2 * M_PI)) * cos(theta) * sin(theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == -1 && l == 2 && m == 1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), phi) * sqrt(5 / M_PI) *
         pow(cos(theta / 2), 2) * (-1 + 2 * cos(theta));
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == -1 && l == 2 && m == 2), complex<T> >
sYlm(const T theta, const T phi) {
  return polar(T(1), 2 * phi) * sqrt(5 / M_PI) * pow(cos(theta / 2), 3) *
         sin(theta / 2);
}

// s = 0

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 0 && m == 0), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / (2 * sqrt(M_PI));
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 1 && m == -1), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / T(2) * polar(T(1), -phi) * sqrt(3 / (2 * M_PI)) * sin(theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 1 && m == 0), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / T(2) * sqrt(3 / M_PI) * cos(theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 1 && m == 1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), phi) * sqrt(3 / (2 * M_PI)) * sin(theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == -2), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / T(4) * polar(T(1), -2 * phi) * sqrt(15 / (2 * M_PI)) *
         pow(sin(theta), 2);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == -1), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / T(4) * polar(T(1), -phi) * sqrt(15 / (2 * M_PI)) * sin(2 * theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == 0), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / T(8) * sqrt(5 / M_PI) * (1 + 3 * cos(2 * theta));
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == 1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), phi) * sqrt(15 / (2 * M_PI)) * cos(theta) *
         sin(theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == 2), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / T(4) * polar(T(1), 2 * phi) * sqrt(15 / (2 * M_PI)) *
         pow(sin(theta), 2);
}

// s = 1

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 1 && l == 1 && m == -1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), -phi) * sqrt(3 / M_PI) *
         pow(cos(theta / 2), 2);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 1 && l == 1 && m == 0), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / T(2) * sqrt(3 / (2 * M_PI)) * sin(theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 1 && l == 1 && m == 1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), phi) * sqrt(3 / M_PI) * pow(sin(theta / 2), 2);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 1 && l == 2 && m == -2), complex<T> >
sYlm(const T theta, const T phi) {
  return -polar(T(1), -2 * phi) * sqrt(5 / M_PI) * pow(cos(theta / 2), 3) *
         sin(theta / 2);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 1 && l == 2 && m == -1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), -phi) * sqrt(5 / M_PI) *
         pow(cos(theta / 2), 2) * (-1 + 2 * cos(theta));
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 1 && l == 2 && m == 0), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / T(2) * sqrt(15 / (2 * M_PI)) * cos(theta) * sin(theta);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 1 && l == 2 && m == 1), complex<T> >
sYlm(const T theta, const T phi) {
  return -1 / T(2) * polar(T(1), phi) * sqrt(5 / M_PI) * (1 + 2 * cos(theta)) *
         pow(sin(theta / 2), 2);
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 1 && l == 2 && m == 2), complex<T> >
sYlm(const T theta, const T phi) {
  return 1 / T(2) * polar(T(1), 2 * phi) * sqrt(5 / M_PI) *
         pow(sin(theta / 2), 2) * sin(theta);
}

// partial derivatives of s = 0: [\partial_theta, 1/\sin\theta \partial_\phi]

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 0 && m == 0), array<complex<T>, 2> >
dsYlm(const T theta, const T phi) {
  return {0, 0};
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 1 && m == -1), array<complex<T>, 2> >
dsYlm(const T theta, const T phi) {
  return {1 / T(2) * polar(T(1), -phi) * sqrt(3 / (2 * M_PI)) * cos(theta),
          -1.0i / T(2) * polar(T(1), -phi) * sqrt(3 / (2 * M_PI))};
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 1 && m == 0), array<complex<T>, 2> >
dsYlm(const T theta, const T phi) {
  return {-1 / T(2) * sqrt(3 / M_PI) * sin(theta), 0};
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 1 && m == 1), array<complex<T>, 2> >
dsYlm(const T theta, const T phi) {
  return {-1 / T(2) * polar(T(1), phi) * sqrt(3 / (2 * M_PI)) * cos(theta),
          -1.0i / T(2) * polar(T(1), phi) * sqrt(3 / (2 * M_PI))};
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == -2), array<complex<T>, 2> >
dsYlm(const T theta, const T phi) {
  return {1 / T(2) * polar(T(1), -2 * phi) * sqrt(15 / (2 * M_PI)) *
              cos(theta) * sin(theta),
          -1.0i / T(2) * polar(T(1), -2 * phi) * sqrt(15 / (2 * M_PI)) *
              sin(theta)};
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == -1), array<complex<T>, 2> >
dsYlm(const T theta, const T phi) {
  return {1 / T(2) * polar(T(1), -phi) * sqrt(15 / (2 * M_PI)) * cos(2 * theta),
          -1.0i / T(2) * polar(T(1), -phi) * sqrt(15 / (2 * M_PI)) *
              cos(theta)};
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == 0), array<complex<T>, 2> >
dsYlm(const T theta, const T phi) {
  return {-3 / T(4) * sqrt(5 / M_PI) * sin(2 * theta), 0};
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == 1), array<complex<T>, 2> >
dsYlm(const T theta, const T phi) {
  return {-1 / T(2) * polar(T(1), phi) * sqrt(15 / (2 * M_PI)) * cos(2 * theta),
          -1.0i / T(2) * polar(T(1), phi) * sqrt(15 / (2 * M_PI)) * cos(theta)};
}

template <int s, int l, int m, typename T>
constexpr enable_if_t<(s == 0 && l == 2 && m == 2), array<complex<T>, 2> >
dsYlm(const T theta, const T phi) {
  return {1 / T(2) * polar(T(1), 2 * phi) * sqrt(15 / (2 * M_PI)) * cos(theta) *
              sin(theta),
          1.0i / T(2) * polar(T(1), 2 * phi) * sqrt(15 / (2 * M_PI)) *
              sin(theta)};
}

//

template <typename T>
constexpr complex<T> sYlm(const int s, const int l, const int m, const T theta,
                          const T phi) {
  if (s == -1 && l == 1 && m == -1)
    return sYlm<-1, 1, -1>(theta, phi);
  if (s == -1 && l == 1 && m == 0)
    return sYlm<-1, 1, 0>(theta, phi);
  if (s == -1 && l == 1 && m == 1)
    return sYlm<-1, 1, 1>(theta, phi);
  if (s == -1 && l == 2 && m == -2)
    return sYlm<-1, 2, -2>(theta, phi);
  if (s == -1 && l == 2 && m == -1)
    return sYlm<-1, 2, -1>(theta, phi);
  if (s == -1 && l == 2 && m == 0)
    return sYlm<-1, 2, 0>(theta, phi);
  if (s == -1 && l == 2 && m == 1)
    return sYlm<-1, 2, 1>(theta, phi);
  if (s == -1 && l == 2 && m == 2)
    return sYlm<-1, 2, 2>(theta, phi);

  if (s == 0 && l == 0 && m == 0)
    return sYlm<0, 0, 0>(theta, phi);
  if (s == 0 && l == 1 && m == -1)
    return sYlm<0, 1, -1>(theta, phi);
  if (s == 0 && l == 1 && m == 0)
    return sYlm<0, 1, 0>(theta, phi);
  if (s == 0 && l == 1 && m == 1)
    return sYlm<0, 1, 1>(theta, phi);
  if (s == 0 && l == 2 && m == -2)
    return sYlm<0, 2, -2>(theta, phi);
  if (s == 0 && l == 2 && m == -1)
    return sYlm<0, 2, -1>(theta, phi);
  if (s == 0 && l == 2 && m == 0)
    return sYlm<0, 2, 0>(theta, phi);
  if (s == 0 && l == 2 && m == 1)
    return sYlm<0, 2, 1>(theta, phi);
  if (s == 0 && l == 2 && m == 2)
    return sYlm<0, 2, 2>(theta, phi);

  if (s == 1 && l == 1 && m == -1)
    return sYlm<1, 1, -1>(theta, phi);
  if (s == 1 && l == 1 && m == 0)
    return sYlm<1, 1, 0>(theta, phi);
  if (s == 1 && l == 1 && m == 1)
    return sYlm<1, 1, 1>(theta, phi);
  if (s == 1 && l == 2 && m == -2)
    return sYlm<1, 2, -2>(theta, phi);
  if (s == 1 && l == 2 && m == -1)
    return sYlm<1, 2, -1>(theta, phi);
  if (s == 1 && l == 2 && m == 0)
    return sYlm<1, 2, 0>(theta, phi);
  if (s == 1 && l == 2 && m == 1)
    return sYlm<1, 2, 1>(theta, phi);
  if (s == 1 && l == 2 && m == 2)
    return sYlm<1, 2, 2>(theta, phi);

  abort();
}

template <typename T>
constexpr array<complex<T>, 2> dsYlm(const int s, const int l, const int m,
                                     const T theta, const T phi) {
  if (s == 0 && l == 0 && m == 0)
    return dsYlm<0, 0, 0>(theta, phi);
  if (s == 0 && l == 1 && m == -1)
    return dsYlm<0, 1, -1>(theta, phi);
  if (s == 0 && l == 1 && m == 0)
    return dsYlm<0, 1, 0>(theta, phi);
  if (s == 0 && l == 1 && m == 1)
    return dsYlm<0, 1, 1>(theta, phi);
  if (s == 0 && l == 2 && m == -2)
    return dsYlm<0, 2, -2>(theta, phi);
  if (s == 0 && l == 2 && m == -1)
    return dsYlm<0, 2, -1>(theta, phi);
  if (s == 0 && l == 2 && m == 0)
    return dsYlm<0, 2, 0>(theta, phi);
  if (s == 0 && l == 2 && m == 1)
    return dsYlm<0, 2, 1>(theta, phi);
  if (s == 0 && l == 2 && m == 2)
    return dsYlm<0, 2, 2>(theta, phi);

  abort();
}

} // namespace AHFinder

#endif //#ifndef SYLM_HXX
