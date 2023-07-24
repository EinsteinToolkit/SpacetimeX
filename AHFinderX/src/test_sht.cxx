#include "sYlm.hxx"

#include <ssht/ssht.h>

#include <cctk.h>
#include <cctk_Arguments.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

namespace AHFinder {
using namespace std;
using namespace std::complex_literals;

extern "C" void AHFinder_test_sht(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_test_sht;

  // Choose resolution

  const int nmodes = 5;

  // Access grid and spectral information

  const int ntheta = nmodes;
  const int nphi = 2 * nmodes - 1;
  const int npoints = ntheta * nphi;
  const auto coord_theta = [&](const int i, const int j) -> CCTK_REAL {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    return ssht_sampling_mw_t2theta(i, nmodes);
  };
  const auto coord_phi = [&](const int i, const int j) -> CCTK_REAL {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    return ssht_sampling_mw_p2phi(j, nmodes);
  };
  const auto gind = [&](const int i, const int j) -> int {
    // 0 <= i < ntheta
    // 0 <= j < nphi
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const int ind = j + nphi * i;
    assert(ind >= 0 && ind < npoints);
    return ind;
  };

  const int lmax = nmodes - 1;
  const int ncoeffs = nmodes * nmodes;
  const auto cind = [&](const int l, const int m) -> int {
    // 0 <= l <= lmax
    // -l <= m <= l
    assert(l >= 0 && l <= lmax);
    assert(m >= -l && m <= l);
    int ind;
    ssht_sampling_elm2ind(&ind, l, m);
    assert(ind >= 0 && ind < ncoeffs);
    return ind;
  };

  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]

  // TODO: Test s_Y_lm* = (-1)^(s+m) -s_Y_l,-m

  // Test complex spherical harmonic functions

  for (int spin = sYlm_smin; spin <= sYlm_smax; ++spin) {
    for (int ll = abs(spin); ll <= sYlm_lmax; ++ll) {
      for (int mm = -ll; mm <= ll; ++mm) {

        // Choose function
        vector<complex<double> > alm(ncoeffs, NAN);
        for (int l = abs(spin); l <= lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            const complex<double> a = l == ll && m == mm ? 1 : 0;
            alm.at(cind(l, m)) = a;
          }
        }

        // Evaluate
        vector<complex<double> > aij(npoints, NAN);
        ssht_core_mw_inverse_sov_sym(aij.data(), alm.data(), nmodes, spin,
                                     method, verbosity);

        for (int i = 0; i < ntheta; ++i) {
          for (int j = 0; j < nphi; ++j) {
            const double theta = coord_theta(i, j);
            const double phi = coord_phi(i, j);
            const complex<double> s = sYlm(spin, ll, mm, theta, phi);
            const complex<double> sc = double(bitsign(spin + mm)) *
                                       conj(sYlm(-spin, ll, -mm, theta, phi));
            const complex<double> a = aij.at(gind(i, j));
            if (!(abs(a - s) <= 1.0e-12))
              CCTK_VINFO("spin=%d,ll=%d,mm=% i=%d,j=%d,theta=%f,phi=%f "
                         "s=(%.17g,%.17g) a=(%.17g,%.17g)",
                         spin, ll, mm, i, j, theta, phi, real(s), imag(s),
                         real(a), imag(a));
            assert(abs(a - s) <= 1.0e-12);
            assert(abs(s - sc) <= 1.0e-12);
          }
        }

        // Expand
        ssht_core_mw_forward_sov_conv_sym(alm.data(), aij.data(), nmodes, spin,
                                          method, verbosity);

        for (int l = abs(spin); l <= lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            const complex<double> a = l == ll && m == mm ? 1 : 0;
            assert(abs(alm.at(cind(l, m)) - a) <= 1.0e-12);
          }
        }

#if 0
        // Check phase convention
        vector<complex<double> > blm(npoints, NAN);
        ssht_core_mw_forward_sov_conv_sym(blm.data(), aij.data(), nmodes, -spin,
                                          method, verbosity);

        for (int l = abs(spin); l <= lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            // const complex<double> b =
            //     double(bitsign(spin + m)) * conj(alm.at(cind(l, -m)));
            const complex<double> b =
                double(bitsign(spin + 1)) * conj(alm.at(cind(l, -m)));
            const complex<double> a = l == ll && m == -mm ? 1 : 0;
            if (!(abs(b - a) <= 1.0e-12))
              CCTK_VINFO("spin=%d,ll=%d,mm=%d l=%d,m=%d "
                         "b=(%.17g,%.17g) a=(%.17g,%.17g)",
                         spin, ll, mm, l, m, real(b), imag(b), real(a),
                         imag(a));
            assert(abs(b - a) <= 1.0e-12);
          }
        }
#endif
      }
    }
  }

  // Test real spherical harmonic functions

  {
    const int spin = 0;
    for (int ll = abs(spin); ll <= sYlm_lmax; ++ll) {
      for (int mm = -ll; mm <= ll; ++mm) {

        vector<complex<double> > alm(ncoeffs, complex<double>{NAN, NAN});
        for (int l = abs(spin); l <= lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            // f*_lm = (-1)^m f_l,-m
            const complex<double> c = mm < 0 ? 1.0i : 1.0;
            complex<double> a{0};
            if (l == ll && m == mm)
              a += c;
            if (l == ll && -m == mm)
              a += double(bitsign(mm)) * conj(c);
            alm.at(cind(l, m)) = a;
          }
        }

        vector<double> aij(npoints, NAN);
        ssht_core_mw_inverse_sov_sym_real(aij.data(), alm.data(), nmodes,
                                          method, verbosity);

        for (int i = 0; i < ntheta; ++i) {
          for (int j = 0; j < nphi; ++j) {
            const double theta = coord_theta(i, j);
            const double phi = coord_phi(i, j);
            const complex<double> c = mm < 0 ? 1.0i : 1.0;
            const complex<double> scp = sYlm(spin, ll, mm, theta, phi);
            const complex<double> scm = sYlm(spin, ll, -mm, theta, phi);
            const complex<double> sc =
                c * scp + double(bitsign(mm)) * conj(c) * scm;
            assert(abs(imag(sc)) <= 1.0e-12);
            const double s = real(sc);
            const double a = aij.at(gind(i, j));
            if (!(abs(a - s) <= 1.0e-12))
              CCTK_VINFO("spin=%d,ll=%d,mm=%d i=%d,j=%d,theta=%f,phi=%f "
                         "s=%.17g a=%.17g",
                         spin, ll, mm, i, j, theta, phi, s, a);
            assert(abs(a - s) <= 1.0e-12);
          }
        }

        fill(alm.begin(), alm.end(), complex<double>{NAN, NAN});
        ssht_core_mw_forward_sov_conv_sym_real(alm.data(), aij.data(), nmodes,
                                               method, verbosity);

        for (int l = abs(spin); l <= lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            const complex<double> c = mm < 0 ? 1.0i : 1.0;
            complex<double> a{0};
            if (l == ll && m == mm)
              a += c;
            if (l == ll && -m == mm)
              a += double(bitsign(mm)) * conj(c);
            assert(abs(alm.at(cind(l, m)) - a) <= 1.0e-12);
          }
        }
      }
    }
  }

  // Test gradient of real spherical harmonic functions

  for (const int spin : {-1, 1}) {
    for (int ll = abs(spin); ll <= sYlm_lmax; ++ll) {
      for (int mm = -ll; mm <= ll; ++mm) {

        for (int i = 0; i < ntheta; ++i) {
          for (int j = 0; j < nphi; ++j) {
            const double theta = coord_theta(i, j);
            const double phi = coord_phi(i, j);
            // ds =[\partial_theta, 1/\sin\theta \partial_\phi]
            const array<complex<double>, 2> ds = dsYlm(0, ll, mm, theta, phi);
            const complex<double> sc = -double(spin) /
                                       sqrt(double(ll * (ll + 1))) *
                                       (ds[0] + double(spin) * 1.0i * ds[1]);
            // const complex<double> ac = aij.at(gind(i, j));
            const complex<double> ac = sYlm(spin, ll, mm, theta, phi);
            if (!(abs(ac - sc) <= 1.0e-12))
              CCTK_VINFO("spin=%d,ll=%d,mm=%d i=%d,j=%d,theta=%f,phi=%f "
                         "sc=(%.17g,%.17g) ac=(%.17g,%.17g)",
                         spin, ll, mm, i, j, theta, phi, real(sc), imag(sc),
                         real(ac), imag(ac));
            assert(abs(ac - sc) <= 1.0e-12);
          }
        }
      }
    }
  }
}

} // namespace AHFinder
