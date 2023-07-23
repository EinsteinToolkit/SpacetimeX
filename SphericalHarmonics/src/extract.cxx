#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <hdf5.h>

#include <ssht/ssht.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <string>
#include <vector>

namespace SphericalHarmonics {

extern "C" void SphericalHarmonics_extract(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_SphericalHarmonics_extract;
  DECLARE_CCTK_PARAMETERS;

  if (!(out_every > 0 && cctk_iteration % out_every == 0))
    return;

  const int psi4re_ind = CCTK_VarIndex("Weyl::Psi4re");
  assert(psi4re_ind >= 0);
  const int psi4im_ind = CCTK_VarIndex("Weyl::Psi4im");
  assert(psi4im_ind >= 0);

  const int nvars = 2;
  const std::array<CCTK_INT, nvars> varinds{psi4re_ind, psi4im_ind};

  const std::array<CCTK_INT, nvars> operations{0, 0}; // no derivatives

  // Coefficients
  const int lmax = 1000; // choice
  const int nmodes = lmax + 1;
  const int ncoeffs = nmodes * nmodes;

  // Grid
  const int ntheta = nmodes;
  const int nphi = 2 * nmodes - 1;
  const int npoints = ntheta * nphi;

  const auto coord_theta = [&](const int i, const int j) {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const double theta = ssht_sampling_mw_t2theta(i, nmodes);
    assert(theta >= 0 && theta <= M_PI);
    return theta;
  };

  const auto coord_phi = [&](const int i, const int j) {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const double phi = ssht_sampling_mw_p2phi(j, nmodes);
    assert(phi >= 0 && phi < 2 * M_PI);
    return phi;
  };

  const auto coord_dtheta = [&](const int i, const int j) {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const double theta_m =
        i == 0 ? 0.0 : (coord_theta(i - 1, j) + coord_theta(i, j)) / 2;
    const double theta_p =
        i == ntheta - 1 ? M_PI
                        : (coord_theta(i, j) + coord_theta(i + 1, j)) / 2;
    return theta_p - theta_m;
  };

  const auto coord_dphi = [&](const int i, const int j) {
    return 2 * M_PI / nphi;
  };

  const auto gind = [&](const int i, const int j) {
    // 0 <= i < ntheta
    // 0 <= j < nphi
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const int ind = j + nphi * i;
    assert(ind >= 0 && ind < npoints);
    return ind;
  };

  const auto cind = [&](const int l, const int m) {
    // 0 <= l <= lmax
    // -l <= m <= l
    assert(l >= 0 && l <= lmax);
    assert(m >= -l && m <= l);
    int ind;
    ssht_sampling_elm2ind(&ind, l, m);
    assert(ind >= 0 && ind < ncoeffs);
    return ind;
  };

  // Coordinates
  const int npoints_local = CCTK_MyProc(0) == 0 ? npoints : 0;
  std::vector<CCTK_REAL> coord_x(npoints_local), coord_y(npoints_local),
      coord_z(npoints_local);
  if (CCTK_MyProc(0) == 0) {
    for (int i = 0; i < ntheta; ++i) {
      for (int j = 0; j < nphi; ++j) {
        const int gi = gind(i, j);
        const double r = 64.0; // choice
        const double theta = coord_theta(i, j);
        const double phi = coord_phi(i, j);
        const double x0 = 0.0; // choice
        const double y0 = 0.0; // choice
        const double z0 = 0.0; // choice
        const double x = x0 + r * sin(theta) * cos(phi);
        const double y = y0 + r * sin(theta) * sin(phi);
        const double z = z0 + r * cos(theta);
        coord_x.at(gi) = x;
        coord_y.at(gi) = y;
        coord_z.at(gi) = z;
      }
    }
  }

  // Data
  std::vector<CCTK_REAL> psi4re(npoints_local), psi4im(npoints_local);
  std::array<CCTK_POINTER, nvars> ptrs{psi4re.data(), psi4im.data()};

  Interpolate(cctkGH, npoints_local, coord_x.data(), coord_y.data(),
              coord_z.data(), nvars, varinds.data(), operations.data(),
              ptrs.data());

  // Calculate and output modes only on a single process
  if (CCTK_MyProc(0) == 0) {

    for (int n = 0; n < npoints; ++n) {
      if (!(isfinite(psi4re.at(n))))
        CCTK_VWARN(CCTK_WARN_ALERT, "psi4re is not finite: n=%d [x,y,z]=[%.17g,%.17g,%.17g]", n,
                   coord_x.at(n), coord_y.at(n), coord_z.at(n));
      // assert(isfinite(psi4re.at(n)));
      // assert(isfinite(psi4im.at(n)));
    }

    // Convert values to complex numbers
    std::vector<std::complex<double> > psi4(npoints);
    for (int n = 0; n < npoints; ++n)
      psi4.at(n) = std::complex<double>(psi4re.at(n), psi4im.at(n));

    // Expand into spherical harmonics
    const int spin = -2; // choice (depends on grid function)
    const ssht_dl_method_t method = SSHT_DL_RISBO;
    const int verbosity = 0; // [0..5]
    std::vector<std::complex<double> > psi4lm(ncoeffs);
    ssht_core_mw_forward_sov_conv_sym(psi4lm.data(), psi4.data(), nmodes, spin,
                                      method, verbosity);

    for (int l = abs(spin); l <= lmax; ++l) {
      for (int m = -l; m <= +l; ++m) {
        // assert(isfinite(real(psi4lm.at(cind(l, m)))));
        // assert(isfinite(imag(psi4lm.at(cind(l, m)))));
      }
    }

    // Output spherical harmonics
    const std::string path_name = out_dir; // choice
    static bool did_create_directory = false;
    if (!did_create_directory) {
      const int mode = 0755;
      const int ierr = CCTK_CreateDirectory(mode, path_name.c_str());
      assert(ierr >= 0);
      did_create_directory = true;
    }

    const std::string sep = "\t";
    const std::string eol = "\n";
    const std::string quote = "\"";
    using std::abs;
    const int lmin_out = abs(spin);
    const int lmax_out = 16;                   // choice
    const std::string file_name = "modes.tsv"; // choice
    const std::string output_name = path_name + "/" + file_name;
    static bool did_create_file = false;
    const std::ios_base::openmode mode =
        (did_create_file ? std::ios_base::app : std::ios_base::out) |
        std::ios_base::ate;
    std::ofstream file(output_name, mode);
    if (!did_create_file) {
      file << "iteration" << sep << "time";
      for (int l = lmin_out; l <= lmax_out; ++l) {
        for (int m = -l; m <= +l; ++m) {
          file << sep << quote << "real(l=" << l << ",m=" << m << ")" << quote;
          file << sep << quote << "imag(l=" << l << ",m=" << m << ")" << quote;
        }
      }
      file << eol;
      did_create_file = true;
    }

    file << std::setprecision(std::numeric_limits<double>::digits10 + 1)
         << std::scientific;
    file << cctk_iteration << sep << cctk_time;
    for (int l = lmin_out; l <= lmax_out; ++l) {
      for (int m = -l; m <= l; ++m) {
        file << sep << real(psi4lm.at(cind(l, m)));
        file << sep << imag(psi4lm.at(cind(l, m)));
      }
    }
    file << eol;

    file.close();

  } // if (CCTK_MyProc(0) == 0)
}

} // namespace SphericalHarmonics
