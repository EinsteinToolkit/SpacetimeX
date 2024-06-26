/* surface.hxx */
/* (c) Liwei Ji 06/2024 */

#ifndef MULTIPOLE_SURFACE_HXX
#define MULTIPOLE_SURFACE_HXX

#include "sphericalharmonic.hxx"

#include <cctk.h>

#include <vector>

namespace Multipole {

// 2D Surface Embedding in 3D Space
class Surface {
public:
  Surface(int nTheta, int nPhi, bool isMidPoint)
      : nTheta_(nTheta), nPhi_(nPhi) {
    constexpr CCTK_REAL PI = std::acos(-1.0);
    const int arraySize = (nTheta + 1) * (nPhi + 1);
    theta_.resize(arraySize);
    phi_.resize(arraySize);
    x_.resize(arraySize);
    y_.resize(arraySize);
    z_.resize(arraySize);
    real_.resize(arraySize);
    imag_.resize(arraySize);

    // Add an offset for midpoint integration.
    const CCTK_REAL is_midpoint = static_cast<CCTK_REAL>(isMidPoint);
    const CCTK_REAL dTheta = PI / (nTheta + is_midpoint);
    const CCTK_REAL dPhi = 2 * PI / (nPhi + is_midpoint);

    for (int it = 0; it <= nTheta; ++it) {
      for (int ip = 0; ip <= nPhi; ++ip) {
        const int i = index2D(it, ip);
        theta_[i] = it * dTheta + 0.5 * dTheta * is_midpoint;
        phi_[i] = ip * dPhi + 0.5 * dPhi * is_midpoint;
      }
    }
  }

  virtual ~Surface() = default; // Use default for trivial destructors

  const std::vector<CCTK_REAL> &getTheta() const { return theta_; }
  const std::vector<CCTK_REAL> &getPhi() const { return phi_; }

  void interp(CCTK_ARGUMENTS, int real_idx, int imag_idx);

protected:
  inline int index2D(int it, int ip) const { return it + (nTheta_ + 1) * ip; }

  const int nTheta_, nPhi_;

  std::vector<CCTK_REAL> theta_, phi_;
  std::vector<CCTK_REAL> x_, y_, z_; // embedding map

  // a complex field to be integrated on the surface
  std::vector<CCTK_REAL> real_, imag_;
};

// 2D Sphere
class Sphere : public Surface {
public:
  Sphere(int nTheta, int nPhi, bool isMidPoint,
         const std::vector<int> &spinWeights, int lmax)
      : Surface(nTheta, nPhi, isMidPoint), spin_weights_(spinWeights),
        lmax_(lmax) {
    xhat_.resize(theta_.size());
    yhat_.resize(theta_.size());
    zhat_.resize(theta_.size());
    for (size_t i = 0; i < theta_.size(); ++i) {
      xhat_[i] = std::cos(phi_[i]) * std::sin(theta_[i]);
      yhat_[i] = std::sin(phi_[i]) * std::sin(theta_[i]);
      zhat_[i] = std::cos(theta_[i]);
    }

    realY_.resize(spinWeights.size());
    imagY_.resize(spinWeights.size());

    for (size_t si = 0; si < spinWeights.size(); ++si) {
      int sw = spinWeights[si];
      realY_[si].resize(lmax + 1);
      imagY_[si].resize(lmax + 1);

      for (int l = 0; l <= lmax; ++l) {
        realY_[si][l].resize(2 * l + 1);
        imagY_[si][l].resize(2 * l + 1);

        for (int m = -l; m <= l; ++m) {
          std::vector<CCTK_REAL> realY_s_l_m(theta_.size());
          std::vector<CCTK_REAL> imagY_s_l_m(theta_.size());
          HarmonicSetup(sw, l, m, theta_.size(), theta_, phi_, realY_s_l_m,
                        imagY_s_l_m);
          realY_[si][l][m + l] = realY_s_l_m;
          imagY_[si][l][m + l] = imagY_s_l_m;
        }
      }
    }
  }

  ~Sphere() override = default;

  void setRadius(CCTK_REAL radius) {
    for (size_t i = 0; i < theta_.size(); ++i) {
      x_[i] = radius * xhat_[i];
      y_[i] = radius * yhat_[i];
      z_[i] = radius * zhat_[i];
    }
  }

private:
  std::vector<CCTK_REAL> xhat_, yhat_, zhat_; // unit sphere
  std::vector<int> spin_weights_;             // spin weights
  // spin weighted spherical harmonics
  int lmax_;
  std::vector<std::vector<std::vector<std::vector<CCTK_REAL> > > > realY_;
  std::vector<std::vector<std::vector<std::vector<CCTK_REAL> > > > imagY_;
};

} // namespace Multipole

#endif // #ifndef MULTIPOLE_SURFACE_HXX
