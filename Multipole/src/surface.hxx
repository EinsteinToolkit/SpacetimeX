/* surface.hxx */
/* (c) Liwei Ji 06/2024 */

#ifndef MULTIPOLE_SURFACE_HXX
#define MULTIPOLE_SURFACE_HXX

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

protected:
  inline int index2D(int it, int ip) { return it + (nTheta_ + 1) * ip; }

  const int nTheta_, nPhi_;
  std::vector<CCTK_REAL> theta_, phi_;
  std::vector<CCTK_REAL> x_, y_, z_; // embedding map
};

// 2D Sphere
class Sphere : public Surface {
public:
  Sphere(int nTheta, int nPhi, bool isMidPoint)
      : Surface(nTheta, nPhi, isMidPoint) {
    xhat_.resize(theta_.size());
    yhat_.resize(theta_.size());
    zhat_.resize(theta_.size());
    for (size_t i = 0; i < theta_.size(); ++i) {
      xhat_[i] = std::cos(phi_[i]) * std::sin(theta_[i]);
      yhat_[i] = std::sin(phi_[i]) * std::sin(theta_[i]);
      zhat_[i] = std::cos(theta_[i]);
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
};

} // namespace Multipole

#endif // #ifndef MULTIPOLE_SURFACE_HXX
