/* surface.cxx */
/* (c) Liwei Ji 06/2024 */

#include "surface.hxx"

namespace Multipole {

Surface::Surface(int nTheta, int nPhi, bool isMidPoint)
    : nTheta_(nTheta), nPhi_(nPhi) {
  constexpr CCTK_REAL PI = acos(-1.0);
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

Sphere::Sphere(int nTheta, int nPhi, bool isMidPoint)
    : Surface(nTheta, nPhi, isMidPoint) {
  xhat_.resize(theta_.size());
  yhat_.resize(theta_.size());
  zhat_.resize(theta_.size());
  for (size_t i = 0; i < theta_.size(); ++i) {
    xhat_[i] = cos(phi_[i]) * sin(theta_[i]);
    yhat_[i] = sin(phi_[i]) * sin(theta_[i]);
    zhat_[i] = cos(theta_[i]);
  }
}

void Sphere::setRadius(CCTK_REAL radius) {
  for (size_t i = 0; i < theta_.size(); ++i) {
    x_[i] = radius * xhat_[i];
    y_[i] = radius * yhat_[i];
    z_[i] = radius * zhat_[i];
  }
}

} // namespace Multipole
