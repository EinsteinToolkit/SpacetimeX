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
  Surface(int nTheta, int nPhi, bool isMidPoint);
  virtual ~Surface() {} // Ready for potential inheritance

  std::vector<CCTK_REAL> getTheta() const { return theta_; }
  std::vector<CCTK_REAL> getPhi() const { return phi_; }

private:
  friend class Sphere; // allow Sphere to access private members

  const int nTheta_, nPhi_;
  std::vector<CCTK_REAL> theta_, phi_;
  std::vector<CCTK_REAL> x_, y_, z_; // embedding information

  inline int index2D(int it, int ip) { return it + (nTheta_ + 1) * ip; }
};

// 2D Sphere
class Sphere : public Surface {
public:
  Sphere(int nTheta, int nPhi, bool isMidPoint);
  ~Sphere() {} // Ready for potential inheritance

  void setRadius(CCTK_REAL radius);

private:
  std::vector<CCTK_REAL> xhat_, yhat_, zhat_; // unit sphere
};

} // namespace Multipole

#endif // #ifndef MULTIPOLE_SURFACE_HXX
