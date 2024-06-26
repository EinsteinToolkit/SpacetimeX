/* surface.hxx */
/* (c) Liwei Ji 06/2024 */

#ifndef MULTIPOLE_SURFACE_HXX
#define MULTIPOLE_SURFACE_HXX

#include <cctk.h>

#include <vector>

namespace Multipole {

class Surface {
public:
  Surface(int nTheta, int nPhi, bool isMidPoint);
  virtual ~Surface() {} // Ready for potential inheritance

  std::vector<CCTK_REAL> getTheta() const { return theta_; }
  std::vector<CCTK_REAL> getPhi() const { return phi_; }

private:
  const int nTheta_, nPhi_;
  std::vector<CCTK_REAL> theta_, phi_;
  std::vector<CCTK_REAL> x_, y_, z_;

  inline int index2D(int it, int ip) { return it + (nTheta_ + 1) * ip; }
};

} // namespace Multipole

#endif // #ifndef MULTIPOLE_SURFACE_HXX
