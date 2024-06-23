#ifndef MULTIPOLE_MULTIPOLE_HXX
#define MULTIPOLE_MULTIPOLE_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

#include <vector>

namespace Multipole {

// information about variable which we decompose
struct variable_desc {
  int index;
  int imag_index;
  CCTK_INT spin_weight;
  std::string name;
};

// a simple array class to hold complex modes abs(m) <= l, l <= lmax, for
// nradii radii for nvars variables
class mode_array {
public:
  mode_array(int nvars_, int nradii_, int lmax_)
      : nvars(nvars_), nradii(nradii_), lmax(lmax_),
        modes(size_t(nvars * nradii * (lmax + 1) * (lmax + 1) * 2)) {}
  virtual ~mode_array() {}
  // default copy and assignment is ok

  CCTK_REAL &operator()(int v, int ri, int l, int m, bool is_im) {
    return modes.at(mode_idx(v, ri, l, m, is_im));
  }

  const CCTK_REAL &operator()(int v, int ri, int l, int m, bool is_im) const {
    return modes.at(mode_idx(v, ri, l, m, is_im));
  }

  int get_nvars() const { return nvars; }
  int get_nradii() const { return nradii; }
  int get_lmax() const { return lmax; }

private:
  size_t mode_idx(int v, int ri, int l, int m, int is_im) const {
    assert(v >= 0 && v < nvars);
    assert(ri >= 0 && ri < nradii);
    assert(l >= 0 && l <= lmax);
    assert(m <= l && -m <= l);
    return size_t(v * nradii * (lmax + 1) * (lmax + 1) * 2 +
                  ri * (lmax + 1) * (lmax + 1) * 2 + (l * l + (m + l)) * 2 +
                  is_im);
  }

  const int nvars;
  const int nradii;
  const int lmax;
  std::vector<CCTK_REAL> modes;
};

} // namespace Multipole

#endif // #ifndef MULTIPOLE_MULTIPOLE_HXX
