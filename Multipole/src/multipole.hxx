#ifndef MULTIPOLE_MULTIPOLE_HXX
#define MULTIPOLE_MULTIPOLE_HXX

#include <cctk.h>

#include <algorithm> // std::find
#include <vector>

namespace Multipole {

// information about variable which we decompose
struct VariableParse {
  int realIndex;
  int imagIndex;
  CCTK_INT spinWeight;
  std::string name;
};

// a simple array class to hold complex modes abs(m) <= l, l <= lmax, for
// nradii radii for nvars variables
class ModeArray {
public:
  ModeArray(int numVars, int numRadii, int maxL)
      : numVars_(numVars), numRadii_(numRadii), maxL_(maxL),
        modes_(static_cast<std::size_t>(numVars * numRadii * (maxL + 1) *
                                        (maxL + 1) * 2)) {}
  // Ready for potential inheritance
  virtual ~ModeArray() {}

  // For writing to modes_
  CCTK_REAL &operator()(int v, int ri, int l, int m, bool isImaginary) {
    return modes_.at(modeIndex(v, ri, l, m, isImaginary));
  }

  // For reading from modes_
  const CCTK_REAL &operator()(int v, int ri, int l, int m,
                              bool isImaginary) const {
    return modes_.at(modeIndex(v, ri, l, m, isImaginary));
  }

  int getNumVars() const { return numVars_; }
  int getNumRadii() const { return numRadii_; }
  int getMaxL() const { return maxL_; }

private:
  inline std::size_t modeIndex(int v, int ri, int l, int m,
                               int isImaginary) const {
    assert(v >= 0 && v < numVars_);
    assert(ri >= 0 && ri < numRadii_);
    assert(l >= 0 && l <= maxL_);
    assert(m <= l && -m <= l);
    return static_cast<std::size_t>(v * numRadii_ * (maxL_ + 1) * (maxL_ + 1) *
                                        2 +
                                    ri * (maxL_ + 1) * (maxL_ + 1) * 2 +
                                    (l * l + (m + l)) * 2 + isImaginary);
  }

  const int numVars_;
  const int numRadii_;
  const int maxL_;
  std::vector<CCTK_REAL> modes_; // 1D array to store the modes
};

inline bool isIntInArray(int a, const std::vector<int> &array) {
  return std::find(array.begin(), array.end(), a) != array.end();
}

inline int findIntInArray(int a, const std::vector<int> &array) {
  auto it = std::find(array.begin(), array.end(), a);
  return (it != array.end()) ? std::distance(array.begin(), it) : -1;
}

} // namespace Multipole

#endif // #ifndef MULTIPOLE_MULTIPOLE_HXX
