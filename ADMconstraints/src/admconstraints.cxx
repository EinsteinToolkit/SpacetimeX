#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace ADMconstraints {
using namespace Arith;
using namespace Loop;

template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T Power(T x, int y) {
  return (y == 2) ? Arith::pow2(x) : Arith::pown(x, y);
}

extern "C" void ADMconstraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMconstraints;
  DECLARE_CCTK_PARAMETERS;
}

} // namespace ADMconstraints
