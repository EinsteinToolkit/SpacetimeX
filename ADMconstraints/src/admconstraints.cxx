#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace ADMconstraints {
using namespace Arith;
using namespace Loop;

extern "C" void ADMconstraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMconstraints;
  DECLARE_CCTK_PARAMETERS;
}

} // namespace ADMconstraints
