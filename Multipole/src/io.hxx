#ifndef MULTIPOLE_IO_HXX
#define MULTIPOLE_IO_HXX

#include "multipole.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>

namespace Multipole {
using namespace std;

FILE *OpenOutputFile(CCTK_ARGUMENTS, const std::string &name);

void OutputComplexToFile(CCTK_ARGUMENTS, const string &name, CCTK_REAL redata,
                         CCTK_REAL imdata);

void OutputComplexToH5File(CCTK_ARGUMENTS, const VariableParse vars[],
                           const CCTK_REAL radii[], const ModeArray &modes);

} // namespace Multipole

#endif // #ifndef MULTIPOLE_IO_HXX
