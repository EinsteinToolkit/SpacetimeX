#ifndef MULTIPOLE_IO_HXX
#define MULTIPOLE_IO_HXX

#include "multipole.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>

namespace Multipole {
using namespace std;

enum mp_coord { mp_theta, mp_phi };

void Output1D(CCTK_ARGUMENTS, const string &name, CCTK_REAL const th[],
              CCTK_REAL const ph[], mp_coord coord, CCTK_REAL const data[]);

void OutputComplexToFile(CCTK_ARGUMENTS, const string &name, CCTK_REAL redata,
                         CCTK_REAL imdata);

void OutputComplexToH5File(CCTK_ARGUMENTS, const VariableParse vars[],
                           const CCTK_REAL radii[], const ModeArray &modes);

} // namespace Multipole

#endif // #ifndef MULTIPOLE_IO_HXX
