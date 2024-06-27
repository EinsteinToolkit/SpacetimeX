#ifndef MULTIPOLE_IO_HXX
#define MULTIPOLE_IO_HXX

#include "multipole.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>

namespace Multipole {
using namespace std;

enum MpCoord { MpTheta, MpPhi };

void Output1D(CCTK_ARGUMENTS, const string &name, CCTK_REAL const th[],
              CCTK_REAL const ph[], MpCoord coord, CCTK_REAL const data[]);

void OutputComplexToFile(CCTK_ARGUMENTS, const string &name, CCTK_REAL redata,
                         CCTK_REAL imdata);

void OutputComplexToH5File(CCTK_ARGUMENTS, const VariableParse vars[],
                           const CCTK_REAL radii[], const ModeArray &modes);

} // namespace Multipole

#endif // #ifndef MULTIPOLE_IO_HXX
