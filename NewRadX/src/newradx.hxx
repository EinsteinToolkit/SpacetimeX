#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cmath>
#include "loop.hxx"
#include "loop_device.hxx"

namespace NewRadX {

// Adapted from BSSN_MoL's files NewRad.F and newrad.h. This code was probably
// originally written by Miguel Alcubierre.
void NewRadX_Apply(const cGH *restrict const cctkGH,
            const Loop::GF3D2<const CCTK_REAL> var, Loop::GF3D2<CCTK_REAL> rhs,
            const CCTK_REAL var0, //!< value at infinity
            const CCTK_REAL v0,    //!< propagation speed
      const CCTK_REAL radpower //!< radial fall-off exponent
);
} // namespace NewRadX
