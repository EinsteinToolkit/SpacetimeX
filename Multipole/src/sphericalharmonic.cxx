#include "sphericalharmonic.hxx"

#include <assert.h>
#include <iostream>
#include <math.h>

#include <loop_device.hxx>

namespace Multipole {
using namespace Loop;
using namespace std;

static const CCTK_REAL PI = acos(-1.0);

static inline double factorial(int n) {
  double returnval = 1;
  for (int i = n; i >= 1; i--) {
    returnval *= i;
  }
  return returnval;
}

static inline double combination(int n, int m) {
  // Binomial coefficient is undefined if these conditions do not hold
  assert(n >= 0);
  assert(m >= 0);
  assert(m <= n);
  return factorial(n) / (factorial(m) * factorial(n - m));
}

static inline int imin(int a, int b) { return a < b ? a : b; }

static inline int imax(int a, int b) { return a > b ? a : b; }

void SphericalHarmonic(int s, int l, int m, CCTK_REAL th, CCTK_REAL ph,
                       CCTK_REAL *reY, CCTK_REAL *imY) {
  double all_coeff = 0, sum = 0;
  all_coeff = pow(-1.0, m);
  all_coeff *= sqrt(factorial(l + m) * factorial(l - m) * (2 * l + 1) /
                    (4. * PI * factorial(l + s) * factorial(l - s)));
  sum = 0.;
  for (int i = imax(m - s, 0); i <= imin(l + m, l - s); i++) {
    double sum_coeff = combination(l - s, i) * combination(l + s, i + s - m);
    sum += sum_coeff * pow(-1.0, l - i - s) * pow(cos(th / 2.), 2 * i + s - m) *
           pow(sin(th / 2.), 2 * (l - i) + m - s);
  }
  *reY = all_coeff * sum * cos(m * ph);
  *imY = all_coeff * sum * sin(m * ph);
}

void HarmonicSetup(int s, int l, int m, int array_size,
                   const std::vector<CCTK_REAL> th,
                   const std::vector<CCTK_REAL> ph, std::vector<CCTK_REAL> &reY,
                   std::vector<CCTK_REAL> &imY) {
  for (int i = 0; i < array_size; i++) {
    SphericalHarmonic(s, l, m, th[i], ph[i], &reY[i], &imY[i]);
  }
}

} // namespace Multipole
