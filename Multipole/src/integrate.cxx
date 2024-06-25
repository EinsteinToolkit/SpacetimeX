#include "integrate.hxx"
#include "utils.hxx"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

namespace Multipole {

/**
 * We will want to integrate functions F(th,ph) from th = 0 to pi, ph = 0
 * to 2 pi with a weighting function sin(th).  Alternatively, we might
 * want to use u = cos(th) as the variable, in which case we will go from
 * u = -1 to 1 and ph = 0 to 2 pi.  For simplicity, we implement an
 * integration routine with a weight function of 1, and require the user
 * to multiply the integrand by their own weight function.  We divide the
 * interval [a,b] into nx subintervals of spacing h = (b-a)/nx.  These
 * have coordinates [x_i-1, xi] where x_i = x_0 + i h. So i runs from 0
 * to nx.  We require the function to integrate at the points F[x_i, y_i].
 * We have x_0 = a and x_n = b.  Check: x_n = x_0 + n (b-a)/n = a + b - a = b.
 * Good.
 *
 * If we are given these points in an array, we also need the width and
 * height of the array.  To get an actual integral, we also need the grid
 * spacing hx and hy, but these are just multiplied by the result to give
 * the integral.
 */

#define idx(xx, yy)                                                            \
  (assert((xx) <= nx), assert((xx) >= 0), assert((yy) <= ny),                  \
   assert((yy) >= 0), ((xx) + (yy) * (nx + 1)))

////////////////////////////////////////////////////////////////////////////////
// Hard coded 2D integrals
////////////////////////////////////////////////////////////////////////////////

CCTK_REAL Midpoint2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                             CCTK_REAL hy) {
  CCTK_REAL integrand_sum = 0.0;
  int ix = 0, iy = 0;

  assert(nx > 0);
  assert(ny > 0);
  assert(f);

  for (iy = 0; iy <= ny; iy++)
    for (ix = 0; ix <= nx; ix++)
      integrand_sum += f[idx(ix, iy)];

  return hx * hy * integrand_sum;
}

CCTK_REAL Trapezoidal2DIntegral(CCTK_REAL const *f, int nx, int ny,
                                CCTK_REAL hx, CCTK_REAL hy) {
  CCTK_REAL integrand_sum = 0.0;
  int ix = 0, iy = 0;

  assert(nx > 0);
  assert(ny > 0);
  assert(f);

  // Corners
  integrand_sum +=
      f[idx(0, 0)] + f[idx(nx, 0)] + f[idx(0, ny)] + f[idx(nx, ny)];

  // Edges
  for (ix = 1; ix <= nx - 1; ix++)
    integrand_sum += 2 * f[idx(ix, 0)] + 2 * f[idx(ix, ny)];

  for (iy = 1; iy <= ny - 1; iy++)
    integrand_sum += 2 * f[idx(0, iy)] + 2 * f[idx(nx, iy)];

  // Interior
  for (iy = 1; iy <= ny - 1; iy++)
    for (ix = 1; ix <= nx - 1; ix++)
      integrand_sum += 4 * f[idx(ix, iy)];

  return (double)1 / (double)4 * hx * hy * integrand_sum;
}

CCTK_REAL Simpson2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                            CCTK_REAL hy) {
  CCTK_REAL integrand_sum = 0;
  int ix = 0, iy = 0;

  assert(nx > 0);
  assert(ny > 0);
  assert(f);
  assert(nx % 2 == 0);
  assert(ny % 2 == 0);

  int px = nx / 2;
  int py = ny / 2;

  // Corners
  integrand_sum +=
      f[idx(0, 0)] + f[idx(nx, 0)] + f[idx(0, ny)] + f[idx(nx, ny)];

  // Edges
  for (iy = 1; iy <= py; iy++)
    integrand_sum += 4 * f[idx(0, 2 * iy - 1)] + 4 * f[idx(nx, 2 * iy - 1)];

  for (iy = 1; iy <= py - 1; iy++)
    integrand_sum += 2 * f[idx(0, 2 * iy)] + 2 * f[idx(nx, 2 * iy)];

  for (ix = 1; ix <= px; ix++)
    integrand_sum += 4 * f[idx(2 * ix - 1, 0)] + 4 * f[idx(2 * ix - 1, ny)];

  for (ix = 1; ix <= px - 1; ix++)
    integrand_sum += 2 * f[idx(2 * ix, 0)] + 2 * f[idx(2 * ix, ny)];

  // Interior
  for (iy = 1; iy <= py; iy++)
    for (ix = 1; ix <= px; ix++)
      integrand_sum += 16 * f[idx(2 * ix - 1, 2 * iy - 1)];

  for (iy = 1; iy <= py - 1; iy++)
    for (ix = 1; ix <= px; ix++)
      integrand_sum += 8 * f[idx(2 * ix - 1, 2 * iy)];

  for (iy = 1; iy <= py; iy++)
    for (ix = 1; ix <= px - 1; ix++)
      integrand_sum += 8 * f[idx(2 * ix, 2 * iy - 1)];

  for (iy = 1; iy <= py - 1; iy++)
    for (ix = 1; ix <= px - 1; ix++)
      integrand_sum += 4 * f[idx(2 * ix, 2 * iy)];

  return ((double)1 / (double)9) * hx * hy * integrand_sum;
}

// See: J.R. Driscoll and D.M. Healy Jr., Computing Fourier transforms
// and convolutions on the 2-sphere. Advances in Applied Mathematics,
// 15(2):202–250, 1994.
CCTK_REAL DriscollHealy2DIntegral(CCTK_REAL const *const f, int const nx,
                                  int const ny, CCTK_REAL const hx,
                                  CCTK_REAL const hy) {
  assert(f);
  assert(nx >= 0);
  assert(ny >= 0);
  assert(nx % 2 == 0);

  CCTK_REAL integrand_sum = 0.0;

  // Skip the poles (ix=0 and ix=nx), since the weight there is zero
  // anyway
#pragma omp parallel for reduction(+ : integrand_sum)
  for (int ix = 1; ix < nx; ++ix) {

    // These weights lead to an almost spectral convergence
    CCTK_REAL const theta = M_PI * ix / nx;
    CCTK_REAL weight = 0.0;
    for (int l = 0; l < nx / 2; ++l) {
      weight += sin((2 * l + 1) * theta) / (2 * l + 1);
    }
    weight *= 4.0 / M_PI;

    CCTK_REAL local_sum = 0.0;
    // Skip the last point (iy=ny), since we assume periodicity and
    // therefor it has the same value as the first point. We don't use
    // weights in this direction, which leads to spectral convergence.
    // (Yay periodicity!)
    for (int iy = 0; iy < ny; ++iy) {
      local_sum += f[idx(ix, iy)];
    }

    integrand_sum += weight * local_sum;
  }

  return hx * hy * integrand_sum;
}

void Integrate(int array_size, int nthetap,
               const std::vector<CCTK_REAL> array1r,
               const std::vector<CCTK_REAL> array1i,
               const std::vector<CCTK_REAL> array2r,
               const std::vector<CCTK_REAL> array2i, CCTK_REAL const th[],
               CCTK_REAL const ph[], CCTK_REAL *outre, CCTK_REAL *outim) {
  DECLARE_CCTK_PARAMETERS;

  int il = Index_2d(0, 0, ntheta);
  int iu = Index_2d(1, 0, ntheta);
  CCTK_REAL dth = th[iu] - th[il];
  iu = Index_2d(0, 1, ntheta);
  CCTK_REAL dph = ph[iu] - ph[il];

  static CCTK_REAL *fr = 0;
  static CCTK_REAL *fi = 0;
  static bool allocated_memory = false;

  // Construct an array for the real integrand
  if (!allocated_memory) {
    fr = new CCTK_REAL[array_size];
    fi = new CCTK_REAL[array_size];
    allocated_memory = true;
  }

  // the below calculations take the integral of conj(array1)*array2*sin(th)
  for (int i = 0; i < array_size; i++) {
    fr[i] = (array1r[i] * array2r[i] + array1i[i] * array2i[i]) * sin(th[i]);
    fi[i] = (array1r[i] * array2i[i] - array1i[i] * array2r[i]) * sin(th[i]);
  }

  if (CCTK_Equals(integration_method, "midpoint")) {
    *outre = Midpoint2DIntegral(fr, ntheta, nphi, dth, dph);
    *outim = Midpoint2DIntegral(fi, ntheta, nphi, dth, dph);
  } else if (CCTK_Equals(integration_method, "trapezoidal")) {
    *outre = Trapezoidal2DIntegral(fr, ntheta, nphi, dth, dph);
    *outim = Trapezoidal2DIntegral(fi, ntheta, nphi, dth, dph);
  } else if (CCTK_Equals(integration_method, "Simpson")) {
    if (nphi % 2 != 0 || ntheta % 2 != 0) {
      CCTK_WARN(
          CCTK_WARN_ABORT,
          "The Simpson integration method requires even ntheta and even nphi");
    }
    *outre = Simpson2DIntegral(fr, ntheta, nphi, dth, dph);
    *outim = Simpson2DIntegral(fi, ntheta, nphi, dth, dph);
  } else if (CCTK_Equals(integration_method, "DriscollHealy")) {
    if (ntheta % 2 != 0) {
      CCTK_WARN(CCTK_WARN_ABORT,
                "The Driscoll&Healy integration method requires even ntheta");
    }
    *outre = DriscollHealy2DIntegral(fr, ntheta, nphi, dth, dph);
    *outim = DriscollHealy2DIntegral(fi, ntheta, nphi, dth, dph);
  } else {
    CCTK_WARN(CCTK_WARN_ABORT, "internal error");
  }
}

////////////////////////////////////////////////////////////////////////////////
// 1D integrals
////////////////////////////////////////////////////////////////////////////////
//
// static CCTK_REAL Simpson1DIntegral(CCTK_REAL const *f, int n, CCTK_REAL h)
//{
//  CCTK_REAL integrand_sum = 0;
//  int i = 0;
//
//  assert(f);
//  assert(n > 0);
//  assert(n % 2 == 0);
//
//  int p = n / 2;
//
//  integrand_sum += f[0] + f[n];
//
//  for (i = 1; i <= p-1; i++)
//    integrand_sum += 4 * f[2*i-1] + 2 * f[2*i];
//
//  integrand_sum += 4 * f[2*p-1];
//
//  return 1.0/3.0 * h * integrand_sum;
//}
//
//// 2D integral built up from 1D
//
// static CCTK_REAL Composite2DIntegral(CCTK_REAL const *f, int nx, int ny,
// CCTK_REAL hx, CCTK_REAL hy)
//{
//  CCTK_REAL integrand_sum = 0;
//
//  assert(nx > 0); assert(ny > 0); assert (f);
//  assert(nx % 2 == 0);
//  assert(ny % 2 == 0);
//
//  CCTK_REAL *g = new CCTK_REAL[ny+1];
//
//  for (int i = 0; i <= ny; i++)
//  {
//    g[i] = Simpson1DIntegral(&f[idx(0,i)], nx, hx);
//  }
//
//  integrand_sum = Simpson1DIntegral(g, ny, hy);
//  delete [] g;
//  return integrand_sum;
//}

} // namespace Multipole
