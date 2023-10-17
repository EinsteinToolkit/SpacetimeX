#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



static
CCTK_REAL min (CCTK_REAL const x, CCTK_REAL const y)
{
  return x < y ? x : y;
}

static
CCTK_REAL max (CCTK_REAL const x, CCTK_REAL const y)
{
  return x > y ? x : y;
}



void SphericalSurface_Set (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_SphericalSurface_Set;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_REAL const pi = 3.1415926535897932384626433832795028841971693993751;
  
  int n;
  int i, j;
  
  
  
  for (n=0; n<nsurfaces; ++n) {
    
    if (set_spherical[n]) {
      
      sf_active[n] = +1;
      sf_valid[n] = +1;
      
      sf_area[n] = 4 * pi * pow (radius[n], 2);
      
      sf_mean_radius[n] = radius[n];
      
      sf_centroid_x[n] = origin_x[n];
      sf_centroid_y[n] = origin_y[n];
      sf_centroid_z[n] = origin_z[n];
      
      sf_quadrupole_xx[n] = 0.0;
      sf_quadrupole_xy[n] = 0.0;
      sf_quadrupole_xz[n] = 0.0;
      sf_quadrupole_yy[n] = 0.0;
      sf_quadrupole_yz[n] = 0.0;
      sf_quadrupole_zz[n] = 0.0;
      
      sf_min_radius[n] = radius[n];
      sf_max_radius[n] = radius[n];
      
      sf_min_x[n] = origin_x[n] - radius[n];
      sf_min_y[n] = origin_y[n] - radius[n];
      sf_min_z[n] = origin_z[n] - radius[n];
      sf_max_x[n] = origin_x[n] + radius[n];
      sf_max_y[n] = origin_y[n] + radius[n];
      sf_max_z[n] = origin_z[n] + radius[n];
      
      for (j=0; j<sf_nphi[n]; ++j) {
        for (i=0; i<sf_ntheta[n]; ++i) {
          int const ind = i + maxntheta * (j + maxnphi * n);
          sf_radius[ind] = radius[n];
        }
      }
      
      sf_origin_x[n] = origin_x[n];
      sf_origin_y[n] = origin_y[n];
      sf_origin_z[n] = origin_z[n];
      
    } else if (set_elliptic[n]) {
      
      /*
       * Equation for an ellipsoid:
       *    x^2 / rx^2 + y^2 / ry^2 + z^2/ rz^2 = 1
       * where rx, ry, rz are the radii in the x, y, and z directions.
       *
       * With
       *    x = r (sin theta) (cos phi)
       *    y = r (sin theta) (sin phi)
       *    z = r (cos theta)
       *
       * we get the solution
       *    1/r^2 = x^2/rx^2 + y^2/ry^2 + z^2/rz^2
       *
       * With the form
       *    x A x = 1
       * we obtain
       *    A = diag [1/rx^2, 1/ry^2, 1/rz^2]
       */
      
      CCTK_REAL const rx2 = pow (radius_x[n], 2);
      CCTK_REAL const ry2 = pow (radius_y[n], 2);
      CCTK_REAL const rz2 = pow (radius_z[n], 2);
      
      sf_active[n] = +1;
      sf_valid[n] = +1;
      
      sf_area[n] = 0 * 4 * pi * pow (radius[n], 2);
      
      sf_mean_radius[n] = 0 * radius[n];
      
      sf_centroid_x[n] = origin_x[n];
      sf_centroid_y[n] = origin_y[n];
      sf_centroid_z[n] = origin_z[n];
      
      sf_quadrupole_xx[n] = 1.0 / rz2;
      sf_quadrupole_xy[n] = 0.0;
      sf_quadrupole_xz[n] = 0.0;
      sf_quadrupole_yy[n] = 1.0 / ry2;
      sf_quadrupole_yz[n] = 0.0;
      sf_quadrupole_zz[n] = 1.0 / rx2;
      
      sf_min_radius[n] = min (radius_x[n], min (radius_y[n], radius_z[n]));
      sf_max_radius[n] = max (radius_x[n], max (radius_y[n], radius_z[n]));
      
      sf_min_x[n] = origin_x[n] - radius_x[n];
      sf_min_y[n] = origin_y[n] - radius_y[n];
      sf_min_z[n] = origin_z[n] - radius_z[n];
      sf_max_x[n] = origin_x[n] + radius_x[n];
      sf_max_y[n] = origin_y[n] + radius_y[n];
      sf_max_z[n] = origin_z[n] + radius_z[n];
      
      for (j=0; j<sf_nphi[n]; ++j) {
        for (i=0; i<sf_ntheta[n]; ++i) {
          int const ind = i + maxntheta * (j + maxnphi * n);
          CCTK_REAL const theta = sf_origin_theta[n] + i * sf_delta_theta[n];
          CCTK_REAL const phi = sf_origin_phi[n] + j * sf_delta_phi[n];
          CCTK_REAL const x2 = pow (sin(theta) * cos(phi), 2);
          CCTK_REAL const y2 = pow (sin(theta) * sin(phi), 2);
          CCTK_REAL const z2 = pow (cos(theta)           , 2);
          sf_radius[ind] = 1.0 / sqrt (x2 / rx2 + y2 / ry2 + z2 / rz2);
        }
      }
      
      sf_origin_x[n] = origin_x[n];
      sf_origin_y[n] = origin_y[n];
      sf_origin_z[n] = origin_z[n];
      
    }
    
  } /* for n */
}
