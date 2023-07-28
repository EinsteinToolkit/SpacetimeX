#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

//#include <mpi.h>


// local functions
static bool all_sf_names_are_unique(CCTK_INT warnlevel);
static int get_reflevel (cGH const * const cctkGH);


static int get_reflevel (cGH const * const cctkGH)
{
  if (CCTK_IsFunctionAliased ("GetRefinementLevel")) {
    return GetRefinementLevel (cctkGH);
  } else {
    return 0;
  }
}


// simple translation function from user-friendly name to numerical Id
// if name is not empty returns the matching name and -1 if no name matches
// if name is empty returns fallbackid
extern "C" CCTK_INT SphericalSurface_IdFromName (CCTK_INT fallbackid, CCTK_POINTER_TO_CONST sfname)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT retval = -1;

  // re-check for unique names since the parameters are steerable
  if(!all_sf_names_are_unique(CCTK_WARN_ALERT))
    CCTK_WARN(CCTK_WARN_ABORT, "Not all spherical surface names are unique."); 

  if (CCTK_Equals(static_cast<CCTK_STRING>(sfname), "")) {

    retval = fallbackid;
  } else {

    for (CCTK_INT n=0; n<nsurfaces; ++n) {
      
      if (CCTK_Equals(static_cast<CCTK_STRING>(sfname), name[n])) {
        
        retval = n;
        break;
      }
    }

    if (retval == -1) {
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Did not find spherical surface with name '%s'.", 
                 static_cast<CCTK_STRING>(sfname));
    }
  }

  return retval;
}

// utility function to check if all SphericalSurface names are unique (or empty)
static bool all_sf_names_are_unique(CCTK_INT warnlevel)
{
  DECLARE_CCTK_PARAMETERS;
  bool all_unique = true;

  for (int i=0; i<nsurfaces; ++i) {
    if (CCTK_Equals(static_cast<CCTK_STRING>(name[i]), ""))
      continue;
    for (int j=i+1; j<nsurfaces; ++j) {
      if (CCTK_Equals(static_cast<CCTK_STRING>(name[i]), static_cast<CCTK_STRING>(name[j]))) {
        CCTK_VWarn(warnlevel, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Duplicate names for spherical surfaces %d and %d: %s == %s", 
                   i,j, name[i], name[j]);
        all_unique = false;
      }
    }
  }

  return all_unique;
}

// check early that all SphericalSurface names are valid
extern "C" void SphericalSurface_ParamCheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_SphericalSurface_ParamCheck;

  if (!all_sf_names_are_unique(CCTK_WARN_ALERT)) {
    CCTK_PARAMWARN("Not all spherical surface names are unique.");
  }
}
      
//TODO:  overlap with finer levels still seems to be incorrect!

extern "C" void SphericalSurface_Setup (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_SphericalSurface_Setup;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_REAL const pi = 4 * atan(CCTK_REAL(1.0));
  
  int ierr;
  
  
  
  if (nsurfaces == 0) return;
  
  
  
  int const group = CCTK_GroupIndex ("SphericalSurface::sf_radius");
  assert (group>=0);
  cGroup groupinfo;
  ierr = CCTK_GroupData (group, &groupinfo);
  assert (!ierr);
  cGroupDynamicData groupdata;
  ierr = CCTK_GroupDynamicData (cctkGH, group, &groupdata);
  assert (!ierr);
  
  
  
  for (int n=0; n<nsurfaces; ++n) {
    
    if (!auto_res[n])
    {
      // set resolution according to given parameters
    
      // internal consistency checks
      assert (groupdata.dim == 2);
      assert (groupdata.gsh[0] >= ntheta[n]);
      assert (groupdata.gsh[1] >= nphi[n]);
      
      assert (ntheta[n] >= 3*nghoststheta[n] && ntheta[n] <= maxntheta);
      assert (nphi[n] >= 3*nghostsphi[n] && nphi[n] <= maxnphi);
      
      
      
      // copy parameters into grid functions
      sf_ntheta[n] = ntheta[n];
      sf_nphi[n] = nphi[n];
      sf_nghoststheta[n] = nghoststheta[n];
      sf_nghostsphi[n] = nghostsphi[n];
    }
    else
    {
      // we already have an approximate delta-spacing calculated in "SetupRes"
      // based on that, we calculate the number of necessary gridpoints...
      CCTK_WARN(1, "AutoRes will NOT work with checkpoint/recovery!");
      
      CCTK_REAL theta_range = pi;
      CCTK_REAL phi_range = 2*pi;
      
      if (symmetric_x[n])
         phi_range /= 2;
         
      if (symmetric_y[n])
         phi_range /= 2;
      
      if (symmetric_z[n])
         theta_range /= 2;
      
      sf_ntheta[n] = lrint(theta_range/sf_delta_theta_estimate[n]) + 2*nghoststheta[n];
      sf_nphi[n]   = lrint(phi_range  /sf_delta_phi_estimate[n]  ) + 2*nghostsphi[n];
      
      // ...make number of theta-gridpoints odd...
      if (sf_ntheta[n] % 2 != 1)
        sf_ntheta[n] += 1;
      
      // ... and make inner phi-gridpoints divisible by 4 since some thorns require that (e.g. IsolatedHorizon)
      sf_nphi[n]   += 4 - (sf_nphi[n]-2*nghostsphi[n]) % 4 + 2*nghostsphi[n];
      
      // consistency check
      if (sf_ntheta[n] < 3*nghoststheta[n]) 
      {
         CCTK_WARN(1, "AutoRes: Determined sf_ntheta[n] < 3*nghoststheta[n]! Setting sf_ntheta[n] = 3*nghoststheta[n]!");
         sf_ntheta[n] = 3*nghoststheta[n];
      }
      if (sf_ntheta[n] > maxntheta)
      {
         CCTK_WARN(1, "AutoRes: Determined sf_ntheta[n] > maxntheta! You should increase maxntheta or auto_res_ratio[n]! Setting sf_ntheta[n] = maxntheta!");
         sf_ntheta[n] = maxntheta;
      }
      if (sf_nphi[n] < 3*nghostsphi[n]) 
      {
         CCTK_WARN(1, "AutoRes: Determined sf_nphi[n] < 3*nghostsphi[n]! Setting sf_nphi[n] = 3*nghostsphi[n]!");
         sf_nphi[n] = 3*nghostsphi[n];
      }
      if (sf_nphi[n] > maxnphi)
      {
         CCTK_WARN(1, "AutoRes: Determined sf_nphi[n] > maxnphi! You should increase maxnphi or auto_res_ratio[n]! Setting sf_nphi[n] = maxnphi!");
         sf_nphi[n] = maxnphi;
      }
      
      // check, if there are enough maxntheta,maxnphi gridpoints to remove symmetries
      // (some thorns require that)
      if (symmetric_z[n])
      {
         if (2*(sf_ntheta[n]-2*nghoststheta[n])+2*nghoststheta[n] > maxntheta)
	    CCTK_WARN(0, "maxntheta is not big enough to remove z-symmetry (some thorns may require that)!");
      }
      
      if (symmetric_x[n])
      {
         if (symmetric_y[n])
         {
            if (4*(sf_nphi[n]-2*nghostsphi[n])+2*nghostsphi[n] > maxnphi)
            {
               CCTK_WARN(0, "maxnphi is not big enough to remove x-symmetry and y-symmetry (some thorns may require that)!");
            }
         }
         else
         {
            if (2*(sf_nphi[n]-2*nghostsphi[n])+2*nghostsphi[n] > maxnphi)
            {
               CCTK_WARN(0, "maxnphi is not big enough to remove x-symmetry and y-symmetry (some thorns may require that)!");
            }
         }
      }
      else
      {
         if (symmetric_y[n])
         {
            if (2*(sf_nphi[n]-2*nghostsphi[n])+2*nghostsphi[n] > maxnphi)
            {
               CCTK_WARN(0, "maxnphi is not big enough to remove x-symmetry and y-symmetry (some thorns may require that)!");
            }
         }
      }

      
      sf_nghoststheta[n] = nghoststheta[n];
      sf_nghostsphi[n] = nghostsphi[n];
      
      if (verbose)
      {
        CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface[%d]: setting sf_ntheta = %d", n, int(sf_ntheta[n]));
        CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface[%d]: setting sf_nphi   = %d", n, int(sf_nphi[n]));
      }
    }
    
    if (verbose)
    {
       CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface[%d]: finest refinement-level that fully contains this surface:  sf_maxreflevel = %d", n, int(sf_maxreflevel[n]));
       CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface[%d]: finest refinement-level that overlaps with this surface :  sf_minreflevel = %d", n, int(sf_minreflevel[n]));
    }
    
    // coordinates in the theta direction
    // avoid_sf_origin_theta = 1
    if (symmetric_z[n]) {
      
      // upper hemisphere: z>=0, theta in (0, pi/2)
      sf_delta_theta[n] = pi/2 / (sf_ntheta[n] - 2*nghoststheta[n] - 0.5);
      sf_origin_theta[n] = - (nghoststheta[n] - 0.5) * sf_delta_theta[n];
      
    } else {
      
      // both hemispheres: theta in (0, pi)
      sf_delta_theta[n] = pi / (sf_ntheta[n] - 2*nghoststheta[n]);
      sf_origin_theta[n] = - (nghoststheta[n] - 0.5) * sf_delta_theta[n];
      
    }
    
    
    
    // coordinates in the phi direction
    // avoid_sf_origin_phi = 0
    if (symmetric_x[n]) {
      if (symmetric_y[n]) {
        
        // one quadrant: x>=0, y>=0, phi in [0, pi/2]
        assert (sf_nphi[n] - 2*nghostsphi[n] >= 1);
        sf_delta_phi[n] = pi/2 / (sf_nphi[n] - 2*nghostsphi[n] - 1);
        sf_origin_phi[n] = - nghostsphi[n] * sf_delta_phi[n];
        
      } else {
        
        // two quadrants: x>=0, phi in [-pi/2, pi/2]
        assert (sf_nphi[n] - 2*nghostsphi[n] >= 2);
        sf_delta_phi[n] = pi / (sf_nphi[n] - 2*nghostsphi[n] - 1);
        sf_origin_phi[n] = - pi/2 - nghostsphi[n] * sf_delta_phi[n];
        
      }
    } else {
      if (symmetric_y[n]) {
        
        // two quadrants: y>=0, phi in [0, pi]
        assert (sf_nphi[n] - 2*nghostsphi[n] >= 2);
        sf_delta_phi[n] = pi / (sf_nphi[n] - 2*nghostsphi[n] - 1);
        sf_origin_phi[n] = - nghostsphi[n] * sf_delta_phi[n];
        
      } else {
        
        // all quadrants: phi in [0, 2pi)
        assert (sf_nphi[n] - 2*nghostsphi[n] >= 4);
        sf_delta_phi[n] = 2*pi / (sf_nphi[n] - 2*nghostsphi[n]);
        sf_origin_phi[n] = - nghostsphi[n] * sf_delta_phi[n];
        
      }
    }
    
    if (verbose) {
      CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface[%d]: sf_delta_theta = %.6f, sf_delta_phi = %.6f\n", n, sf_delta_theta[n], sf_delta_phi[n]);
    }
    
    // mark surface as uninitialised
    sf_active[n] = 0;
    sf_valid[n] = 0;
    
  } // for n
}



extern "C" void SphericalSurface_SetupRes (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_SphericalSurface_SetupRes;
  DECLARE_CCTK_PARAMETERS;
  
  static bool first_call = true;
  
  int min_reduction_handle = CCTK_ReductionArrayHandle("minimum");
  int max_reduction_handle = CCTK_ReductionArrayHandle("maximum");

  if (min_reduction_handle < 0)
     CCTK_WARN(0, "Cannot get reduction handle for minimum operation.");
  if (max_reduction_handle < 0)
     CCTK_WARN(0, "Cannot get reduction handle for maximum operation.");
  
  for (int n=0; n<nsurfaces; ++n) 
  {
    
    // try to find out with which Carpet::reflevel completely contains
    // the current surface and will intersect and and store the
    // result so that other thorns such as the kick-thorn
    // can decide which reflevel to use
     
    // set resolution according to Radius and Cartesian resolution

    if (!auto_res[n])
        continue;
      
    CCTK_REAL my_radius;
      
    if (!set_spherical[n] && !set_elliptic[n])
    {
       // we have a problem...without a radius we cannot estimate a proper resolution
       // Here, we simply assume r=1 if radius[n] == 0.
       // If you want a better estimate of the resolution, set the radius parameter to an initial guess!
       
       if (radius[n] == 0)
          my_radius = 1;
       else
	  my_radius = radius[n];
    }
    else
    {
       my_radius = radius[n];
    }
      
    if (verbose) {
      CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface: myradius = %f", double(my_radius));
    }
    // get theta/phi delta-spacing by looking at Cartesian grid
      
    // first test, if we have an overlap with a certain reflevel
    // we do this by calculating a sphere based on the corners of the Cartesian local patch
    // on this processor.
    //unsigned long l[2];
      
    //l[0] = CCTK_GFINDEX3D(cctkGH, 0,             0,             0);
    //l[1] = CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1);
      
    bool is_overlapping = false;
    CCTK_REAL cart_radius = 0;    // radius of Cart. grid
    CCTK_REAL cart_origin[3];     // origin of Cart. grid
    
    //CCTK_ORIGIN_SPACE(0) + (cctk_lbnd[0]) * CCTK_DELTA_SPACE(0) - (CCTK_ORIGIN_SPACE(0) + (cctk_lbnd[0] + cctk_lsh[0]) * CCTK_DELTA_SPACE(0))
    
    // calculate radius of Cartesian grid
    cart_origin[0] = 0.5 * cctk_lsh[0] * CCTK_DELTA_SPACE(0); //(x[l[1]] - x[l[0]]) / 2.0;
    cart_origin[1] = 0.5 * cctk_lsh[1] * CCTK_DELTA_SPACE(1); //(y[l[1]] - y[l[0]]) / 2.0;
    cart_origin[2] = 0.5 * cctk_lsh[2] * CCTK_DELTA_SPACE(2); //(z[l[1]] - z[l[0]]) / 2.0;
      
    cart_radius = sqrt(cart_origin[0]*cart_origin[0] + cart_origin[1]*cart_origin[1] + cart_origin[2]*cart_origin[2]);
      
    if (verbose) {
      CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface: cart_radius = %f, cart_origin = %f, %f, %f", double(cart_radius), double(cart_origin[0]), double(cart_origin[1]), double(cart_origin[2]));
    }

    // set origin of Cartesian grid
    cart_origin[0] += CCTK_ORIGIN_SPACE(0) + cctk_lbnd[0] * CCTK_DELTA_SPACE(0); //x[l[0]];
    cart_origin[1] += CCTK_ORIGIN_SPACE(1) + cctk_lbnd[1] * CCTK_DELTA_SPACE(1); //y[l[0]];
    cart_origin[2] += CCTK_ORIGIN_SPACE(2) + cctk_lbnd[2] * CCTK_DELTA_SPACE(2); //z[l[0]];
      
    // norm of displacement vector between origin of Cart. grid and origin of SphericalSurface
    CCTK_REAL dist = sqrt( (cart_origin[0]-origin_x[n]) * (cart_origin[0]-origin_x[n])
                         + (cart_origin[1]-origin_y[n]) * (cart_origin[1]-origin_y[n])
                         + (cart_origin[2]-origin_z[n]) * (cart_origin[2]-origin_z[n]));
      
    // does our spherical shell overlap with the processor's local grid patch? 
    if (my_radius+cart_radius >= dist &&  // Cart. grid and SphericalSurface overlap
        cart_radius+dist > my_radius  &&  // Cart. grid patch is not fully contained in this SphericalSurface
        -cart_radius+dist < my_radius)    // Cart. grid patch is not outside of the spherical shell
       is_overlapping = true;
         
    // get resolution
    if (is_overlapping)
    {
       if (CCTK_Equals(auto_res_grid[n], "overlap"))
       {
          sf_delta_theta_estimate[n] = auto_res_ratio[n] / my_radius * sqrt(CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0) 
                                                                          + CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1) 
                                                                          + CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2));
         
          sf_delta_phi_estimate[n] = auto_res_ratio[n] / my_radius * sqrt(CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0) 
                                                                        + CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1) 
                                                                        + CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2));
       }
       
       // store reflevel
       sf_minreflevel[n] = get_reflevel(cctkGH);
    }
    else if (first_call)
    {
       // set to very high value initially (otherwise would be zero
       // and reported as minimum value in following MPI communication)
       sf_delta_theta_estimate[n] = 10000;
       sf_delta_phi_estimate[n] = 10000; 
    }
    
    if (CCTK_Equals(auto_res_grid[n], "overlap"))
    {
      // find minimum value of "sf_delta_phi/theta" over all processors
      CCTK_REAL dummy = 0;
      //MPI_Allreduce(&sf_delta_theta_estimate[n], &dummy, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      CCTK_ReduceLocScalar(cctkGH, -1, min_reduction_handle, &sf_delta_theta_estimate[n], &dummy, CCTK_VARIABLE_REAL);
      sf_delta_theta_estimate[n] = dummy;
         
      dummy = 0;
      //MPI_Allreduce(&sf_delta_phi_estimate[n], &dummy, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      CCTK_ReduceLocScalar(cctkGH, -1, min_reduction_handle, &sf_delta_phi_estimate[n], &dummy, CCTK_VARIABLE_REAL);
      sf_delta_phi_estimate[n] = dummy;
      
      /*
      // input arrays and output values
      const CCTK_INT input_array_variable_indices[2]
        = { CCTK_VarIndex("SphericalSurface::sf_delta_theta_estimate"),
            CCTK_VarIndex("SphericalSurface::sf_delta_phi_estimate") };
      const CCTK_INT output_value_type_codes[2]
        = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL };
      void *const output_numbers[2]
        = { (void *) output_for_real_values,
            (void *) output_for_real_values };
      
      const int status
         = CCTK_ReduceGridArrays(GH,
                                 0,
                                 2, input_array_variable_indices,
                                 2, output_value_type_codes,
                                 output_values);
      */
    }
    
    // find maximum value of "sf_minreflevel" over all processors
    int idummy;
    //MPI_Allreduce(&sf_minreflevel[n], &idummy, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    CCTK_ReduceLocScalar(cctkGH, -1, max_reduction_handle, &sf_minreflevel[n], &idummy, CCTK_VARIABLE_INT);
    sf_minreflevel[n] = idummy;



    // find reflevel that completely contains this surface
    // FIXME: This algorithm does not take into account multiple
    //        disconnected refinement regions.
    
    
    //safety margin of 6 points around sphere
    CCTK_INT const np = 6;
   
    CCTK_REAL max_radius = my_radius + ( sqrt(  np*np*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0) 
                                              + np*np*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1)
                                              + np*np*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2) ) ); 

    if (verbose) {
    //CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface: calculated max_radius = %.12f\n", max_radius);
      CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface: buffer_radius = %.12f", double (sqrt(  np*np*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0) 
                                                                                            + np*np*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1)
                                                                                            + np*np*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2) )));
    }

   /*
    unsigned long k[8];
      
    k[0] = CCTK_GFINDEX3D(cctkGH, 0,             0,             0);
    k[1] = CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, 0,             0);
    k[2] = CCTK_GFINDEX3D(cctkGH, 0,             cctk_lsh[1]-1, 0);
    k[3] = CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, 0);
    k[4] = CCTK_GFINDEX3D(cctkGH, 0,             0,             cctk_lsh[2]-1);
    k[5] = CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, 0,             cctk_lsh[2]-1);
    k[6] = CCTK_GFINDEX3D(cctkGH, 0,             cctk_lsh[1]-1, cctk_lsh[2]-1);
    k[7] = CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1);
    */
    
    unsigned long kx[8], ky[8], kz[8];
      
    kx[0] = 0;             ky[0] = 0;             kz[0] = 0;
    kx[1] = cctk_lsh[0]-1, ky[1] = 0;             kz[1] = 0;
    kx[2] = 0;             ky[2] = cctk_lsh[1]-1; kz[2] = 0;
    kx[3] = cctk_lsh[0]-1; ky[3] = cctk_lsh[1]-1; kz[3] = 0;
    kx[4] = 0;             ky[4] = 0;             kz[4] = cctk_lsh[2]-1;
    kx[5] = cctk_lsh[0]-1; ky[5] = 0;             kz[5] = cctk_lsh[2]-1;
    kx[6] = 0;             ky[6] = cctk_lsh[1]-1; kz[6] = cctk_lsh[2]-1;
    kx[7] = cctk_lsh[0]-1; ky[7] = cctk_lsh[1]-1; kz[7] = cctk_lsh[2]-1;
    
    
    int i;
    bool is_contained = false;
    for (i=0; i < 8; i++)
    {
       //if (max_radius < sqrt(x[k[i]]*x[k[i]] + y[k[i]]*y[k[i]] + z[k[i]]*z[k[i]]))
       //   is_contained = true;
       CCTK_REAL xx = CCTK_ORIGIN_SPACE(0) + (cctk_lbnd[0] + kx[i]) * CCTK_DELTA_SPACE(0);
       CCTK_REAL yy = CCTK_ORIGIN_SPACE(1) + (cctk_lbnd[1] + ky[i]) * CCTK_DELTA_SPACE(1);
       CCTK_REAL zz = CCTK_ORIGIN_SPACE(2) + (cctk_lbnd[2] + kz[i]) * CCTK_DELTA_SPACE(2);
       if (max_radius < sqrt(xx*xx + yy*yy + zz*zz))
          is_contained = true;
    }
   
    int max_rl = 0;
    if (is_contained)
    {
       if (CCTK_Equals(auto_res_grid[n], "fully contained"))
       {
          sf_delta_theta_estimate[n] = auto_res_ratio[n] / my_radius * sqrt(CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0) 
                                                                          + CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1) 
                                                                          + CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2));
         
          sf_delta_phi_estimate[n] = auto_res_ratio[n] / my_radius * sqrt(CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0) 
                                                                        + CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1) 
                                                                        + CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2));
       }
    
       max_rl = get_reflevel(cctkGH);
    }
    else if (first_call)
    {
       // set to zero initially
       sf_delta_theta_estimate[n] = 0;
       sf_delta_phi_estimate[n] = 0; 
    }
    
    
    // find minimum value of "sf_maxreflevel" over all processors
    //MPI_Allreduce(&max_rl, &idummy, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    CCTK_ReduceLocScalar(cctkGH, -1, min_reduction_handle, &max_rl, &idummy, CCTK_VARIABLE_INT);
    if (idummy != 0 || first_call)
	sf_maxreflevel[n] = idummy;
    
    // this condition can become true because of this approximate algorithm
    // which does not take into account disconnected ref. regions
    if (sf_maxreflevel[n] > sf_minreflevel[n])
        sf_maxreflevel[n] = sf_minreflevel[n];   // can still be incorrect!
	
    
    if (CCTK_Equals(auto_res_grid[n], "fully contained"))
    {
      // find maximum value of "sf_delta_phi/theta" over all processors
      CCTK_REAL dummy = 0;
      //MPI_Allreduce(&sf_delta_theta_estimate[n], &dummy, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      CCTK_ReduceLocScalar(cctkGH, -1, max_reduction_handle, &sf_delta_theta_estimate[n], &dummy, CCTK_VARIABLE_REAL);
      sf_delta_theta_estimate[n] = dummy;
         
      dummy = 0;
      //MPI_Allreduce(&sf_delta_phi_estimate[n], &dummy, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      CCTK_ReduceLocScalar(cctkGH, -1, max_reduction_handle, &sf_delta_phi_estimate[n], &dummy, CCTK_VARIABLE_REAL);
      sf_delta_phi_estimate[n] = dummy;
    }
    
    if (verbose) {
      CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface: sf_maxreflevel[%d] = %d,  sf_minreflevel[%d] = %d", n, int(sf_maxreflevel[n]), n, int(sf_minreflevel[n]));
      CCTK_VInfo (CCTK_THORNSTRING, "SphericalSurface: sf_delta_theta[%d] = %.6f, sf_delta_phi[%d] = %.6f\n", n, sf_delta_theta_estimate[n], n, sf_delta_phi_estimate[n]);
    }

  }
  
  first_call = false;
  
}
