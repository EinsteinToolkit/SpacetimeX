/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"

namespace WeylScal4 {


static void WeylScal4_invars_calc_4th_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = cctk_time;
  const CCTK_REAL cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_TIME;
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  /* Initialize predefined quantities */
  const CCTK_REAL p1o12dx CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o12dy CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o12dz CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o144dxdy CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o144dxdz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o144dydz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o180dx2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dx,-2);
  const CCTK_REAL p1o180dy2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dy,-2);
  const CCTK_REAL p1o180dz2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dz,-2);
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dx,-1);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dy,-1);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dz,-1);
  const CCTK_REAL p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o3600dydz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dxdy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o4dxdz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dydz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dx,-2);
  const CCTK_REAL p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dy,-2);
  const CCTK_REAL p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dz,-2);
  const CCTK_REAL p1o60dx CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o60dy CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o60dz CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o705600dydz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o840dx CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dx,-1);
  const CCTK_REAL p1o840dy CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dy,-1);
  const CCTK_REAL p1o840dz CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dz,-1);
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  //#pragma omp parallel
  CCTK_LOOP3(WeylScal4_invars_calc_4th,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL Psi0iL CCTK_ATTRIBUTE_UNUSED = Psi0i[index];
    CCTK_REAL Psi0rL CCTK_ATTRIBUTE_UNUSED = Psi0r[index];
    CCTK_REAL Psi1iL CCTK_ATTRIBUTE_UNUSED = Psi1i[index];
    CCTK_REAL Psi1rL CCTK_ATTRIBUTE_UNUSED = Psi1r[index];
    CCTK_REAL Psi2iL CCTK_ATTRIBUTE_UNUSED = Psi2i[index];
    CCTK_REAL Psi2rL CCTK_ATTRIBUTE_UNUSED = Psi2r[index];
    CCTK_REAL Psi3iL CCTK_ATTRIBUTE_UNUSED = Psi3i[index];
    CCTK_REAL Psi3rL CCTK_ATTRIBUTE_UNUSED = Psi3r[index];
    CCTK_REAL Psi4iL CCTK_ATTRIBUTE_UNUSED = Psi4i[index];
    CCTK_REAL Psi4rL CCTK_ATTRIBUTE_UNUSED = Psi4r[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    
    switch (fdOrder)
    {
      case 2:
      {
        break;
      }
      
      case 4:
      {
        break;
      }
      
      case 6:
      {
        break;
      }
      
      case 8:
      {
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL curvIrL CCTK_ATTRIBUTE_UNUSED = 4*Psi1iL*Psi3iL - 
      4*Psi1rL*Psi3rL - Psi0iL*Psi4iL + Psi0rL*Psi4rL - 3*pow(Psi2iL,2) + 
      3*pow(Psi2rL,2);
    
    CCTK_REAL curvIiL CCTK_ATTRIBUTE_UNUSED = 6*Psi2iL*Psi2rL - 
      4*(Psi1rL*Psi3iL + Psi1iL*Psi3rL) + Psi0rL*Psi4iL + Psi0iL*Psi4rL;
    
    CCTK_REAL curvJrL CCTK_ATTRIBUTE_UNUSED = 2*Psi0iL*Psi3iL*Psi3rL + 
      2*Psi1iL*Psi1rL*Psi4iL - Psi2iL*(2*Psi1rL*Psi3iL + 2*Psi1iL*Psi3rL + 
      Psi0rL*Psi4iL + Psi0iL*Psi4rL) + Psi4rL*(pow(Psi1iL,2) - pow(Psi1rL,2)) 
      + Psi2rL*(-2*Psi1iL*Psi3iL + 2*Psi1rL*Psi3rL - Psi0iL*Psi4iL + 
      Psi0rL*Psi4rL + 3*pow(Psi2iL,2)) - pow(Psi2rL,3) + 
      Psi0rL*(pow(Psi3iL,2) - pow(Psi3rL,2));
    
    CCTK_REAL curvJiL CCTK_ATTRIBUTE_UNUSED = 2*Psi2rL*(Psi1rL*Psi3iL + 
      Psi1iL*Psi3rL) + Psi0rL*(-2*Psi3iL*Psi3rL + Psi2rL*Psi4iL) - 
      2*Psi1iL*Psi1rL*Psi4rL + Psi4iL*(pow(Psi1iL,2) - pow(Psi1rL,2)) + 
      pow(Psi2iL,3) + Psi2iL*(-2*Psi1iL*Psi3iL + 2*Psi1rL*Psi3rL - 
      Psi0iL*Psi4iL + Psi0rL*Psi4rL - 3*pow(Psi2rL,2)) + 
      Psi0iL*(Psi2rL*Psi4rL + pow(Psi3iL,2) - pow(Psi3rL,2));
    
    CCTK_REAL curvJ1L CCTK_ATTRIBUTE_UNUSED = -16*(-4*Psi1iL*Psi3iL + 
      4*Psi1rL*Psi3rL + Psi0iL*Psi4iL - Psi0rL*Psi4rL + 3*pow(Psi2iL,2) - 
      3*pow(Psi2rL,2));
    
    CCTK_REAL curvJ2L CCTK_ATTRIBUTE_UNUSED = 96*(-2*(Psi0iL*Psi3iL*Psi3rL 
      + Psi1iL*Psi1rL*Psi4iL) + Psi2iL*(2*(Psi1rL*Psi3iL + Psi1iL*Psi3rL) + 
      Psi0rL*Psi4iL + Psi0iL*Psi4rL) + Psi4rL*(-pow(Psi1iL,2) + 
      pow(Psi1rL,2)) + Psi2rL*(2*Psi1iL*Psi3iL - 2*Psi1rL*Psi3rL + 
      Psi0iL*Psi4iL - Psi0rL*Psi4rL - 3*pow(Psi2iL,2)) + pow(Psi2rL,3) + 
      Psi0rL*(-pow(Psi3iL,2) + pow(Psi3rL,2)));
    
    CCTK_REAL curvJ3L CCTK_ATTRIBUTE_UNUSED = 64*(((8*Psi0iL*Psi1iL - 
      8*Psi0rL*Psi1rL)*Psi3rL - 4*Psi0iL*Psi0rL*Psi4iL)*Psi4rL + 
      12*Psi2iL*Psi2rL*(4*(Psi1rL*Psi3iL + Psi1iL*Psi3rL) - Psi0rL*Psi4iL - 
      Psi0iL*Psi4rL) + Psi1iL*Psi3iL*(-64*Psi1rL*Psi3rL - 8*Psi0iL*Psi4iL + 
      8*Psi0rL*Psi4rL) + 8*(Psi0rL*(Psi1rL*Psi3iL + Psi1iL*Psi3rL)*Psi4iL + 
      Psi0iL*Psi1rL*(Psi3rL*Psi4iL + Psi3iL*Psi4rL)) + 6*(4*Psi1iL*Psi3iL - 
      4*Psi1rL*Psi3rL - Psi0iL*Psi4iL + Psi0rL*Psi4rL)*pow(Psi2rL,2) - 
      6*pow(Psi2iL,2)*(4*Psi1iL*Psi3iL - 4*Psi1rL*Psi3rL - Psi0iL*Psi4iL + 
      Psi0rL*Psi4rL + 9*pow(Psi2rL,2)) + 9*(pow(Psi2iL,4) + pow(Psi2rL,4)) + 
      (16*pow(Psi1iL,2) - 16*pow(Psi1rL,2))*pow(Psi3iL,2) + 
      (-16*pow(Psi1iL,2) + 16*pow(Psi1rL,2))*pow(Psi3rL,2) + (pow(Psi0iL,2) - 
      pow(Psi0rL,2))*pow(Psi4iL,2) + (-pow(Psi0iL,2) + 
      pow(Psi0rL,2))*pow(Psi4rL,2));
    
    CCTK_REAL curvJ4L CCTK_ATTRIBUTE_UNUSED = 
      -640*(Psi3iL*(Psi4iL*(-2*Psi3rL*pow(Psi0iL,2) + 
      12*Psi1rL*pow(Psi1iL,2)) + 4*Psi4rL*pow(Psi1iL,3)) - 2*(5*Psi1iL*Psi3iL 
      - 5*Psi1rL*Psi3rL + Psi0iL*Psi4iL - Psi0rL*Psi4rL)*pow(Psi2rL,3) - 
      3*pow(Psi2rL,5) + ((12*Psi0iL*Psi1iL - 12*Psi0rL*Psi1rL)*Psi3rL - 
      2*Psi0iL*Psi0rL*Psi4iL - Psi4rL*pow(Psi0iL,2))*pow(Psi3iL,2) + 
      Psi4rL*(Psi4iL*(4*Psi0rL*Psi1iL*Psi1rL - 2*Psi0iL*pow(Psi1iL,2)) + 
      Psi3rL*(4*Psi0iL*Psi0rL*Psi3iL - 12*Psi1rL*pow(Psi1iL,2)) + 
      pow(Psi0rL,2)*(pow(Psi3iL,2) - pow(Psi3rL,2)) + 
      pow(Psi0iL,2)*pow(Psi3rL,2)) - 12*(Psi1iL*Psi3rL*Psi4iL*pow(Psi1rL,2) + 
      Psi3iL*(Psi1iL*Psi4rL*pow(Psi1rL,2) + (Psi0rL*Psi1iL + 
      Psi0iL*Psi1rL)*pow(Psi3rL,2))) + 
      3*(pow(Psi2rL,2)*(2*(Psi0iL*Psi3iL*Psi3rL + Psi1iL*Psi1rL*Psi4iL) + 
      Psi4rL*(pow(Psi1iL,2) - pow(Psi1rL,2)) + Psi0rL*(pow(Psi3iL,2) - 
      pow(Psi3rL,2))) + pow(Psi2iL,2)*(-2*(Psi0iL*Psi3iL*Psi3rL + 
      Psi1iL*Psi1rL*Psi4iL) + 2*Psi2rL*(5*Psi1iL*Psi3iL - 5*Psi1rL*Psi3rL + 
      Psi0iL*Psi4iL - Psi0rL*Psi4rL) + Psi4rL*(-pow(Psi1iL,2) + 
      pow(Psi1rL,2)) + 10*pow(Psi2rL,3) + Psi0rL*(-pow(Psi3iL,2) + 
      pow(Psi3rL,2)))) - 4*(Psi3iL*Psi4iL*pow(Psi1rL,3) + 
      Psi0iL*Psi1iL*pow(Psi3rL,3)) + 4*(Psi3rL*(Psi4iL*pow(Psi1iL,3) + 
      Psi4rL*pow(Psi1rL,3)) + (Psi0rL*Psi1iL + Psi0iL*Psi1rL)*pow(Psi3iL,3) + 
      Psi0rL*Psi1rL*pow(Psi3rL,3)) + (-2*Psi0iL*Psi1iL*Psi1rL + 
      Psi0rL*(-pow(Psi1iL,2) + pow(Psi1rL,2)))*pow(Psi4iL,2) + 
      (2*Psi0iL*Psi1iL*Psi1rL + Psi0rL*(pow(Psi1iL,2) - 
      pow(Psi1rL,2)))*pow(Psi4rL,2) + Psi2rL*(-4*Psi0iL*Psi0rL*Psi4iL*Psi4rL 
      + 2*(Psi1rL*((Psi0rL*Psi3iL + Psi0iL*Psi3rL)*Psi4iL + (Psi0iL*Psi3iL - 
      Psi0rL*Psi3rL)*Psi4rL) + Psi1iL*(Psi3iL*(16*Psi1rL*Psi3rL - 
      Psi0iL*Psi4iL) + Psi0iL*Psi3rL*Psi4rL + Psi0rL*(Psi3rL*Psi4iL + 
      Psi3iL*Psi4rL))) - 15*pow(Psi2iL,4) + (-8*pow(Psi1iL,2) + 
      8*pow(Psi1rL,2))*(pow(Psi3iL,2) - pow(Psi3rL,2)) + (pow(Psi0iL,2) - 
      pow(Psi0rL,2))*pow(Psi4iL,2) + (-pow(Psi0iL,2) + 
      pow(Psi0rL,2))*pow(Psi4rL,2)) + 2*((5*(Psi1rL*Psi3iL + Psi1iL*Psi3rL) + 
      Psi0rL*Psi4iL + Psi0iL*Psi4rL)*pow(Psi2iL,3) + 
      Psi4iL*(Psi3iL*Psi3rL*pow(Psi0rL,2) + Psi0iL*(Psi4rL*pow(Psi1rL,2) + 
      Psi0rL*pow(Psi3rL,2))) + Psi2iL*((8*Psi3iL*Psi3rL + 
      3*Psi2rL*Psi4iL)*pow(Psi1rL,2) + Psi3rL*(Psi3iL*(6*Psi0rL*Psi2rL - 
      8*pow(Psi1iL,2)) + Psi1iL*(Psi0rL*Psi4rL - 15*pow(Psi2rL,2))) + 
      Psi4iL*(-(Psi0rL*Psi1iL*Psi3iL) + Psi4rL*(pow(Psi0iL,2) - 
      pow(Psi0rL,2)) - 3*(Psi2rL*pow(Psi1iL,2) + Psi0rL*pow(Psi2rL,2))) + 
      Psi1rL*(Psi0iL*Psi3rL*Psi4rL + Psi0rL*(Psi3rL*Psi4iL + Psi3iL*Psi4rL) + 
      Psi3iL*(-(Psi0iL*Psi4iL) - 15*pow(Psi2rL,2)) + Psi1iL*(6*Psi2rL*Psi4rL 
      - 8*pow(Psi3iL,2) + 8*pow(Psi3rL,2))) - Psi0iL*(Psi1iL*(Psi3rL*Psi4iL + 
      Psi3iL*Psi4rL) + 3*Psi4rL*pow(Psi2rL,2) + Psi2rL*(3*pow(Psi3iL,2) - 
      3*pow(Psi3rL,2)) + Psi0rL*(-pow(Psi4iL,2) + pow(Psi4rL,2))))));
    /* Copy local copies back to grid functions */
    curvIi[index] = curvIiL;
    curvIr[index] = curvIrL;
    curvJ1[index] = curvJ1L;
    curvJ2[index] = curvJ2L;
    curvJ3[index] = curvJ3L;
    curvJ4[index] = curvJ4L;
    curvJi[index] = curvJiL;
    curvJr[index] = curvJrL;
  }
  CCTK_ENDLOOP3(WeylScal4_invars_calc_4th);
}
extern "C" void WeylScal4_invars_calc_4th(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_WeylScal4_invars_calc_4th
  DECLARE_CCTK_ARGUMENTS_CHECKED(WeylScal4_invars_calc_4th);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WeylScal4_invars_calc_4th_Body");
  }
  if (cctk_iteration % WeylScal4_invars_calc_4th_calc_every != WeylScal4_invars_calc_4th_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "WeylScal4::curvIi_group",
    "WeylScal4::curvIr_group",
    "WeylScal4::curvJ1_group",
    "WeylScal4::curvJ2_group",
    "WeylScal4::curvJ3_group",
    "WeylScal4::curvJ4_group",
    "WeylScal4::curvJi_group",
    "WeylScal4::curvJr_group",
    "WeylScal4::Psi0i_group",
    "WeylScal4::Psi0r_group",
    "WeylScal4::Psi1i_group",
    "WeylScal4::Psi1r_group",
    "WeylScal4::Psi2i_group",
    "WeylScal4::Psi2r_group",
    "WeylScal4::Psi3i_group",
    "WeylScal4::Psi3r_group",
    "WeylScal4::Psi4i_group",
    "WeylScal4::Psi4r_group"};
  AssertGroupStorage(cctkGH, "WeylScal4_invars_calc_4th", 18, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      break;
    }
    
    case 4:
    {
      break;
    }
    
    case 6:
    {
      break;
    }
    
    case 8:
    {
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverEverything(cctkGH, WeylScal4_invars_calc_4th_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving WeylScal4_invars_calc_4th_Body");
  }
}

} // namespace WeylScal4
