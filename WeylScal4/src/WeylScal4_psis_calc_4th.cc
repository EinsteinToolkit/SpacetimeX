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

extern "C" void WeylScal4_psis_calc_4th_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_WeylScal4_psis_calc_4th_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(WeylScal4_psis_calc_4th_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % WeylScal4_psis_calc_4th_calc_every != WeylScal4_psis_calc_4th_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi0i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi0i_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi0r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi0r_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi1i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi1i_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi1r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi1r_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi2i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi2i_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi2r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi2r_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi3i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi3i_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi3r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi3r_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi4i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi4i_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi4r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi4r_group.");
  return;
}

static void WeylScal4_psis_calc_4th_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(WeylScal4_psis_calc_4th,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = kxx[index];
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = kxy[index];
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = kxz[index];
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = kyy[index];
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = kyz[index];
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = kzz[index];
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = vcoordx[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = vcoordy[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = vcoordz[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandard4th1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th22gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th33gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th23gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th33gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th12gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th13gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th23gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th22gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th12gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th13gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th23gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th11gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th33gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th13gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th11gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th12gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th13gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th23gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th11gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th22gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th12gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1kzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2kzz CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandard4th1gxx = PDstandard4th1(&gxx[index]);
        PDstandard4th2gxx = PDstandard4th2(&gxx[index]);
        PDstandard4th3gxx = PDstandard4th3(&gxx[index]);
        PDstandard4th22gxx = PDstandard4th22(&gxx[index]);
        PDstandard4th33gxx = PDstandard4th33(&gxx[index]);
        PDstandard4th23gxx = PDstandard4th23(&gxx[index]);
        PDstandard4th1gxy = PDstandard4th1(&gxy[index]);
        PDstandard4th2gxy = PDstandard4th2(&gxy[index]);
        PDstandard4th3gxy = PDstandard4th3(&gxy[index]);
        PDstandard4th33gxy = PDstandard4th33(&gxy[index]);
        PDstandard4th12gxy = PDstandard4th12(&gxy[index]);
        PDstandard4th13gxy = PDstandard4th13(&gxy[index]);
        PDstandard4th23gxy = PDstandard4th23(&gxy[index]);
        PDstandard4th1gxz = PDstandard4th1(&gxz[index]);
        PDstandard4th2gxz = PDstandard4th2(&gxz[index]);
        PDstandard4th3gxz = PDstandard4th3(&gxz[index]);
        PDstandard4th22gxz = PDstandard4th22(&gxz[index]);
        PDstandard4th12gxz = PDstandard4th12(&gxz[index]);
        PDstandard4th13gxz = PDstandard4th13(&gxz[index]);
        PDstandard4th23gxz = PDstandard4th23(&gxz[index]);
        PDstandard4th1gyy = PDstandard4th1(&gyy[index]);
        PDstandard4th2gyy = PDstandard4th2(&gyy[index]);
        PDstandard4th3gyy = PDstandard4th3(&gyy[index]);
        PDstandard4th11gyy = PDstandard4th11(&gyy[index]);
        PDstandard4th33gyy = PDstandard4th33(&gyy[index]);
        PDstandard4th13gyy = PDstandard4th13(&gyy[index]);
        PDstandard4th1gyz = PDstandard4th1(&gyz[index]);
        PDstandard4th2gyz = PDstandard4th2(&gyz[index]);
        PDstandard4th3gyz = PDstandard4th3(&gyz[index]);
        PDstandard4th11gyz = PDstandard4th11(&gyz[index]);
        PDstandard4th12gyz = PDstandard4th12(&gyz[index]);
        PDstandard4th13gyz = PDstandard4th13(&gyz[index]);
        PDstandard4th23gyz = PDstandard4th23(&gyz[index]);
        PDstandard4th1gzz = PDstandard4th1(&gzz[index]);
        PDstandard4th2gzz = PDstandard4th2(&gzz[index]);
        PDstandard4th3gzz = PDstandard4th3(&gzz[index]);
        PDstandard4th11gzz = PDstandard4th11(&gzz[index]);
        PDstandard4th22gzz = PDstandard4th22(&gzz[index]);
        PDstandard4th12gzz = PDstandard4th12(&gzz[index]);
        PDstandard4th2kxx = PDstandard4th2(&kxx[index]);
        PDstandard4th3kxx = PDstandard4th3(&kxx[index]);
        PDstandard4th1kxy = PDstandard4th1(&kxy[index]);
        PDstandard4th2kxy = PDstandard4th2(&kxy[index]);
        PDstandard4th3kxy = PDstandard4th3(&kxy[index]);
        PDstandard4th1kxz = PDstandard4th1(&kxz[index]);
        PDstandard4th2kxz = PDstandard4th2(&kxz[index]);
        PDstandard4th3kxz = PDstandard4th3(&kxz[index]);
        PDstandard4th1kyy = PDstandard4th1(&kyy[index]);
        PDstandard4th3kyy = PDstandard4th3(&kyy[index]);
        PDstandard4th1kyz = PDstandard4th1(&kyz[index]);
        PDstandard4th2kyz = PDstandard4th2(&kyz[index]);
        PDstandard4th3kyz = PDstandard4th3(&kyz[index]);
        PDstandard4th1kzz = PDstandard4th1(&kzz[index]);
        PDstandard4th2kzz = PDstandard4th2(&kzz[index]);
        break;
      }
      
      case 4:
      {
        PDstandard4th1gxx = PDstandard4th1(&gxx[index]);
        PDstandard4th2gxx = PDstandard4th2(&gxx[index]);
        PDstandard4th3gxx = PDstandard4th3(&gxx[index]);
        PDstandard4th22gxx = PDstandard4th22(&gxx[index]);
        PDstandard4th33gxx = PDstandard4th33(&gxx[index]);
        PDstandard4th23gxx = PDstandard4th23(&gxx[index]);
        PDstandard4th1gxy = PDstandard4th1(&gxy[index]);
        PDstandard4th2gxy = PDstandard4th2(&gxy[index]);
        PDstandard4th3gxy = PDstandard4th3(&gxy[index]);
        PDstandard4th33gxy = PDstandard4th33(&gxy[index]);
        PDstandard4th12gxy = PDstandard4th12(&gxy[index]);
        PDstandard4th13gxy = PDstandard4th13(&gxy[index]);
        PDstandard4th23gxy = PDstandard4th23(&gxy[index]);
        PDstandard4th1gxz = PDstandard4th1(&gxz[index]);
        PDstandard4th2gxz = PDstandard4th2(&gxz[index]);
        PDstandard4th3gxz = PDstandard4th3(&gxz[index]);
        PDstandard4th22gxz = PDstandard4th22(&gxz[index]);
        PDstandard4th12gxz = PDstandard4th12(&gxz[index]);
        PDstandard4th13gxz = PDstandard4th13(&gxz[index]);
        PDstandard4th23gxz = PDstandard4th23(&gxz[index]);
        PDstandard4th1gyy = PDstandard4th1(&gyy[index]);
        PDstandard4th2gyy = PDstandard4th2(&gyy[index]);
        PDstandard4th3gyy = PDstandard4th3(&gyy[index]);
        PDstandard4th11gyy = PDstandard4th11(&gyy[index]);
        PDstandard4th33gyy = PDstandard4th33(&gyy[index]);
        PDstandard4th13gyy = PDstandard4th13(&gyy[index]);
        PDstandard4th1gyz = PDstandard4th1(&gyz[index]);
        PDstandard4th2gyz = PDstandard4th2(&gyz[index]);
        PDstandard4th3gyz = PDstandard4th3(&gyz[index]);
        PDstandard4th11gyz = PDstandard4th11(&gyz[index]);
        PDstandard4th12gyz = PDstandard4th12(&gyz[index]);
        PDstandard4th13gyz = PDstandard4th13(&gyz[index]);
        PDstandard4th23gyz = PDstandard4th23(&gyz[index]);
        PDstandard4th1gzz = PDstandard4th1(&gzz[index]);
        PDstandard4th2gzz = PDstandard4th2(&gzz[index]);
        PDstandard4th3gzz = PDstandard4th3(&gzz[index]);
        PDstandard4th11gzz = PDstandard4th11(&gzz[index]);
        PDstandard4th22gzz = PDstandard4th22(&gzz[index]);
        PDstandard4th12gzz = PDstandard4th12(&gzz[index]);
        PDstandard4th2kxx = PDstandard4th2(&kxx[index]);
        PDstandard4th3kxx = PDstandard4th3(&kxx[index]);
        PDstandard4th1kxy = PDstandard4th1(&kxy[index]);
        PDstandard4th2kxy = PDstandard4th2(&kxy[index]);
        PDstandard4th3kxy = PDstandard4th3(&kxy[index]);
        PDstandard4th1kxz = PDstandard4th1(&kxz[index]);
        PDstandard4th2kxz = PDstandard4th2(&kxz[index]);
        PDstandard4th3kxz = PDstandard4th3(&kxz[index]);
        PDstandard4th1kyy = PDstandard4th1(&kyy[index]);
        PDstandard4th3kyy = PDstandard4th3(&kyy[index]);
        PDstandard4th1kyz = PDstandard4th1(&kyz[index]);
        PDstandard4th2kyz = PDstandard4th2(&kyz[index]);
        PDstandard4th3kyz = PDstandard4th3(&kyz[index]);
        PDstandard4th1kzz = PDstandard4th1(&kzz[index]);
        PDstandard4th2kzz = PDstandard4th2(&kzz[index]);
        break;
      }
      
      case 6:
      {
        PDstandard4th1gxx = PDstandard4th1(&gxx[index]);
        PDstandard4th2gxx = PDstandard4th2(&gxx[index]);
        PDstandard4th3gxx = PDstandard4th3(&gxx[index]);
        PDstandard4th22gxx = PDstandard4th22(&gxx[index]);
        PDstandard4th33gxx = PDstandard4th33(&gxx[index]);
        PDstandard4th23gxx = PDstandard4th23(&gxx[index]);
        PDstandard4th1gxy = PDstandard4th1(&gxy[index]);
        PDstandard4th2gxy = PDstandard4th2(&gxy[index]);
        PDstandard4th3gxy = PDstandard4th3(&gxy[index]);
        PDstandard4th33gxy = PDstandard4th33(&gxy[index]);
        PDstandard4th12gxy = PDstandard4th12(&gxy[index]);
        PDstandard4th13gxy = PDstandard4th13(&gxy[index]);
        PDstandard4th23gxy = PDstandard4th23(&gxy[index]);
        PDstandard4th1gxz = PDstandard4th1(&gxz[index]);
        PDstandard4th2gxz = PDstandard4th2(&gxz[index]);
        PDstandard4th3gxz = PDstandard4th3(&gxz[index]);
        PDstandard4th22gxz = PDstandard4th22(&gxz[index]);
        PDstandard4th12gxz = PDstandard4th12(&gxz[index]);
        PDstandard4th13gxz = PDstandard4th13(&gxz[index]);
        PDstandard4th23gxz = PDstandard4th23(&gxz[index]);
        PDstandard4th1gyy = PDstandard4th1(&gyy[index]);
        PDstandard4th2gyy = PDstandard4th2(&gyy[index]);
        PDstandard4th3gyy = PDstandard4th3(&gyy[index]);
        PDstandard4th11gyy = PDstandard4th11(&gyy[index]);
        PDstandard4th33gyy = PDstandard4th33(&gyy[index]);
        PDstandard4th13gyy = PDstandard4th13(&gyy[index]);
        PDstandard4th1gyz = PDstandard4th1(&gyz[index]);
        PDstandard4th2gyz = PDstandard4th2(&gyz[index]);
        PDstandard4th3gyz = PDstandard4th3(&gyz[index]);
        PDstandard4th11gyz = PDstandard4th11(&gyz[index]);
        PDstandard4th12gyz = PDstandard4th12(&gyz[index]);
        PDstandard4th13gyz = PDstandard4th13(&gyz[index]);
        PDstandard4th23gyz = PDstandard4th23(&gyz[index]);
        PDstandard4th1gzz = PDstandard4th1(&gzz[index]);
        PDstandard4th2gzz = PDstandard4th2(&gzz[index]);
        PDstandard4th3gzz = PDstandard4th3(&gzz[index]);
        PDstandard4th11gzz = PDstandard4th11(&gzz[index]);
        PDstandard4th22gzz = PDstandard4th22(&gzz[index]);
        PDstandard4th12gzz = PDstandard4th12(&gzz[index]);
        PDstandard4th2kxx = PDstandard4th2(&kxx[index]);
        PDstandard4th3kxx = PDstandard4th3(&kxx[index]);
        PDstandard4th1kxy = PDstandard4th1(&kxy[index]);
        PDstandard4th2kxy = PDstandard4th2(&kxy[index]);
        PDstandard4th3kxy = PDstandard4th3(&kxy[index]);
        PDstandard4th1kxz = PDstandard4th1(&kxz[index]);
        PDstandard4th2kxz = PDstandard4th2(&kxz[index]);
        PDstandard4th3kxz = PDstandard4th3(&kxz[index]);
        PDstandard4th1kyy = PDstandard4th1(&kyy[index]);
        PDstandard4th3kyy = PDstandard4th3(&kyy[index]);
        PDstandard4th1kyz = PDstandard4th1(&kyz[index]);
        PDstandard4th2kyz = PDstandard4th2(&kyz[index]);
        PDstandard4th3kyz = PDstandard4th3(&kyz[index]);
        PDstandard4th1kzz = PDstandard4th1(&kzz[index]);
        PDstandard4th2kzz = PDstandard4th2(&kzz[index]);
        break;
      }
      
      case 8:
      {
        PDstandard4th1gxx = PDstandard4th1(&gxx[index]);
        PDstandard4th2gxx = PDstandard4th2(&gxx[index]);
        PDstandard4th3gxx = PDstandard4th3(&gxx[index]);
        PDstandard4th22gxx = PDstandard4th22(&gxx[index]);
        PDstandard4th33gxx = PDstandard4th33(&gxx[index]);
        PDstandard4th23gxx = PDstandard4th23(&gxx[index]);
        PDstandard4th1gxy = PDstandard4th1(&gxy[index]);
        PDstandard4th2gxy = PDstandard4th2(&gxy[index]);
        PDstandard4th3gxy = PDstandard4th3(&gxy[index]);
        PDstandard4th33gxy = PDstandard4th33(&gxy[index]);
        PDstandard4th12gxy = PDstandard4th12(&gxy[index]);
        PDstandard4th13gxy = PDstandard4th13(&gxy[index]);
        PDstandard4th23gxy = PDstandard4th23(&gxy[index]);
        PDstandard4th1gxz = PDstandard4th1(&gxz[index]);
        PDstandard4th2gxz = PDstandard4th2(&gxz[index]);
        PDstandard4th3gxz = PDstandard4th3(&gxz[index]);
        PDstandard4th22gxz = PDstandard4th22(&gxz[index]);
        PDstandard4th12gxz = PDstandard4th12(&gxz[index]);
        PDstandard4th13gxz = PDstandard4th13(&gxz[index]);
        PDstandard4th23gxz = PDstandard4th23(&gxz[index]);
        PDstandard4th1gyy = PDstandard4th1(&gyy[index]);
        PDstandard4th2gyy = PDstandard4th2(&gyy[index]);
        PDstandard4th3gyy = PDstandard4th3(&gyy[index]);
        PDstandard4th11gyy = PDstandard4th11(&gyy[index]);
        PDstandard4th33gyy = PDstandard4th33(&gyy[index]);
        PDstandard4th13gyy = PDstandard4th13(&gyy[index]);
        PDstandard4th1gyz = PDstandard4th1(&gyz[index]);
        PDstandard4th2gyz = PDstandard4th2(&gyz[index]);
        PDstandard4th3gyz = PDstandard4th3(&gyz[index]);
        PDstandard4th11gyz = PDstandard4th11(&gyz[index]);
        PDstandard4th12gyz = PDstandard4th12(&gyz[index]);
        PDstandard4th13gyz = PDstandard4th13(&gyz[index]);
        PDstandard4th23gyz = PDstandard4th23(&gyz[index]);
        PDstandard4th1gzz = PDstandard4th1(&gzz[index]);
        PDstandard4th2gzz = PDstandard4th2(&gzz[index]);
        PDstandard4th3gzz = PDstandard4th3(&gzz[index]);
        PDstandard4th11gzz = PDstandard4th11(&gzz[index]);
        PDstandard4th22gzz = PDstandard4th22(&gzz[index]);
        PDstandard4th12gzz = PDstandard4th12(&gzz[index]);
        PDstandard4th2kxx = PDstandard4th2(&kxx[index]);
        PDstandard4th3kxx = PDstandard4th3(&kxx[index]);
        PDstandard4th1kxy = PDstandard4th1(&kxy[index]);
        PDstandard4th2kxy = PDstandard4th2(&kxy[index]);
        PDstandard4th3kxy = PDstandard4th3(&kxy[index]);
        PDstandard4th1kxz = PDstandard4th1(&kxz[index]);
        PDstandard4th2kxz = PDstandard4th2(&kxz[index]);
        PDstandard4th3kxz = PDstandard4th3(&kxz[index]);
        PDstandard4th1kyy = PDstandard4th1(&kyy[index]);
        PDstandard4th3kyy = PDstandard4th3(&kyy[index]);
        PDstandard4th1kyz = PDstandard4th1(&kyz[index]);
        PDstandard4th2kyz = PDstandard4th2(&kyz[index]);
        PDstandard4th3kyz = PDstandard4th3(&kyz[index]);
        PDstandard4th1kzz = PDstandard4th1(&kzz[index]);
        PDstandard4th2kzz = PDstandard4th2(&kzz[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*gxyL*gxzL*gyzL - 
      gzzL*pow(gxyL,2) + gyyL*(gxxL*gzzL - pow(gxzL,2)) - gxxL*pow(gyzL,2);
    
    CCTK_REAL invdetg CCTK_ATTRIBUTE_UNUSED = pow(detg,-1);
    
    CCTK_REAL gInv11 CCTK_ATTRIBUTE_UNUSED = invdetg*(gyyL*gzzL - 
      pow(gyzL,2));
    
    CCTK_REAL gInv12 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*invdetg;
    
    CCTK_REAL gInv13 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*invdetg;
    
    CCTK_REAL gInv21 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*invdetg;
    
    CCTK_REAL gInv22 CCTK_ATTRIBUTE_UNUSED = invdetg*(gxxL*gzzL - 
      pow(gxzL,2));
    
    CCTK_REAL gInv23 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*invdetg;
    
    CCTK_REAL gInv31 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*invdetg;
    
    CCTK_REAL gInv32 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*invdetg;
    
    CCTK_REAL gInv33 CCTK_ATTRIBUTE_UNUSED = invdetg*(gxxL*gyyL - 
      pow(gxyL,2));
    
    CCTK_REAL gamma111 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv11*PDstandard4th1gxx + 2*(gInv12*PDstandard4th1gxy + 
      gInv13*PDstandard4th1gxz) - gInv12*PDstandard4th2gxx - 
      gInv13*PDstandard4th3gxx);
    
    CCTK_REAL gamma211 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv21*PDstandard4th1gxx + 2*(gInv22*PDstandard4th1gxy + 
      gInv23*PDstandard4th1gxz) - gInv22*PDstandard4th2gxx - 
      gInv23*PDstandard4th3gxx);
    
    CCTK_REAL gamma311 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv31*PDstandard4th1gxx + 2*(gInv32*PDstandard4th1gxy + 
      gInv33*PDstandard4th1gxz) - gInv32*PDstandard4th2gxx - 
      gInv33*PDstandard4th3gxx);
    
    CCTK_REAL gamma121 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv12*PDstandard4th1gyy + gInv11*PDstandard4th2gxx + 
      gInv13*(PDstandard4th1gyz + PDstandard4th2gxz - PDstandard4th3gxy));
    
    CCTK_REAL gamma221 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv22*PDstandard4th1gyy + gInv21*PDstandard4th2gxx + 
      gInv23*(PDstandard4th1gyz + PDstandard4th2gxz - PDstandard4th3gxy));
    
    CCTK_REAL gamma321 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv32*PDstandard4th1gyy + gInv31*PDstandard4th2gxx + 
      gInv33*(PDstandard4th1gyz + PDstandard4th2gxz - PDstandard4th3gxy));
    
    CCTK_REAL gamma131 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv13*PDstandard4th1gzz + gInv11*PDstandard4th3gxx + 
      gInv12*(PDstandard4th1gyz - PDstandard4th2gxz + PDstandard4th3gxy));
    
    CCTK_REAL gamma231 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv23*PDstandard4th1gzz + gInv21*PDstandard4th3gxx + 
      gInv22*(PDstandard4th1gyz - PDstandard4th2gxz + PDstandard4th3gxy));
    
    CCTK_REAL gamma331 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv33*PDstandard4th1gzz + gInv31*PDstandard4th3gxx + 
      gInv32*(PDstandard4th1gyz - PDstandard4th2gxz + PDstandard4th3gxy));
    
    CCTK_REAL gamma122 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv11*(-PDstandard4th1gyy + 2*PDstandard4th2gxy) + 
      gInv12*PDstandard4th2gyy + gInv13*(2*PDstandard4th2gyz - 
      PDstandard4th3gyy));
    
    CCTK_REAL gamma222 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv21*(-PDstandard4th1gyy + 2*PDstandard4th2gxy) + 
      gInv22*PDstandard4th2gyy + gInv23*(2*PDstandard4th2gyz - 
      PDstandard4th3gyy));
    
    CCTK_REAL gamma322 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv31*(-PDstandard4th1gyy + 2*PDstandard4th2gxy) + 
      gInv32*PDstandard4th2gyy + gInv33*(2*PDstandard4th2gyz - 
      PDstandard4th3gyy));
    
    CCTK_REAL gamma132 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv13*PDstandard4th2gzz + gInv11*(-PDstandard4th1gyz + 
      PDstandard4th2gxz + PDstandard4th3gxy) + gInv12*PDstandard4th3gyy);
    
    CCTK_REAL gamma232 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv23*PDstandard4th2gzz + gInv21*(-PDstandard4th1gyz + 
      PDstandard4th2gxz + PDstandard4th3gxy) + gInv22*PDstandard4th3gyy);
    
    CCTK_REAL gamma332 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv33*PDstandard4th2gzz + gInv31*(-PDstandard4th1gyz + 
      PDstandard4th2gxz + PDstandard4th3gxy) + gInv32*PDstandard4th3gyy);
    
    CCTK_REAL gamma133 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv11*(-PDstandard4th1gzz + 2*PDstandard4th3gxz) + 
      gInv12*(-PDstandard4th2gzz + 2*PDstandard4th3gyz) + 
      gInv13*PDstandard4th3gzz);
    
    CCTK_REAL gamma233 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv21*(-PDstandard4th1gzz + 2*PDstandard4th3gxz) + 
      gInv22*(-PDstandard4th2gzz + 2*PDstandard4th3gyz) + 
      gInv23*PDstandard4th3gzz);
    
    CCTK_REAL gamma333 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv31*(-PDstandard4th1gzz + 2*PDstandard4th3gxz) + 
      gInv32*(-PDstandard4th2gzz + 2*PDstandard4th3gyz) + 
      gInv33*PDstandard4th3gzz);
    
    CCTK_REAL xmoved CCTK_ATTRIBUTE_UNUSED = xL - xorig;
    
    CCTK_REAL ymoved CCTK_ATTRIBUTE_UNUSED = yL - yorig;
    
    CCTK_REAL zmoved CCTK_ATTRIBUTE_UNUSED = zL - zorig;
    
    CCTK_REAL va1 CCTK_ATTRIBUTE_UNUSED = -ymoved;
    
    CCTK_REAL va2 CCTK_ATTRIBUTE_UNUSED = offset + xmoved;
    
    CCTK_REAL va3 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL vb1 CCTK_ATTRIBUTE_UNUSED = offset + xmoved;
    
    CCTK_REAL vb2 CCTK_ATTRIBUTE_UNUSED = ymoved;
    
    CCTK_REAL vb3 CCTK_ATTRIBUTE_UNUSED = zmoved;
    
    CCTK_REAL vc1 CCTK_ATTRIBUTE_UNUSED = ((-(gInv13*va2) + 
      gInv12*va3)*vb1 + (gInv13*va1 - gInv11*va3)*vb2 + (-(gInv12*va1) + 
      gInv11*va2)*vb3)*pow(detg,0.5);
    
    CCTK_REAL vc2 CCTK_ATTRIBUTE_UNUSED = ((-(gInv23*va2) + 
      gInv22*va3)*vb1 + (gInv23*va1 - gInv21*va3)*vb2 + (-(gInv22*va1) + 
      gInv21*va2)*vb3)*pow(detg,0.5);
    
    CCTK_REAL vc3 CCTK_ATTRIBUTE_UNUSED = ((-(gInv33*va2) + 
      gInv32*va3)*vb1 + (gInv33*va1 - gInv31*va3)*vb2 + (-(gInv32*va1) + 
      gInv31*va2)*vb3)*pow(detg,0.5);
    
    CCTK_REAL wa1 CCTK_ATTRIBUTE_UNUSED = va1;
    
    CCTK_REAL wa2 CCTK_ATTRIBUTE_UNUSED = va2;
    
    CCTK_REAL wa3 CCTK_ATTRIBUTE_UNUSED = va3;
    
    CCTK_REAL omega11 CCTK_ATTRIBUTE_UNUSED = 2*(gyzL*wa2*wa3 + 
      wa1*(gxyL*wa2 + gxzL*wa3)) + gxxL*pow(wa1,2) + gyyL*pow(wa2,2) + 
      gzzL*pow(wa3,2);
    
    CCTK_REAL ea1 CCTK_ATTRIBUTE_UNUSED = wa1*pow(omega11,-0.5);
    
    CCTK_REAL ea2 CCTK_ATTRIBUTE_UNUSED = wa2*pow(omega11,-0.5);
    
    CCTK_REAL ea3 CCTK_ATTRIBUTE_UNUSED = wa3*pow(omega11,-0.5);
    
    CCTK_REAL omega12 CCTK_ATTRIBUTE_UNUSED = ea1*(gxxL*vb1 + gxyL*vb2 + 
      gxzL*vb3) + ea2*(gxyL*vb1 + gyyL*vb2 + gyzL*vb3) + ea3*(gxzL*vb1 + 
      gyzL*vb2 + gzzL*vb3);
    
    CCTK_REAL wb1 CCTK_ATTRIBUTE_UNUSED = -(ea1*omega12) + vb1;
    
    CCTK_REAL wb2 CCTK_ATTRIBUTE_UNUSED = -(ea2*omega12) + vb2;
    
    CCTK_REAL wb3 CCTK_ATTRIBUTE_UNUSED = -(ea3*omega12) + vb3;
    
    CCTK_REAL omega22 CCTK_ATTRIBUTE_UNUSED = 2*(gyzL*wb2*wb3 + 
      wb1*(gxyL*wb2 + gxzL*wb3)) + gxxL*pow(wb1,2) + gyyL*pow(wb2,2) + 
      gzzL*pow(wb3,2);
    
    CCTK_REAL eb1 CCTK_ATTRIBUTE_UNUSED = wb1*pow(omega22,-0.5);
    
    CCTK_REAL eb2 CCTK_ATTRIBUTE_UNUSED = wb2*pow(omega22,-0.5);
    
    CCTK_REAL eb3 CCTK_ATTRIBUTE_UNUSED = wb3*pow(omega22,-0.5);
    
    CCTK_REAL omega13 CCTK_ATTRIBUTE_UNUSED = ea1*(gxxL*vc1 + gxyL*vc2 + 
      gxzL*vc3) + ea2*(gxyL*vc1 + gyyL*vc2 + gyzL*vc3) + ea3*(gxzL*vc1 + 
      gyzL*vc2 + gzzL*vc3);
    
    CCTK_REAL omega23 CCTK_ATTRIBUTE_UNUSED = eb1*(gxxL*vc1 + gxyL*vc2 + 
      gxzL*vc3) + eb2*(gxyL*vc1 + gyyL*vc2 + gyzL*vc3) + eb3*(gxzL*vc1 + 
      gyzL*vc2 + gzzL*vc3);
    
    CCTK_REAL wc1 CCTK_ATTRIBUTE_UNUSED = -(ea1*omega13) - eb1*omega23 + 
      vc1;
    
    CCTK_REAL wc2 CCTK_ATTRIBUTE_UNUSED = -(ea2*omega13) - eb2*omega23 + 
      vc2;
    
    CCTK_REAL wc3 CCTK_ATTRIBUTE_UNUSED = -(ea3*omega13) - eb3*omega23 + 
      vc3;
    
    CCTK_REAL omega33 CCTK_ATTRIBUTE_UNUSED = 2*(gyzL*wc2*wc3 + 
      wc1*(gxyL*wc2 + gxzL*wc3)) + gxxL*pow(wc1,2) + gyyL*pow(wc2,2) + 
      gzzL*pow(wc3,2);
    
    CCTK_REAL ec1 CCTK_ATTRIBUTE_UNUSED = wc1*pow(omega33,-0.5);
    
    CCTK_REAL ec2 CCTK_ATTRIBUTE_UNUSED = wc2*pow(omega33,-0.5);
    
    CCTK_REAL ec3 CCTK_ATTRIBUTE_UNUSED = wc3*pow(omega33,-0.5);
    
    CCTK_REAL isqrt2 CCTK_ATTRIBUTE_UNUSED = 0.707106781186547524;
    
    CCTK_REAL ltet1 CCTK_ATTRIBUTE_UNUSED = eb1*isqrt2;
    
    CCTK_REAL ltet2 CCTK_ATTRIBUTE_UNUSED = eb2*isqrt2;
    
    CCTK_REAL ltet3 CCTK_ATTRIBUTE_UNUSED = eb3*isqrt2;
    
    CCTK_REAL n1 CCTK_ATTRIBUTE_UNUSED = -(eb1*isqrt2);
    
    CCTK_REAL n2 CCTK_ATTRIBUTE_UNUSED = -(eb2*isqrt2);
    
    CCTK_REAL n3 CCTK_ATTRIBUTE_UNUSED = -(eb3*isqrt2);
    
    CCTK_REAL rm1 CCTK_ATTRIBUTE_UNUSED = ec1*isqrt2;
    
    CCTK_REAL rm2 CCTK_ATTRIBUTE_UNUSED = ec2*isqrt2;
    
    CCTK_REAL rm3 CCTK_ATTRIBUTE_UNUSED = ec3*isqrt2;
    
    CCTK_REAL im1 CCTK_ATTRIBUTE_UNUSED = ea1*isqrt2;
    
    CCTK_REAL im2 CCTK_ATTRIBUTE_UNUSED = ea2*isqrt2;
    
    CCTK_REAL im3 CCTK_ATTRIBUTE_UNUSED = ea3*isqrt2;
    
    CCTK_REAL rmbar1 CCTK_ATTRIBUTE_UNUSED = ec1*isqrt2;
    
    CCTK_REAL rmbar2 CCTK_ATTRIBUTE_UNUSED = ec2*isqrt2;
    
    CCTK_REAL rmbar3 CCTK_ATTRIBUTE_UNUSED = ec3*isqrt2;
    
    CCTK_REAL imbar1 CCTK_ATTRIBUTE_UNUSED = -(ea1*isqrt2);
    
    CCTK_REAL imbar2 CCTK_ATTRIBUTE_UNUSED = -(ea2*isqrt2);
    
    CCTK_REAL imbar3 CCTK_ATTRIBUTE_UNUSED = -(ea3*isqrt2);
    
    CCTK_REAL nn CCTK_ATTRIBUTE_UNUSED = isqrt2;
    
    CCTK_REAL R1212 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(-2*(gamma122*(gxxL*gamma111 + gxyL*gamma211 + gxzL*gamma311) + 
      gamma222*(gxyL*gamma111 + gyyL*gamma211 + gyzL*gamma311) + 
      (gxzL*gamma111 + gyzL*gamma211 + gzzL*gamma311)*gamma322) - 
      PDstandard4th11gyy + 2*(gamma121*(gxxL*gamma121 + gxyL*gamma221 + 
      gxzL*gamma321) + gamma221*(gxyL*gamma121 + gyyL*gamma221 + 
      gyzL*gamma321) + gamma321*(gxzL*gamma121 + gyzL*gamma221 + 
      gzzL*gamma321) + PDstandard4th12gxy) - PDstandard4th22gxx);
    
    CCTK_REAL R1213 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(2*(gamma121*(gxxL*gamma131 + gxyL*gamma231 + gxzL*gamma331) + 
      gamma221*(gxyL*gamma131 + gyyL*gamma231 + gyzL*gamma331) + 
      gamma321*(gxzL*gamma131 + gyzL*gamma231 + gzzL*gamma331)) - 
      2*(gamma132*(gxxL*gamma111 + gxyL*gamma211 + gxzL*gamma311) + 
      gamma232*(gxyL*gamma111 + gyyL*gamma211 + gyzL*gamma311) + 
      (gxzL*gamma111 + gyzL*gamma211 + gzzL*gamma311)*gamma332) - 
      PDstandard4th11gyz + PDstandard4th12gxz + PDstandard4th13gxy - 
      PDstandard4th23gxx);
    
    CCTK_REAL R1223 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(2*(gamma122*(gxxL*gamma131 + gxyL*gamma231 + gxzL*gamma331) + 
      gamma222*(gxyL*gamma131 + gyyL*gamma231 + gyzL*gamma331) + 
      gamma322*(gxzL*gamma131 + gyzL*gamma231 + gzzL*gamma331)) - 
      2*(gamma132*(gxxL*gamma121 + gxyL*gamma221 + gxzL*gamma321) + 
      gamma232*(gxyL*gamma121 + gyyL*gamma221 + gyzL*gamma321) + 
      (gxzL*gamma121 + gyzL*gamma221 + gzzL*gamma321)*gamma332) - 
      PDstandard4th12gyz + PDstandard4th13gyy + PDstandard4th22gxz - 
      PDstandard4th23gxy);
    
    CCTK_REAL R1313 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(-2*(gamma133*(gxxL*gamma111 + gxyL*gamma211 + gxzL*gamma311) + 
      gamma233*(gxyL*gamma111 + gyyL*gamma211 + gyzL*gamma311) + 
      (gxzL*gamma111 + gyzL*gamma211 + gzzL*gamma311)*gamma333) - 
      PDstandard4th11gzz + 2*(gamma131*(gxxL*gamma131 + gxyL*gamma231 + 
      gxzL*gamma331) + gamma231*(gxyL*gamma131 + gyyL*gamma231 + 
      gyzL*gamma331) + gamma331*(gxzL*gamma131 + gyzL*gamma231 + 
      gzzL*gamma331) + PDstandard4th13gxz) - PDstandard4th33gxx);
    
    CCTK_REAL R1323 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(2*(gamma132*(gxxL*gamma131 + gxyL*gamma231 + gxzL*gamma331) + 
      gamma232*(gxyL*gamma131 + gyyL*gamma231 + gyzL*gamma331) + 
      (gxzL*gamma131 + gyzL*gamma231 + gzzL*gamma331)*gamma332) - 
      2*(gamma133*(gxxL*gamma121 + gxyL*gamma221 + gxzL*gamma321) + 
      gamma233*(gxyL*gamma121 + gyyL*gamma221 + gyzL*gamma321) + 
      (gxzL*gamma121 + gyzL*gamma221 + gzzL*gamma321)*gamma333) - 
      PDstandard4th12gzz + PDstandard4th13gyz + PDstandard4th23gxz - 
      PDstandard4th33gxy);
    
    CCTK_REAL R2323 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(-2*(gamma133*(gxxL*gamma122 + gxyL*gamma222 + gxzL*gamma322) + 
      gamma233*(gxyL*gamma122 + gyyL*gamma222 + gyzL*gamma322) + 
      (gxzL*gamma122 + gyzL*gamma222 + gzzL*gamma322)*gamma333) - 
      PDstandard4th22gzz + 2*(gamma132*(gxxL*gamma132 + gxyL*gamma232 + 
      gxzL*gamma332) + gamma232*(gxyL*gamma132 + gyyL*gamma232 + 
      gyzL*gamma332) + gamma332*(gxzL*gamma132 + gyzL*gamma232 + 
      gzzL*gamma332) + PDstandard4th23gyz) - PDstandard4th33gyy);
    
    CCTK_REAL R4p1212 CCTK_ATTRIBUTE_UNUSED = kxxL*kyyL + R1212 - 
      pow(kxyL,2);
    
    CCTK_REAL R4p1213 CCTK_ATTRIBUTE_UNUSED = -(kxyL*kxzL) + kxxL*kyzL + 
      R1213;
    
    CCTK_REAL R4p1223 CCTK_ATTRIBUTE_UNUSED = -(kxzL*kyyL) + kxyL*kyzL + 
      R1223;
    
    CCTK_REAL R4p1313 CCTK_ATTRIBUTE_UNUSED = kxxL*kzzL + R1313 - 
      pow(kxzL,2);
    
    CCTK_REAL R4p1323 CCTK_ATTRIBUTE_UNUSED = -(kxzL*kyzL) + kxyL*kzzL + 
      R1323;
    
    CCTK_REAL R4p2323 CCTK_ATTRIBUTE_UNUSED = kyyL*kzzL + R2323 - 
      pow(kyzL,2);
    
    CCTK_REAL Ro111 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro112 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma121 - kyyL*gamma211 
      + kxyL*(-gamma111 + gamma221) - kyzL*gamma311 + kxzL*gamma321 + 
      PDstandard4th1kxy - PDstandard4th2kxx;
    
    CCTK_REAL Ro113 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma131 - kyzL*gamma211 
      + kxyL*gamma231 - kzzL*gamma311 + kxzL*(-gamma111 + gamma331) + 
      PDstandard4th1kxz - PDstandard4th3kxx;
    
    CCTK_REAL Ro121 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma121) + 
      kyyL*gamma211 + kxyL*(gamma111 - gamma221) + kyzL*gamma311 - 
      kxzL*gamma321 - PDstandard4th1kxy + PDstandard4th2kxx;
    
    CCTK_REAL Ro122 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro123 CCTK_ATTRIBUTE_UNUSED = -(kxzL*gamma121) + 
      kxyL*gamma131 + kyyL*gamma231 - kzzL*gamma321 + kyzL*(-gamma221 + 
      gamma331) + PDstandard4th2kxz - PDstandard4th3kxy;
    
    CCTK_REAL Ro131 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma131) + 
      kyzL*gamma211 - kxyL*gamma231 + kzzL*gamma311 + kxzL*(gamma111 - 
      gamma331) - PDstandard4th1kxz + PDstandard4th3kxx;
    
    CCTK_REAL Ro132 CCTK_ATTRIBUTE_UNUSED = kxzL*gamma121 - kxyL*gamma131 
      - kyyL*gamma231 + kzzL*gamma321 + kyzL*(gamma221 - gamma331) - 
      PDstandard4th2kxz + PDstandard4th3kxy;
    
    CCTK_REAL Ro133 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro211 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro212 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma122 - kyyL*gamma221 
      + kxyL*(-gamma121 + gamma222) - kyzL*gamma321 + kxzL*gamma322 + 
      PDstandard4th1kyy - PDstandard4th2kxy;
    
    CCTK_REAL Ro213 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma132 - kyzL*gamma221 
      + kxyL*gamma232 - kzzL*gamma321 + kxzL*(-gamma121 + gamma332) + 
      PDstandard4th1kyz - PDstandard4th3kxy;
    
    CCTK_REAL Ro221 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma122) + 
      kyyL*gamma221 + kxyL*(gamma121 - gamma222) + kyzL*gamma321 - 
      kxzL*gamma322 - PDstandard4th1kyy + PDstandard4th2kxy;
    
    CCTK_REAL Ro222 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro223 CCTK_ATTRIBUTE_UNUSED = -(kxzL*gamma122) + 
      kxyL*gamma132 + kyyL*gamma232 - kzzL*gamma322 + kyzL*(-gamma222 + 
      gamma332) + PDstandard4th2kyz - PDstandard4th3kyy;
    
    CCTK_REAL Ro231 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma132) + 
      kyzL*gamma221 - kxyL*gamma232 + kzzL*gamma321 + kxzL*(gamma121 - 
      gamma332) - PDstandard4th1kyz + PDstandard4th3kxy;
    
    CCTK_REAL Ro232 CCTK_ATTRIBUTE_UNUSED = kxzL*gamma122 - kxyL*gamma132 
      - kyyL*gamma232 + kzzL*gamma322 + kyzL*(gamma222 - gamma332) - 
      PDstandard4th2kyz + PDstandard4th3kyy;
    
    CCTK_REAL Ro233 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro311 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro312 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma132 - kyyL*gamma231 
      + kxyL*(-gamma131 + gamma232) - kyzL*gamma331 + kxzL*gamma332 + 
      PDstandard4th1kyz - PDstandard4th2kxz;
    
    CCTK_REAL Ro313 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma133 - kyzL*gamma231 
      + kxyL*gamma233 - kzzL*gamma331 + kxzL*(-gamma131 + gamma333) + 
      PDstandard4th1kzz - PDstandard4th3kxz;
    
    CCTK_REAL Ro321 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma132) + 
      kyyL*gamma231 + kxyL*(gamma131 - gamma232) + kyzL*gamma331 - 
      kxzL*gamma332 - PDstandard4th1kyz + PDstandard4th2kxz;
    
    CCTK_REAL Ro322 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro323 CCTK_ATTRIBUTE_UNUSED = -(kxzL*gamma132) + 
      kxyL*gamma133 + kyyL*gamma233 - kzzL*gamma332 + kyzL*(-gamma232 + 
      gamma333) + PDstandard4th2kzz - PDstandard4th3kyz;
    
    CCTK_REAL Ro331 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma133) + 
      kyzL*gamma231 - kxyL*gamma233 + kzzL*gamma331 + kxzL*(gamma131 - 
      gamma333) - PDstandard4th1kzz + PDstandard4th3kxz;
    
    CCTK_REAL Ro332 CCTK_ATTRIBUTE_UNUSED = kxzL*gamma132 - kxyL*gamma133 
      - kyyL*gamma233 + kzzL*gamma332 + kyzL*(gamma232 - gamma333) - 
      PDstandard4th2kzz + PDstandard4th3kyz;
    
    CCTK_REAL Ro333 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Rojo11 CCTK_ATTRIBUTE_UNUSED = (gInv23 + 
      gInv32)*(-(kxyL*kxzL) + kxxL*kyzL + R1213) + gInv22*(kxxL*kyyL + R1212 
      - pow(kxyL,2)) + gInv33*(kxxL*kzzL + R1313 - pow(kxzL,2));
    
    CCTK_REAL Rojo12 CCTK_ATTRIBUTE_UNUSED = (kxyL*kxzL - 
      kxxL*kyzL)*gInv13 + (-(kxzL*kyyL) + kxyL*kyzL)*gInv32 - gInv21*R1212 - 
      gInv31*R1213 + gInv23*R1223 + gInv33*(-(kxzL*kyzL) + kxyL*kzzL + R1323) 
      + gInv12*(-(kxxL*kyyL) + pow(kxyL,2));
    
    CCTK_REAL Rojo13 CCTK_ATTRIBUTE_UNUSED = (kxyL*kxzL - 
      kxxL*kyzL)*gInv12 + (kxzL*kyzL - kxyL*kzzL)*gInv23 - gInv21*R1213 + 
      gInv22*(kxzL*kyyL - kxyL*kyzL - R1223) - gInv31*R1313 - gInv32*R1323 + 
      gInv13*(-(kxxL*kzzL) + pow(kxzL,2));
    
    CCTK_REAL Rojo21 CCTK_ATTRIBUTE_UNUSED = (-(kxzL*kyyL) + 
      kxyL*kyzL)*gInv23 + (kxyL*kxzL - kxxL*kyzL)*gInv31 - gInv12*R1212 - 
      gInv13*R1213 + gInv32*R1223 + gInv33*(-(kxzL*kyzL) + kxyL*kzzL + R1323) 
      + gInv21*(-(kxxL*kyyL) + pow(kxyL,2));
    
    CCTK_REAL Rojo22 CCTK_ATTRIBUTE_UNUSED = (gInv13 + gInv31)*(kxzL*kyyL 
      - kxyL*kyzL - R1223) + gInv11*(kxxL*kyyL + R1212 - pow(kxyL,2)) + 
      gInv33*(kyyL*kzzL + R2323 - pow(kyzL,2));
    
    CCTK_REAL Rojo23 CCTK_ATTRIBUTE_UNUSED = (kxzL*kyzL - 
      kxyL*kzzL)*gInv13 + (-(kxzL*kyyL) + kxyL*kyzL)*gInv21 + 
      gInv11*(-(kxyL*kxzL) + kxxL*kyzL + R1213) + gInv12*R1223 - gInv31*R1323 
      - gInv32*R2323 + gInv23*(-(kyyL*kzzL) + pow(kyzL,2));
    
    CCTK_REAL Rojo31 CCTK_ATTRIBUTE_UNUSED = (kxyL*kxzL - 
      kxxL*kyzL)*gInv21 + (kxzL*kyzL - kxyL*kzzL)*gInv32 - gInv12*R1213 + 
      gInv22*(kxzL*kyyL - kxyL*kyzL - R1223) - gInv13*R1313 - gInv23*R1323 + 
      gInv31*(-(kxxL*kzzL) + pow(kxzL,2));
    
    CCTK_REAL Rojo32 CCTK_ATTRIBUTE_UNUSED = (-(kxzL*kyyL) + 
      kxyL*kyzL)*gInv12 + (kxzL*kyzL - kxyL*kzzL)*gInv31 + 
      gInv11*(-(kxyL*kxzL) + kxxL*kyzL + R1213) + gInv21*R1223 - gInv13*R1323 
      - gInv23*R2323 + gInv32*(-(kyyL*kzzL) + pow(kyzL,2));
    
    CCTK_REAL Rojo33 CCTK_ATTRIBUTE_UNUSED = (gInv12 + 
      gInv21)*(-(kxzL*kyzL) + kxyL*kzzL + R1323) + gInv11*(kxxL*kzzL + R1313 
      - pow(kxzL,2)) + gInv22*(kyyL*kzzL + R2323 - pow(kyzL,2));
    
    CCTK_REAL Psi4rL CCTK_ATTRIBUTE_UNUSED = 2*(n1*(n2*R4p1212 + 
      n3*R4p1213) - n3*(n2*R4p1223 + n3*R4p1323))*(imbar1*imbar2 - 
      rmbar1*rmbar2) + 2*(-(imbar2*imbar3) + rmbar2*rmbar3)*(n1*(n2*R4p1223 - 
      n3*R4p1323) - n2*n3*R4p2323 + R4p1213*pow(n1,2)) + 2*(imbar1*imbar3 - 
      rmbar1*rmbar3)*(n1*n2*R4p1213 + n1*n3*R4p1313 + n2*n3*R4p1323 + 
      R4p1223*pow(n2,2)) - (2*n2*n3*R4p1213 + R4p1212*pow(n2,2) + 
      R4p1313*pow(n3,2))*(pow(imbar1,2) - pow(rmbar1,2)) - (-2*n1*n3*R4p1223 
      + R4p1212*pow(n1,2) + R4p2323*pow(n3,2))*(pow(imbar2,2) - 
      pow(rmbar2,2)) - pow(nn,2)*((imbar2*imbar3 - rmbar2*rmbar3)*Rojo23 + 
      imbar1*(imbar2*(Rojo12 + Rojo21) + imbar3*(Rojo13 + Rojo31)) - 
      rmbar1*(rmbar2*(Rojo12 + Rojo21) + rmbar3*(Rojo13 + Rojo31)) + 
      (imbar2*imbar3 - rmbar2*rmbar3)*Rojo32 + Rojo11*(pow(imbar1,2) - 
      pow(rmbar1,2)) + Rojo22*(pow(imbar2,2) - pow(rmbar2,2)) + 
      Rojo33*(pow(imbar3,2) - pow(rmbar3,2))) - (2*n1*n2*R4p1323 + 
      R4p1313*pow(n1,2) + R4p2323*pow(n2,2))*(pow(imbar3,2) - pow(rmbar3,2)) 
      + 2*nn*((-(imbar1*imbar2) + rmbar1*rmbar2)*(n1*Ro112 + n2*Ro122 + 
      n3*Ro132) + (-(imbar1*imbar3) + rmbar1*rmbar3)*(n1*Ro113 + n2*Ro123 + 
      n3*Ro133) + (-(imbar1*imbar2) + rmbar1*rmbar2)*(n1*Ro211 + n2*Ro221 + 
      n3*Ro231) + (-(imbar2*imbar3) + rmbar2*rmbar3)*(n1*Ro213 + n2*Ro223 + 
      n3*Ro233) + (-(imbar1*imbar3) + rmbar1*rmbar3)*(n1*Ro311 + n2*Ro321 + 
      n3*Ro331) + (-(imbar2*imbar3) + rmbar2*rmbar3)*(n1*Ro312 + n2*Ro322 + 
      n3*Ro332) + (n1*Ro111 + n2*Ro121 + n3*Ro131)*(-pow(imbar1,2) + 
      pow(rmbar1,2)) + (n1*Ro212 + n2*Ro222 + n3*Ro232)*(-pow(imbar2,2) + 
      pow(rmbar2,2)) + (n1*Ro313 + n2*Ro323 + n3*Ro333)*(-pow(imbar3,2) + 
      pow(rmbar3,2)));
    
    CCTK_REAL Psi4iL CCTK_ATTRIBUTE_UNUSED = 2*((n1*(n2*R4p1212 + 
      n3*R4p1213) - n3*(n2*R4p1223 + n3*R4p1323))*(im2*rm1 + im1*rm2) + 
      nn*((im2*rm1 + im1*rm2)*(n1*(-Ro112 - Ro211) + n2*(-Ro122 - Ro221) + 
      n3*(-Ro132 - Ro231)) - 2*(im1*rm1*(n1*Ro111 + n2*Ro121 + n3*Ro131) + 
      im2*rm2*(n1*Ro212 + n2*Ro222 + n3*Ro232)) + (im3*rm1 + 
      im1*rm3)*(n1*(-Ro113 - Ro311) + n2*(-Ro123 - Ro321) + n3*(-Ro133 - 
      Ro331)) + (im3*rm2 + im2*rm3)*(n1*(-Ro213 - Ro312) + n2*(-Ro223 - 
      Ro322) + n3*(-Ro233 - Ro332)) - 2*im3*rm3*(n1*Ro313 + n2*Ro323 + 
      n3*Ro333)) + (im3*rm1 + im1*rm3)*(n1*(n2*R4p1213 + n3*R4p1313) + 
      n2*n3*R4p1323 + R4p1223*pow(n2,2))) - 2*((im3*rm2 + 
      im2*rm3)*(n1*(n2*R4p1223 - n3*R4p1323) - n2*n3*R4p2323 + 
      R4p1213*pow(n1,2)) + im3*rm3*(2*n1*n2*R4p1323 + R4p1313*pow(n1,2) + 
      R4p2323*pow(n2,2)) + im1*rm1*(2*n2*n3*R4p1213 + R4p1212*pow(n2,2) + 
      R4p1313*pow(n3,2)) + im2*rm2*(-2*n1*n3*R4p1223 + R4p1212*pow(n1,2) + 
      R4p2323*pow(n3,2))) - (im1*(2*rm1*Rojo11 + rm2*(Rojo12 + Rojo21) + 
      rm3*(Rojo13 + Rojo31)) + im2*(rm1*(Rojo12 + Rojo21) + 2*rm2*Rojo22 + 
      rm3*(Rojo23 + Rojo32)) + im3*(rm1*(Rojo13 + Rojo31) + rm2*(Rojo23 + 
      Rojo32) + 2*rm3*Rojo33))*pow(nn,2);
    
    CCTK_REAL Psi3rL CCTK_ATTRIBUTE_UNUSED = rm1*(n1*nn*((-2*ltet1 + 
      n1)*Ro111 - ltet2*Ro121 - ltet3*Ro131 + (-ltet2 + n2)*Ro211 + (-ltet3 + 
      n3)*Ro311) + n2*((-(ltet2*n1) + ltet1*n2)*R4p1212 + (-(ltet3*n1) + 
      ltet1*n3)*R4p1213 + (-(ltet3*n2) + ltet2*n3)*R4p1223 + nn*((-2*ltet1 + 
      n1)*Ro112 - ltet2*Ro122 - ltet3*Ro132 + (-ltet2 + n2)*Ro212 + (-ltet3 + 
      n3)*Ro312)) + n3*((-(ltet2*n1) + ltet1*n2)*R4p1213 + (-(ltet3*n1) + 
      ltet1*n3)*R4p1313 + (-(ltet3*n2) + ltet2*n3)*R4p1323 + nn*((-2*ltet1 + 
      n1)*Ro113 - ltet2*Ro123 - ltet3*Ro133 + (-ltet2 + n2)*Ro213 + (-ltet3 + 
      n3)*Ro313))) + rm2*(n1*((ltet2*n1 - ltet1*n2)*R4p1212 + (ltet3*n1 - 
      ltet1*n3)*R4p1213 + (ltet3*n2 - ltet2*n3)*R4p1223 + nn*((-ltet1 + 
      n1)*Ro121 - ltet1*Ro211 + (-2*ltet2 + n2)*Ro221 - ltet3*Ro231 + (-ltet3 
      + n3)*Ro321)) + n2*nn*((-ltet1 + n1)*Ro122 - ltet1*Ro212 + (-2*ltet2 + 
      n2)*Ro222 - ltet3*Ro232 + (-ltet3 + n3)*Ro322) + n3*((-(ltet2*n1) + 
      ltet1*n2)*R4p1223 + (-(ltet3*n1) + ltet1*n3)*R4p1323 + (-(ltet3*n2) + 
      ltet2*n3)*R4p2323 + nn*((-ltet1 + n1)*Ro123 - ltet1*Ro213 + (-2*ltet2 + 
      n2)*Ro223 - ltet3*Ro233 + (-ltet3 + n3)*Ro323))) + rm3*(n1*((ltet2*n1 - 
      ltet1*n2)*R4p1213 + (ltet3*n1 - ltet1*n3)*R4p1313 + (ltet3*n2 - 
      ltet2*n3)*R4p1323 + nn*((-ltet1 + n1)*Ro131 + (-ltet2 + n2)*Ro231 - 
      ltet1*Ro311 - ltet2*Ro321 + (-2*ltet3 + n3)*Ro331)) + n2*((ltet2*n1 - 
      ltet1*n2)*R4p1223 + (ltet3*n1 - ltet1*n3)*R4p1323 + (ltet3*n2 - 
      ltet2*n3)*R4p2323 + nn*((-ltet1 + n1)*Ro132 + (-ltet2 + n2)*Ro232 - 
      ltet1*Ro312 - ltet2*Ro322 + (-2*ltet3 + n3)*Ro332)) + n3*nn*((-ltet1 + 
      n1)*Ro133 + (-ltet2 + n2)*Ro233 - ltet1*Ro313 - ltet2*Ro323 + (-2*ltet3 
      + n3)*Ro333)) + ((ltet1 - n1)*(rm1*Rojo11 + rm2*Rojo12 + rm3*Rojo13) + 
      (ltet2 - n2)*(rm1*Rojo21 + rm2*Rojo22 + rm3*Rojo23) + (ltet3 - 
      n3)*(rm1*Rojo31 + rm2*Rojo32 + rm3*Rojo33))*pow(nn,2);
    
    CCTK_REAL Psi3iL CCTK_ATTRIBUTE_UNUSED = (im3*n2 - 
      im2*n3)*((-(ltet2*n1) + ltet1*n2)*R4p1223 + (-(ltet3*n1) + 
      ltet1*n3)*R4p1323 + (-(ltet3*n2) + ltet2*n3)*R4p2323) + 
      n1*nn*(im1*ltet1*Ro111 - im2*(-ltet1 + n1)*Ro121) + 
      im1*ltet2*n2*nn*Ro122 + im2*(ltet1 - n1)*n2*nn*Ro122 + 
      im1*ltet2*n3*nn*Ro123 + im2*(ltet1 - n1)*n3*nn*Ro123 + 
      im1*ltet3*n1*nn*Ro131 + n1*(-(im2*((ltet2*n1 - ltet1*n2)*R4p1212 + 
      (ltet3*n1 - ltet1*n3)*R4p1213 + (ltet3*n2 - ltet2*n3)*R4p1223)) - 
      im1*(-ltet1 + n1)*nn*Ro111 + im3*((-(ltet2*n1) + ltet1*n2)*R4p1213 + 
      (-(ltet3*n1) + ltet1*n3)*R4p1313 + (-(ltet3*n2) + ltet2*n3)*R4p1323 - 
      (-ltet1 + n1)*nn*Ro131)) + im1*ltet3*n2*nn*Ro132 + im3*(ltet1 - 
      n1)*n2*nn*Ro132 + im1*ltet3*n3*nn*Ro133 + im3*(ltet1 - n1)*n3*nn*Ro133 
      + im2*ltet1*n1*nn*Ro211 + im1*n1*(ltet2 - n2)*nn*Ro211 + 
      im2*ltet1*n2*nn*Ro212 + im1*(-(n3*((-(ltet2*n1) + ltet1*n2)*R4p1213 + 
      (-(ltet3*n1) + ltet1*n3)*R4p1313 + (-(ltet3*n2) + ltet2*n3)*R4p1323)) + 
      n2*((ltet2*n1 - ltet1*n2)*R4p1212 + (ltet3*n1 - ltet1*n3)*R4p1213 + 
      (ltet3*n2 - ltet2*n3)*R4p1223 - (-ltet2 + n2)*nn*Ro212)) + 
      im2*ltet1*n3*nn*Ro213 + im1*(ltet2 - n2)*n3*nn*Ro213 + 
      im2*ltet2*n1*nn*Ro221 + im2*n1*(ltet2 - n2)*nn*Ro221 + 
      im2*ltet2*n2*nn*Ro222 + n2*nn*(im1*ltet1*Ro112 - im2*(-ltet2 + 
      n2)*Ro222) + im2*ltet2*n3*nn*Ro223 + im2*(ltet2 - n2)*n3*nn*Ro223 + 
      im2*ltet3*n1*nn*Ro231 + im3*n1*(ltet2 - n2)*nn*Ro231 + 
      im2*ltet3*n2*nn*Ro232 + n2*nn*(im1*(ltet1 - n1)*Ro112 - im3*(-ltet2 + 
      n2)*Ro232) + im2*ltet3*n3*nn*Ro233 + im3*(ltet2 - n2)*n3*nn*Ro233 + 
      im3*ltet1*n1*nn*Ro311 + im1*n1*(ltet3 - n3)*nn*Ro311 + 
      im3*ltet1*n2*nn*Ro312 + im1*n2*(ltet3 - n3)*nn*Ro312 + 
      im3*ltet1*n3*nn*Ro313 + im1*n3*nn*(ltet1*Ro113 - (-ltet3 + n3)*Ro313) + 
      im3*ltet2*n1*nn*Ro321 + im2*n1*(ltet3 - n3)*nn*Ro321 + 
      im3*ltet2*n2*nn*Ro322 + im2*n2*(ltet3 - n3)*nn*Ro322 + 
      im3*ltet2*n3*nn*Ro323 + n3*nn*(im1*(ltet1 - n1)*Ro113 - im2*(-ltet3 + 
      n3)*Ro323) + im3*ltet3*n1*nn*Ro331 + im3*n1*(ltet3 - n3)*nn*Ro331 + 
      im3*ltet3*n2*nn*Ro332 + im3*n2*(ltet3 - n3)*nn*Ro332 + 
      im3*ltet3*n3*nn*Ro333 + nn*(im1*ltet2*n1*Ro121 - im3*n3*(-ltet3 + 
      n3)*Ro333) + (-ltet1 + n1)*(im1*Rojo11 + im2*Rojo12 + 
      im3*Rojo13)*pow(nn,2) + (-ltet2 + n2)*(im1*Rojo21 + im2*Rojo22 + 
      im3*Rojo23)*pow(nn,2) + (-ltet3 + n3)*(im1*Rojo31 + im2*Rojo32 + 
      im3*Rojo33)*pow(nn,2);
    
    CCTK_REAL Psi2rL CCTK_ATTRIBUTE_UNUSED = (ltet2*(n1*R4p1212 - 
      n3*R4p1223) + ltet3*(n1*R4p1213 - n3*R4p1323))*(im1*im2 + rm1*rm2) + 
      (ltet1*(n2*R4p1212 + n3*R4p1213) - ltet3*(n2*R4p1223 + 
      n3*R4p1323))*(im1*im2 + rm1*rm2) + (ltet2*n1*R4p1213 + ltet2*n2*R4p1223 
      + ltet3*n1*R4p1313 + ltet3*n2*R4p1323)*(im1*im3 + rm1*rm3) + 
      (ltet1*n2*R4p1213 + ltet2*n2*R4p1223 + ltet1*n3*R4p1313 + 
      ltet2*n3*R4p1323)*(im1*im3 + rm1*rm3) + (n1*(-(ltet1*R4p1213) - 
      ltet2*R4p1223) + ltet1*(-(n1*R4p1213) - n2*R4p1223 + n3*R4p1323) + 
      ltet2*n3*R4p2323 + ltet3*(n1*R4p1323 + n2*R4p2323))*(im2*im3 + rm2*rm3) 
      - (ltet2*n2*R4p1212 + ltet3*n2*R4p1213 + ltet2*n3*R4p1213 + 
      ltet3*n3*R4p1313)*(pow(im1,2) + pow(rm1,2)) - (n1*(ltet1*R4p1212 - 
      ltet3*R4p1223) + n3*(-(ltet1*R4p1223) + ltet3*R4p2323))*(pow(im2,2) + 
      pow(rm2,2)) - (ltet1*n1*R4p1313 + ltet2*n1*R4p1323 + ltet1*n2*R4p1323 + 
      ltet2*n2*R4p2323)*(pow(im3,2) + pow(rm3,2)) - 
      pow(nn,2)*(im1*im2*(Rojo12 + Rojo21) + rm1*rm2*(Rojo12 + Rojo21) + 
      im2*im3*Rojo23 + rm2*rm3*Rojo23 + im1*im3*(Rojo13 + Rojo31) + 
      rm1*rm3*(Rojo13 + Rojo31) + im2*im3*Rojo32 + rm2*rm3*Rojo32 + 
      Rojo11*pow(im1,2) + Rojo22*pow(im2,2) + Rojo33*pow(im3,2) + 
      Rojo11*pow(rm1,2) + Rojo22*pow(rm2,2) + Rojo33*pow(rm3,2)) + 
      nn*(n3*rm1*rm2*Ro213 + rm1*(-(ltet2*rm3*Ro123) + rm2*(n2*Ro122 - 
      ltet2*Ro221)) + rm3*(n1*rm1*Ro131 - ltet2*rm2*Ro223) + n1*rm2*rm3*Ro231 
      + rm1*(n2*rm3*Ro132 - ltet3*rm2*Ro231) + n2*rm2*rm3*Ro232 + 
      n3*rm2*rm3*Ro233 + im2*im3*(n1*Ro231 - ltet3*Ro233) + rm3*(n3*rm1*Ro133 
      - ltet3*rm2*Ro233) + n1*rm1*rm3*Ro311 + n2*rm1*rm3*Ro312 + 
      n3*rm1*rm3*Ro313 + im2*im3*n1*Ro321 + n1*rm2*rm3*Ro321 + 
      rm1*(rm2*(n1*Ro121 - ltet2*Ro122 - ltet3*Ro132) - ltet2*rm3*Ro321) + 
      im2*im3*n2*Ro322 + n2*rm2*rm3*Ro322 + im2*im3*(n2*Ro232 - ltet2*Ro322) 
      + rm2*(n1*rm1*Ro211 - ltet2*rm3*Ro322) + n3*rm2*rm3*Ro323 + 
      rm1*(n3*rm2*Ro123 + ltet3*rm3*(-Ro133 - Ro331)) + im1*(im2*((-ltet2 + 
      n2)*Ro122 + n3*Ro123 - ltet1*(Ro112 + Ro211) + n1*(Ro121 + Ro211) + 
      n2*Ro212 + n3*Ro213 - ltet2*Ro221 + ltet3*(-Ro132 - Ro231)) + 
      im3*(n1*Ro131 + n2*Ro132 + (-ltet3 + n3)*Ro133 + n1*Ro311 - 
      ltet1*(Ro113 + Ro311) + n2*Ro312 + n3*Ro313 - ltet2*(Ro123 + Ro321) - 
      ltet3*Ro331)) + im2*im3*(n3*Ro233 - ltet3*Ro332) + rm2*(n2*rm1*Ro212 - 
      ltet3*rm3*Ro332) + ((-ltet1 + n1)*Ro111 + n2*Ro112 + n3*Ro113 - 
      ltet2*Ro121 - ltet3*Ro131)*pow(im1,2) + (n1*Ro221 - 
      ltet2*Ro222)*pow(im2,2) + (n2*Ro222 - ltet3*Ro232)*pow(im2,2) + 
      Ro223*(-(im2*im3*ltet2) + n3*pow(im2,2)) + n2*Ro332*pow(im3,2) + 
      n3*Ro333*pow(im3,2) + (n1*Ro331 - ltet3*Ro333)*pow(im3,2) + 
      n3*Ro113*pow(rm1,2) + (n1*Ro111 - ltet2*Ro121)*pow(rm1,2) + (n2*Ro112 - 
      ltet3*Ro131)*pow(rm1,2) + (n1*Ro221 - ltet2*Ro222)*pow(rm2,2) + 
      n3*Ro223*pow(rm2,2) + (n2*Ro222 - ltet3*Ro232)*pow(rm2,2) + 
      Ro323*(im2*im3*n3 + ltet2*(-pow(im3,2) - pow(rm3,2))) + 
      n2*Ro332*pow(rm3,2) + n3*Ro333*pow(rm3,2) + (n1*Ro331 - 
      ltet3*Ro333)*pow(rm3,2) - ltet1*(rm1*rm2*(Ro112 + Ro211) + 
      rm2*rm3*Ro213 + rm1*rm3*(Ro113 + Ro311) + rm2*rm3*Ro312 + 
      im2*im3*(Ro213 + Ro312) + Ro212*pow(im2,2) + Ro313*pow(im3,2) + 
      Ro111*pow(rm1,2) + Ro212*pow(rm2,2) + Ro313*pow(rm3,2)));
    
    CCTK_REAL Psi2iL CCTK_ATTRIBUTE_UNUSED = (n1*(-(ltet2*R4p1212) - 
      ltet3*R4p1213) + ltet1*(n2*R4p1212 + n3*R4p1213) + (-(ltet3*n2) + 
      ltet2*n3)*R4p1223)*(im2*rm1 - im1*rm2) + ((-(ltet2*n1) + 
      ltet1*n2)*R4p1213 + (-(ltet3*n1) + ltet1*n3)*R4p1313 + (-(ltet3*n2) + 
      ltet2*n3)*R4p1323)*(im3*rm1 - im1*rm3) + (n1*(-(ltet1*R4p1213) - 
      ltet2*R4p1223) + ltet1*(n1*R4p1213 + n2*R4p1223) + ltet1*n3*R4p1323 + 
      ltet2*n3*R4p2323 - ltet3*(n1*R4p1323 + n2*R4p2323))*(im3*rm2 - im2*rm3) 
      - nn*(im1*(rm2*((-ltet2 - n2)*Ro122 - n3*Ro123 - ltet3*Ro132) + 
      rm3*(-(ltet2*Ro123) - n2*Ro132 + (-ltet3 - n3)*Ro133) + 
      ltet1*rm2*(-Ro112 + Ro211) + n1*rm2*(-Ro121 + Ro211) + n2*rm2*Ro212 + 
      n3*rm2*Ro213 + ltet2*rm2*Ro221 + ltet3*rm2*Ro231 + ltet1*rm3*(-Ro113 + 
      Ro311) + n1*rm3*(-Ro131 + Ro311) + n2*rm3*Ro312 + n3*rm3*Ro313 + 
      ltet2*rm3*Ro321 + ltet3*rm3*Ro331) + im2*(rm1*(ltet1*(Ro112 - Ro211) + 
      n1*(Ro121 - Ro211) + n2*(Ro122 - Ro212) + n3*(Ro123 - Ro213) + 
      ltet2*(Ro122 - Ro221) + ltet3*(Ro132 - Ro231)) + rm3*(-(ltet2*Ro223) - 
      n2*Ro232 + (-ltet3 - n3)*Ro233) + ltet1*rm3*(-Ro213 + Ro312) + 
      n1*rm3*(-Ro231 + Ro321) + ltet2*rm3*Ro322 + n2*rm3*Ro322 + n3*rm3*Ro323 
      + ltet3*rm3*Ro332) + im3*(rm1*((ltet3 + n3)*Ro133 + ltet1*(Ro113 - 
      Ro311) + n1*(Ro131 - Ro311) + n2*(Ro132 - Ro312) - n3*Ro313 + 
      ltet2*(Ro123 - Ro321) - ltet3*Ro331) + rm2*((ltet3 + n3)*Ro233 + 
      ltet1*(Ro213 - Ro312) + n1*(Ro231 - Ro321) + ltet2*(Ro223 - Ro322) + 
      n2*(Ro232 - Ro322) - n3*Ro323 - ltet3*Ro332))) + (rm1*(im2*(Rojo12 - 
      Rojo21) + im3*(Rojo13 - Rojo31)) + im1*(rm2*(-Rojo12 + Rojo21) + 
      rm3*(-Rojo13 + Rojo31)) + im3*rm2*(Rojo23 - Rojo32) + im2*rm3*(-Rojo23 
      + Rojo32))*pow(nn,2);
    
    CCTK_REAL Psi1rL CCTK_ATTRIBUTE_UNUSED = ltet1*(ltet1*(n2*R4p1212*rm2 
      + n3*R4p1213*rm2 + n2*R4p1213*rm3 + n3*R4p1313*rm3) - 
      ltet3*(n1*R4p1213*rm2 + n2*R4p1223*rm2 + n1*R4p1313*rm3 + 
      n2*R4p1323*rm3) + ltet2*(-(n1*(R4p1212*rm2 + R4p1213*rm3)) + 
      n3*(R4p1223*rm2 + R4p1323*rm3)) + nn*(-(n2*rm1) + (ltet1 - 
      n1)*rm2)*Ro121 + nn*(-(n3*rm1) + (ltet1 - n1)*rm3)*Ro131 + nn*((ltet2 - 
      n2)*rm1 - n1*rm2)*Ro211 + (ltet2 - 2*n2)*nn*rm2*Ro221 + nn*((ltet1 - 
      2*n1)*rm1*Ro111 - (n3*rm2 + (-ltet2 + n2)*rm3)*Ro231) + nn*((ltet3 - 
      n3)*rm1 - n1*rm3)*Ro311 + nn*((ltet3 - n3)*rm2 - n2*rm3)*Ro321 + (ltet3 
      - 2*n3)*nn*rm3*Ro331) + ltet2*(ltet1*((-(n2*R4p1212) - n3*R4p1213)*rm1 
      + n2*R4p1223*rm3 + n3*R4p1323*rm3) + ltet3*((n1*R4p1213 + 
      n2*R4p1223)*rm1 + (-(n1*R4p1323) - n2*R4p2323)*rm3) + 
      ltet2*((n1*R4p1212 - n3*R4p1223)*rm1 + (-(n1*R4p1223) + 
      n3*R4p2323)*rm3) + nn*((ltet1 - 2*n1)*rm1*Ro112 - (n2*rm1 + (-ltet1 + 
      n1)*rm2)*Ro122 - (n3*rm1 + (-ltet1 + n1)*rm3)*Ro132 + ((ltet2 - n2)*rm1 
      - n1*rm2)*Ro212 + (ltet2 - 2*n2)*rm2*Ro222 + (-(n3*rm2) + (ltet2 - 
      n2)*rm3)*Ro232 + ((ltet3 - n3)*rm1 - n1*rm3)*Ro312 + ((ltet3 - n3)*rm2 
      - n2*rm3)*Ro322 + (ltet3 - 2*n3)*rm3*Ro332)) + 
      ltet3*(-(ltet1*(n2*R4p1213*rm1 + n3*R4p1313*rm1 + n2*R4p1223*rm2 + 
      n3*R4p1323*rm2)) + ltet3*(n1*R4p1313*rm1 + n2*R4p1323*rm1 + 
      n1*R4p1323*rm2 + n2*R4p2323*rm2) + ltet2*(n1*(R4p1213*rm1 + 
      R4p1223*rm2) - n3*(R4p1323*rm1 + R4p2323*rm2)) + nn*((ltet2 - n2)*rm1 - 
      n1*rm2)*Ro213 + (ltet2 - 2*n2)*nn*rm2*Ro223 + nn*((ltet1 - 
      2*n1)*rm1*Ro113 - (n2*rm1 + (-ltet1 + n1)*rm2)*Ro123 - (n3*rm1 + 
      (-ltet1 + n1)*rm3)*Ro133 - (n3*rm2 + (-ltet2 + n2)*rm3)*Ro233) + 
      nn*((ltet3 - n3)*rm1 - n1*rm3)*Ro313 + nn*((ltet3 - n3)*rm2 - 
      n2*rm3)*Ro323 + (ltet3 - 2*n3)*nn*rm3*Ro333) - ((ltet1 - 
      n1)*(rm1*Rojo11 + rm2*Rojo12 + rm3*Rojo13) + rm1*((ltet2 - n2)*Rojo21 + 
      (ltet3 - n3)*Rojo31) + rm2*((ltet2 - n2)*Rojo22 + (ltet3 - n3)*Rojo32) 
      + rm3*((ltet2 - n2)*Rojo23 + (ltet3 - n3)*Rojo33))*pow(nn,2);
    
    CCTK_REAL Psi1iL CCTK_ATTRIBUTE_UNUSED = ltet1*(im2*((-(ltet2*n1) + 
      ltet1*n2)*R4p1212 + (-(ltet3*n1) + ltet1*n3)*R4p1213 + (-(ltet3*n2) + 
      ltet2*n3)*R4p1223 + (ltet2 - 2*n2)*nn*Ro221) + nn*(im1*(ltet1 - 
      2*n1)*Ro111 + (im2*(ltet1 - n1) - im1*n2)*Ro121 + (im3*(ltet1 - n1) - 
      im1*n3)*Ro131 - (im2*n1 + im1*(-ltet2 + n2))*Ro211 + (im3*(ltet2 - n2) 
      - im2*n3)*Ro231 - (im3*n1 + im1*(-ltet3 + n3))*Ro311 - (im3*n2 + 
      im2*(-ltet3 + n3))*Ro321) + im3*((-(ltet2*n1) + ltet1*n2)*R4p1213 + 
      (-(ltet3*n1) + ltet1*n3)*R4p1313 + (-(ltet3*n2) + ltet2*n3)*R4p1323 + 
      (ltet3 - 2*n3)*nn*Ro331)) + ltet2*(im1*((ltet2*n1 - ltet1*n2)*R4p1212 + 
      (ltet3*n1 - ltet1*n3)*R4p1213 + (ltet3*n2 - ltet2*n3)*R4p1223 + (ltet1 
      - 2*n1)*nn*Ro112) + nn*((im2*(ltet1 - n1) - im1*n2)*Ro122 + (im3*(ltet1 
      - n1) - im1*n3)*Ro132 + (-(im2*n1) + im1*(ltet2 - n2))*Ro212 + 
      im2*(ltet2 - 2*n2)*Ro222 + (im3*(ltet2 - n2) - im2*n3)*Ro232 - (im3*n1 
      + im1*(-ltet3 + n3))*Ro312 - (im3*n2 + im2*(-ltet3 + n3))*Ro322) + 
      im3*((-(ltet2*n1) + ltet1*n2)*R4p1223 + (-(ltet3*n1) + 
      ltet1*n3)*R4p1323 + (-(ltet3*n2) + ltet2*n3)*R4p2323 + (ltet3 - 
      2*n3)*nn*Ro332)) + ltet3*(im1*((ltet2*n1 - ltet1*n2)*R4p1213 + 
      (ltet3*n1 - ltet1*n3)*R4p1313 + (ltet3*n2 - ltet2*n3)*R4p1323 + (ltet1 
      - 2*n1)*nn*Ro113) + im2*((ltet2*n1 - ltet1*n2)*R4p1223 + (ltet3*n1 - 
      ltet1*n3)*R4p1323 + (ltet3*n2 - ltet2*n3)*R4p2323 + (ltet2 - 
      2*n2)*nn*Ro223) + nn*((im2*(ltet1 - n1) - im1*n2)*Ro123 + (im3*(ltet1 - 
      n1) - im1*n3)*Ro133 - (im2*n1 + im1*(-ltet2 + n2))*Ro213 + (im3*(ltet2 
      - n2) - im2*n3)*Ro233 + (-(im3*n1) + im1*(ltet3 - n3))*Ro313 + 
      (-(im3*n2) + im2*(ltet3 - n3))*Ro323 + im3*(ltet3 - 2*n3)*Ro333)) - 
      (im1*((ltet1 - n1)*Rojo11 + (ltet2 - n2)*Rojo21 + (ltet3 - n3)*Rojo31) 
      + im2*((ltet1 - n1)*Rojo12 + (ltet2 - n2)*Rojo22 + (ltet3 - n3)*Rojo32) 
      + im3*((ltet1 - n1)*Rojo13 + (ltet2 - n2)*Rojo23 + (ltet3 - 
      n3)*Rojo33))*pow(nn,2);
    
    CCTK_REAL Psi0rL CCTK_ATTRIBUTE_UNUSED = 2*(ltet1*(ltet2*R4p1212 + 
      ltet3*R4p1213) - ltet3*(ltet2*R4p1223 + ltet3*R4p1323))*(im1*im2 - 
      rm1*rm2) + 2*(-(im2*im3) + rm2*rm3)*(ltet1*(ltet2*R4p1223 - 
      ltet3*R4p1323) - ltet2*ltet3*R4p2323 + R4p1213*pow(ltet1,2)) + 
      2*(im1*im3 - rm1*rm3)*(ltet1*ltet2*R4p1213 + ltet1*ltet3*R4p1313 + 
      ltet2*ltet3*R4p1323 + R4p1223*pow(ltet2,2)) - (2*ltet2*ltet3*R4p1213 + 
      R4p1212*pow(ltet2,2) + R4p1313*pow(ltet3,2))*(pow(im1,2) - pow(rm1,2)) 
      - (-2*ltet1*ltet3*R4p1223 + R4p1212*pow(ltet1,2) + 
      R4p2323*pow(ltet3,2))*(pow(im2,2) - pow(rm2,2)) - pow(nn,2)*((im2*im3 - 
      rm2*rm3)*Rojo23 + im1*(im2*(Rojo12 + Rojo21) + im3*(Rojo13 + Rojo31)) - 
      rm1*(rm2*(Rojo12 + Rojo21) + rm3*(Rojo13 + Rojo31)) + (im2*im3 - 
      rm2*rm3)*Rojo32 + Rojo11*(pow(im1,2) - pow(rm1,2)) + Rojo22*(pow(im2,2) 
      - pow(rm2,2)) + Rojo33*(pow(im3,2) - pow(rm3,2))) - 
      (2*ltet1*ltet2*R4p1323 + R4p1313*pow(ltet1,2) + 
      R4p2323*pow(ltet2,2))*(pow(im3,2) - pow(rm3,2)) + 2*nn*((-(im1*im2) + 
      rm1*rm2)*(ltet1*Ro112 + ltet2*Ro122 + ltet3*Ro132) + (-(im1*im3) + 
      rm1*rm3)*(ltet1*Ro113 + ltet2*Ro123 + ltet3*Ro133) + (-(im1*im2) + 
      rm1*rm2)*(ltet1*Ro211 + ltet2*Ro221 + ltet3*Ro231) + (-(im2*im3) + 
      rm2*rm3)*(ltet1*Ro213 + ltet2*Ro223 + ltet3*Ro233) + (-(im1*im3) + 
      rm1*rm3)*(ltet1*Ro311 + ltet2*Ro321 + ltet3*Ro331) + (-(im2*im3) + 
      rm2*rm3)*(ltet1*Ro312 + ltet2*Ro322 + ltet3*Ro332) + (ltet1*Ro111 + 
      ltet2*Ro121 + ltet3*Ro131)*(-pow(im1,2) + pow(rm1,2)) + (ltet1*Ro212 + 
      ltet2*Ro222 + ltet3*Ro232)*(-pow(im2,2) + pow(rm2,2)) + (ltet1*Ro313 + 
      ltet2*Ro323 + ltet3*Ro333)*(-pow(im3,2) + pow(rm3,2)));
    
    CCTK_REAL Psi0iL CCTK_ATTRIBUTE_UNUSED = -2*(im3*rm1 + 
      im1*rm3)*(ltet1*(ltet2*R4p1213 + ltet3*R4p1313) + ltet2*ltet3*R4p1323 + 
      R4p1223*pow(ltet2,2)) + 2*((-(ltet1*(ltet2*R4p1212 + ltet3*R4p1213)) + 
      ltet3*(ltet2*R4p1223 + ltet3*R4p1323))*(im2*rm1 + im1*rm2) + 
      nn*((im2*rm1 + im1*rm2)*(ltet1*(Ro112 + Ro211) + ltet2*(Ro122 + Ro221) 
      + ltet3*(Ro132 + Ro231)) + (im3*rm1 + im1*rm3)*(ltet1*(Ro113 + Ro311) + 
      ltet2*(Ro123 + Ro321) + ltet3*(Ro133 + Ro331)) + (im3*rm2 + 
      im2*rm3)*(ltet1*(Ro213 + Ro312) + ltet2*(Ro223 + Ro322) + ltet3*(Ro233 
      + Ro332)) + 2*(im1*rm1*(ltet1*Ro111 + ltet2*Ro121 + ltet3*Ro131) + 
      im2*rm2*(ltet1*Ro212 + ltet2*Ro222 + ltet3*Ro232) + 
      im3*rm3*(ltet1*Ro313 + ltet2*Ro323 + ltet3*Ro333))) + (im3*rm2 + 
      im2*rm3)*(ltet1*(ltet2*R4p1223 - ltet3*R4p1323) - ltet2*ltet3*R4p2323 + 
      R4p1213*pow(ltet1,2)) + im3*rm3*(2*ltet1*ltet2*R4p1323 + 
      R4p1313*pow(ltet1,2) + R4p2323*pow(ltet2,2)) + 
      im1*rm1*(2*ltet2*ltet3*R4p1213 + R4p1212*pow(ltet2,2) + 
      R4p1313*pow(ltet3,2)) + im2*rm2*(-2*ltet1*ltet3*R4p1223 + 
      R4p1212*pow(ltet1,2) + R4p2323*pow(ltet3,2))) + (im1*(2*rm1*Rojo11 + 
      rm2*(Rojo12 + Rojo21) + rm3*(Rojo13 + Rojo31)) + im2*(rm1*(Rojo12 + 
      Rojo21) + 2*rm2*Rojo22 + rm3*(Rojo23 + Rojo32)) + im3*(rm1*(Rojo13 + 
      Rojo31) + rm2*(Rojo23 + Rojo32) + 2*rm3*Rojo33))*pow(nn,2);
    /* Copy local copies back to grid functions */
    Psi0i[index] = Psi0iL;
    Psi0r[index] = Psi0rL;
    Psi1i[index] = Psi1iL;
    Psi1r[index] = Psi1rL;
    Psi2i[index] = Psi2iL;
    Psi2r[index] = Psi2rL;
    Psi3i[index] = Psi3iL;
    Psi3r[index] = Psi3rL;
    Psi4i[index] = Psi4iL;
    Psi4r[index] = Psi4rL;
  }
  CCTK_ENDLOOP3(WeylScal4_psis_calc_4th);
}
extern "C" void WeylScal4_psis_calc_4th(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_WeylScal4_psis_calc_4th
  DECLARE_CCTK_ARGUMENTS_CHECKED(WeylScal4_psis_calc_4th);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WeylScal4_psis_calc_4th_Body");
  }
  if (cctk_iteration % WeylScal4_psis_calc_4th_calc_every != WeylScal4_psis_calc_4th_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::curv",
    "admbase::metric",
    "grid::coordinates",
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
  AssertGroupStorage(cctkGH, "WeylScal4_psis_calc_4th", 13, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psis_calc_4th", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psis_calc_4th", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psis_calc_4th", 2, 2, 2);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psis_calc_4th", 2, 2, 2);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, WeylScal4_psis_calc_4th_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving WeylScal4_psis_calc_4th_Body");
  }
}

} // namespace WeylScal4
