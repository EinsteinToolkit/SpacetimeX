//
// This function sets up an ellipsoid in the gridfn h, using the
// formulas in "ellipsoid.maple" and the Maple-generated C code in
// "ellipsoid.c":
//
// ellipsoid has center (A,B,C), radius (a,b,c)
// angular coordinate system has center (U,V,W)
//
// direction cosines wrt angular coordinate center are (xcos,ycos,zcos)
// i.e. a point has coordinates (U+xcos*r, V+ycos*r, W+zcos*r)
//
// then the equation of the ellipsoid is
//	(U+xcos*r - A)^2     (V+ycos*r - B)^2     (W+zcos*r - C)^2
//	-----------------  +  ----------------  +  -----------------  =  1
//	        a^2                  b^2                   c^2
//
// to solve this, we introduce intermediate variables
//	AU = A - U
//	BV = B - V
//	CW = C - W
//
namespace {
void setup_ellipsoid_initial_guess
	(patch_system& ps,
	 fp global_center_x, fp global_center_y, fp global_center_z,
	 fp radius_x, fp radius_y, fp radius_z)
{
CCTK_VInfo(CCTK_THORNSTRING,
	   "   h = ellipsoid: global_center=(%g,%g,%g)",
	   global_center_x, global_center_y, global_center_z);
CCTK_VInfo(CCTK_THORNSTRING,
	   "                  radius=(%g,%g,%g)",
	   radius_x, radius_y, radius_z);

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);
		fp xcos, ycos, zcos;
		p.xyzcos_of_rho_sigma(rho,sigma, xcos,ycos,zcos);

		// set up variables used by Maple-generated code
		const fp AU = global_center_x - ps.origin_x();
		const fp BV = global_center_y - ps.origin_y();
		const fp CW = global_center_z - ps.origin_z();
		const fp a = radius_x;
		const fp b = radius_y;
		const fp c = radius_z;

		// compute the solutions r_plus and r_minus
		fp r_plus, r_minus;
		#include "ellipsoid.c"

		// exactly one of the solutions (call it r) should be positive
		fp r;
		if      ((r_plus > 0.0) && (r_minus < 0.0))
		   then r = r_plus;
		else if ((r_plus < 0.0) && (r_minus > 0.0))
		   then r = r_minus;
		else    CCTK_VWarn(-1, __LINE__, __FILE__, CCTK_THORNSTRING,
				   "\n"
"   expected exactly one r>0 solution, got 0 or 2!\n"
"   %s patch (irho,isigma)=(%d,%d) ==> (rho,sigma)=(%g,%g)\n"
"   direction cosines (xcos,ycos,zcos)=(%g,%g,%g)\n"
"   ==> r_plus=%g r_minus=%g\n"
				   ,
				   p.name(), irho, isigma,
				   double(rho), double(sigma),
				   double(xcos), double(ycos), double(zcos),
				   double(r_plus), double(r_minus));
		   						/*NOTREACHED*/

		// r = horizon radius at this grid point
		p.ghosted_gridfn(ghosted_gfns::gfn__h, irho,isigma) = r;
		}
		}
	}
}
	  }
