// BH_diagnostics.cc -- compute/print BH diagnostics
// $Header$
//
// BH_diagnostics::BH_diagnostics - initialize a  struct BH_diagnostics
//
// BH_diagnostics::copy_to_buffer - copy diagnostics to buffer
// BH_diagnostics::copy_from_buffer - copy buffer to diagnostics
//
// BH_diagnostics::compute - compute BH diagnostics after an AH has been found
// BH_diagnostics::surface_integral - integrate gridfn over the 2-sphere
//
// print - print a line or two summarizing the diagnostics
// setup_output_file - create/open output file, write header describing fields
// output - write a (long) line of all the diagnostics
// store - copy the surface and the diagnostics into the SphericalSurface
//         variables
// save - copy the surface and the diagnostics into the Cactus variables
// load - set the surface and the diagnostics from the Cactus variables
//

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <vector>

#include "util_Table.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "config.h"
#include "stdc.h"
#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"
using jtutil::error_exit;

#include "../patch/coords.hh"
#include "../patch/grid.hh"
#include "../patch/fd_grid.hh"
#include "../patch/patch.hh"
#include "../patch/patch_edge.hh"
#include "../patch/patch_interp.hh"
#include "../patch/ghost_zone.hh"
#include "../patch/patch_system.hh"

#include "../elliptic/Jacobian.hh"

#include "../gr/gfns.hh"
#include "../gr/gr.hh"

#include "horizon_sequence.hh"
#include "BH_diagnostics.hh"
#include "driver.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** access to persistent data *****
//
extern struct state state;

//******************************************************************************

//
// This function initializes a  struct BH_diagnostics  to all zeros.
//
BH_diagnostics::BH_diagnostics()
	: origin_x(0.0), origin_y(0.0), origin_z(0.0),
	  centroid_x(0.0), centroid_y(0.0), centroid_z(0.0),
          quadrupole_xx(0.0), quadrupole_xy(0.0), quadrupole_xz(0.0),
          quadrupole_yy(0.0), quadrupole_yz(0.0), quadrupole_zz(0.0),
	  min_radius(0.0), max_radius(0.0), mean_radius(0.0),
	  min_x(0.0), max_x(0.0),
	  min_y(0.0), max_y(0.0),
	  min_z(0.0), max_z(0.0),
	  circumference_xy(0.0), circumference_xz(0.0), circumference_yz(0.0),
	  area(1.0),
          expansion(0.0),
          inner_expansion(0.0),
          product_expansion(0.0),
          mean_curvature(0.0),
          area_gradient(0.0),
          expansion_gradient(0.0),
          inner_expansion_gradient(0.0),
          product_expansion_gradient(0.0),
          mean_curvature_gradient(0.0),
          mean_curvature_minimum(0.0),
          mean_curvature_maximum(0.0),
          mean_curvature_integral(0.0)
{ }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function copies the diagnostics to a user-supplied buffer.
//
void BH_diagnostics::copy_to_buffer(CCTK_REAL buffer[N_buffer])
	const
{
buffer[posn__origin_x] = this->origin_x;
buffer[posn__origin_y] = this->origin_y;
buffer[posn__origin_z] = this->origin_z;

buffer[posn__centroid_x] = centroid_x;
buffer[posn__centroid_y] = centroid_y;
buffer[posn__centroid_z] = centroid_z;

buffer[posn__quadrupole_xx] = quadrupole_xx;
buffer[posn__quadrupole_xy] = quadrupole_xy;
buffer[posn__quadrupole_xz] = quadrupole_xz;
buffer[posn__quadrupole_yy] = quadrupole_yy;
buffer[posn__quadrupole_xz] = quadrupole_yz;
buffer[posn__quadrupole_zz] = quadrupole_zz;

buffer[posn__min_radius]  = min_radius;
buffer[posn__max_radius]  = max_radius;
buffer[posn__mean_radius] = mean_radius;

buffer[posn__min_x] = min_x;
buffer[posn__max_x] = max_x;
buffer[posn__min_y] = min_y;
buffer[posn__max_y] = max_y;
buffer[posn__min_z] = min_z;
buffer[posn__max_z] = max_z;

buffer[posn__circumference_xy]  = circumference_xy;
buffer[posn__circumference_xz]  = circumference_xz;
buffer[posn__circumference_yz]  = circumference_yz;
buffer[posn__area]              = area;
buffer[posn__expansion]         = expansion;
buffer[posn__inner_expansion]   = inner_expansion;
buffer[posn__product_expansion] = product_expansion;
buffer[posn__mean_curvature]    = mean_curvature;

buffer[posn__area_gradient]              = area_gradient;
buffer[posn__expansion_gradient]         = expansion_gradient;
buffer[posn__inner_expansion_gradient]   = inner_expansion_gradient;
buffer[posn__product_expansion_gradient] = product_expansion_gradient;
buffer[posn__mean_curvature_gradient]    = mean_curvature_gradient;

buffer[posn__mean_curvature_minimum]    = mean_curvature_minimum;
buffer[posn__mean_curvature_maximum]    = mean_curvature_maximum;
buffer[posn__mean_curvature_integral]   = mean_curvature_integral;
}

//******************************************************************************

//
// This function copies a user-supplied buffer to the diagnostics.
//
void BH_diagnostics::copy_from_buffer(const CCTK_REAL buffer[N_buffer])
{
this->origin_x = buffer[posn__origin_x];
this->origin_y = buffer[posn__origin_y];
this->origin_z = buffer[posn__origin_z];

centroid_x = buffer[posn__centroid_x];
centroid_y = buffer[posn__centroid_y];
centroid_z = buffer[posn__centroid_z];

quadrupole_xx = buffer[posn__quadrupole_xx];
quadrupole_xy = buffer[posn__quadrupole_xy];
quadrupole_xz = buffer[posn__quadrupole_xz];
quadrupole_yy = buffer[posn__quadrupole_yy];
quadrupole_yz = buffer[posn__quadrupole_yz];
quadrupole_zz = buffer[posn__quadrupole_zz];

 min_radius = buffer[posn__min_radius];
 max_radius = buffer[posn__max_radius];
mean_radius = buffer[posn__mean_radius];

min_x = buffer[posn__min_x];
max_x = buffer[posn__max_x];
min_y = buffer[posn__min_y];
max_y = buffer[posn__max_y];
min_z = buffer[posn__min_z];
max_z = buffer[posn__max_z];

 circumference_xy = buffer[posn__circumference_xy];
 circumference_xz = buffer[posn__circumference_xz];
 circumference_yz = buffer[posn__circumference_yz];
             area = buffer[posn__area];
        expansion = buffer[posn__expansion];
  inner_expansion = buffer[posn__inner_expansion];
product_expansion = buffer[posn__product_expansion];
   mean_curvature = buffer[posn__mean_curvature];

             area_gradient = buffer[posn__area_gradient];
        expansion_gradient = buffer[posn__expansion_gradient];
  inner_expansion_gradient = buffer[posn__inner_expansion_gradient];
product_expansion_gradient = buffer[posn__product_expansion_gradient];
   mean_curvature_gradient = buffer[posn__mean_curvature_gradient];

   mean_curvature_minimum  = buffer[posn__mean_curvature_minimum];
   mean_curvature_maximum  = buffer[posn__mean_curvature_maximum];
   mean_curvature_integral = buffer[posn__mean_curvature_integral];
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// Given that an apparent horizon has been found, this function computes
// various black hole diagnostics.
//
// Inputs (gridfns)
// h		# ghosted
// one		# nominal
// global_[xyz]	# nominal
//
// Bugs:
// The computation is rather inefficient -- we make many passes over the
// angular grid, instead of doing everything in one pass.
//
void BH_diagnostics::compute
	(const patch_system& ps,
         const fp the_area,
         const fp mean_expansion,
         const fp mean_inner_expansion,
         const fp mean_product_expansion,
         const fp mean_mean_curvature,
         const fp the_area_gradient,
         const fp mean_expansion_gradient,
         const fp mean_inner_expansion_gradient,
         const fp mean_product_expansion_gradient,
         const fp mean_mean_curvature_gradient,
	 const struct BH_diagnostics_info& BH_diagnostics_info)
{
//
// min/max radius of horizon
//
jtutil::norm<fp> h_norms;
ps.ghosted_gridfn_norms(gfns::gfn__h, h_norms);
min_radius = h_norms.min_abs_value();
max_radius = h_norms.max_abs_value();


//
// xyz bounding box of horizon
//

// compute bounding box of nominal grid
// ... this is only the stored part of the horizon if there are symmetries
jtutil::norm<fp> x_norms;
ps.gridfn_norms(gfns::gfn__global_x, x_norms);
min_x = x_norms.min_value();
max_x = x_norms.max_value();

jtutil::norm<fp> y_norms;
ps.gridfn_norms(gfns::gfn__global_y, y_norms);
min_y = y_norms.min_value();
max_y = y_norms.max_value();

jtutil::norm<fp> z_norms;
ps.gridfn_norms(gfns::gfn__global_z, z_norms);
min_z = z_norms.min_value();
max_z = z_norms.max_value();

// adjust the bounding box for the symmetries
#define REFLECT(origin_, max_)	(origin_ - (max_ - origin_))
switch	(ps.type())
	{
case patch_system::patch_system__full_sphere:
	break;
case patch_system::patch_system__plus_z_hemisphere:
	min_z = REFLECT(ps.origin_z(), max_z);
	break;
case patch_system::patch_system__plus_xy_quadrant_mirrored:
case patch_system::patch_system__plus_xy_quadrant_rotating:
	min_x = REFLECT(ps.origin_x(), max_x);
	min_y = REFLECT(ps.origin_y(), max_y);
	break;
case patch_system::patch_system__plus_xz_quadrant_mirrored:
case patch_system::patch_system__plus_xz_quadrant_rotating:
	min_x = REFLECT(ps.origin_x(), max_x);
	min_z = REFLECT(ps.origin_z(), max_z);
	break;
case patch_system::patch_system__plus_xyz_octant_mirrored:
case patch_system::patch_system__plus_xyz_octant_rotating:
	min_x = REFLECT(ps.origin_x(), max_x);
	min_y = REFLECT(ps.origin_y(), max_y);
	min_z = REFLECT(ps.origin_z(), max_z);
	break;
default:
	error_exit(PANIC_EXIT,
"***** BH_diagnostics::compute(): unknown patch system type()=(int)%d!\n"
"                                 (this should never happen!)\n",
		   int(ps.type()));				/*NOTREACHED*/
	}


//
// surface integrals
//
const fp integral_one = surface_integral(ps,
					 gfns::gfn__one, true, true, true,
					 BH_diagnostics_info.integral_method);
const fp integral_h = surface_integral(ps,
				       gfns::gfn__h, true, true, true,
				       BH_diagnostics_info.integral_method);
const fp integral_x = surface_integral(ps,
				       gfns::gfn__global_x, true, true, false,
				       BH_diagnostics_info.integral_method);
const fp integral_y = surface_integral(ps,
				       gfns::gfn__global_y, true, false, true,
				       BH_diagnostics_info.integral_method);
const fp integral_z = surface_integral(ps,
				       gfns::gfn__global_z, false, true, true,
				       BH_diagnostics_info.integral_method);
const fp integral_xx = surface_integral(ps,
                                        gfns::gfn__global_xx, true, true, true,
                                        BH_diagnostics_info.integral_method);
const fp integral_xy = surface_integral(ps,
                                        gfns::gfn__global_xy, true, false, false,
                                        BH_diagnostics_info.integral_method);
const fp integral_xz = surface_integral(ps,
                                        gfns::gfn__global_xz, false, true, false,
                                        BH_diagnostics_info.integral_method);
const fp integral_yy = surface_integral(ps,
                                        gfns::gfn__global_yy, true, true, true,
                                        BH_diagnostics_info.integral_method);
const fp integral_yz = surface_integral(ps,
                                        gfns::gfn__global_yz, false, false, true,
                                        BH_diagnostics_info.integral_method);
const fp integral_zz = surface_integral(ps,
                                        gfns::gfn__global_zz, true, true, true,
                                        BH_diagnostics_info.integral_method);


//
// originds
//
this->origin_x = ps.origin_x();
this->origin_y = ps.origin_y();
this->origin_z = ps.origin_z();


//
// centroids
//
centroid_x = integral_x / integral_one;
centroid_y = integral_y / integral_one;
centroid_z = integral_z / integral_one;


//
// quadrupoles
//
quadrupole_xx = integral_xx / integral_one;
quadrupole_xy = integral_xy / integral_one;
quadrupole_xz = integral_xz / integral_one;
quadrupole_yy = integral_yy / integral_one;
quadrupole_yz = integral_yz / integral_one;
quadrupole_zz = integral_zz / integral_one;


//
// area, mean radius, and mass
//
mean_radius = integral_h / integral_one;


//
// expansion
//
area              = the_area;
expansion         = mean_expansion;
inner_expansion   = mean_inner_expansion;
product_expansion = mean_product_expansion;
mean_curvature    = mean_mean_curvature;

area_gradient              = the_area_gradient;
expansion_gradient         = mean_expansion_gradient;
inner_expansion_gradient   = mean_inner_expansion_gradient;
product_expansion_gradient = mean_product_expansion_gradient;
mean_curvature_gradient    = mean_mean_curvature_gradient;

//
// minimum, maximum and the integral of the mean curvature
//
jtutil::norm<fp> mean_curvature_norms;
ps.gridfn_norms(gfns::gfn__mean_curvature, mean_curvature_norms);
mean_curvature_minimum = mean_curvature_norms.min_value();
mean_curvature_maximum = mean_curvature_norms.max_value();
mean_curvature_integral = surface_integral(ps,
                                           gfns::gfn__mean_curvature,
                                           true, true, true,
                                           BH_diagnostics_info.integral_method);


//
// circumferences
//
circumference_xy
  = ps.circumference("xy", gfns::gfn__h,
		     gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
					 gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
							     gfns::gfn__g_dd_33,
		     BH_diagnostics_info.integral_method);
circumference_xz
  = ps.circumference("xz", gfns::gfn__h,
		     gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
					 gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
							     gfns::gfn__g_dd_33,
		     BH_diagnostics_info.integral_method);
circumference_yz
  = ps.circumference("yz", gfns::gfn__h,
		     gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
					 gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
							     gfns::gfn__g_dd_33,
		     BH_diagnostics_info.integral_method);
}

//******************************************************************************

//
// This function computes the surface integral of a gridfn over the
// horizon.
//
//static
  fp BH_diagnostics::surface_integral
	(const patch_system& ps,
	 int src_gfn, bool src_gfn_is_even_across_xy_plane,
		      bool src_gfn_is_even_across_xz_plane,
		      bool src_gfn_is_even_across_yz_plane,
	 enum patch::integration_method method)
{
return ps.integrate_gridfn
	   (src_gfn, src_gfn_is_even_across_xy_plane,
		     src_gfn_is_even_across_xz_plane,
		     src_gfn_is_even_across_yz_plane,
	    gfns::gfn__h,
	    gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
				gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
						    gfns::gfn__g_dd_33,
	    method);
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function prints a line or two summarizing the diagnostics,
// using CCTK_VInfo().
//
void BH_diagnostics::print(int N_horizons, int hn)
	const
{
const fp m_irreducible = sqrt(area / (16*PI));
CCTK_VInfo(CCTK_THORNSTRING,
	   "AH %d/%d: r=%g at (%f,%f,%f)",
	   hn, N_horizons,
	   double(mean_radius),
	   double(centroid_x), double(centroid_y), double(centroid_z));
CCTK_VInfo(CCTK_THORNSTRING,
	   "AH %d/%d: area=%.10g m_irreducible=%.10g",
	   hn, N_horizons,
	   double(area), double(m_irreducible));
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function creates a BH-diagnostics output file, writes a suitable
// header comment identifying the fields to be written by  output() ,
// flushes the stream (to help in examining the output while Cactus is
// still running), and finally returns a stdio file pointer which can be
// used by  output()  to output data to the file.
//
FILE* BH_diagnostics::setup_output_file(const struct IO_info& IO_info,
					int N_horizons, int hn)
	const
{
char file_name_buffer[IO_info::file_name_buffer_size];

const char* directory = IO_info.BH_diagnostics_directory;
const int status = CCTK_CreateDirectory(IO_info.default_directory_permission,
					directory);
if (status < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   BH_diagnostics::setup_output_file():\n"
"        error %d trying to create output directory\n"
"        \"%s\"!"
		   ,
		   status,
		   directory);					/*NOTREACHED*/

snprintf(file_name_buffer, IO_info::file_name_buffer_size,
	 "%s/%s.ah%d.%s",
	 directory, IO_info.BH_diagnostics_base_file_name,
	 hn, IO_info.BH_diagnostics_file_name_extension);

const char *openMode;
if (IO_TruncateOutputFiles(state.cgi.GH) == 1)
   openMode = "w";
else
   openMode = "a";

FILE *fileptr = fopen(file_name_buffer, openMode);
if (fileptr == NULL)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   BH_diagnostics::setup_output_file():\n"
"        can't open BH-diagnostics output file\n"
"        \"%s\"!"
		   ,
		   file_name_buffer);				/*NOTREACHED*/

fprintf(fileptr, "# apparent horizon %d/%d\n", hn, N_horizons);
fprintf(fileptr, "#\n");
fprintf(fileptr, "# column  1 = cctk_iteration\n");
fprintf(fileptr, "# column  2 = cctk_time\n");
fprintf(fileptr, "# column  3 = centroid_x\n");
fprintf(fileptr, "# column  4 = centroid_y\n");
fprintf(fileptr, "# column  5 = centroid_z\n");
fprintf(fileptr, "# column  6 = min radius\n");
fprintf(fileptr, "# column  7 = max radius\n");
fprintf(fileptr, "# column  8 = mean radius\n");
fprintf(fileptr, "# column  9 = quadrupole_xx\n");
fprintf(fileptr, "# column 10 = quadrupole_xy\n");
fprintf(fileptr, "# column 11 = quadrupole_xz\n");
fprintf(fileptr, "# column 12 = quadrupole_yy\n");
fprintf(fileptr, "# column 13 = quadrupole_yz\n");
fprintf(fileptr, "# column 14 = quadrupole_zz\n");
fprintf(fileptr, "# column 15 = min x\n");
fprintf(fileptr, "# column 16 = max x\n");
fprintf(fileptr, "# column 17 = min y\n");
fprintf(fileptr, "# column 18 = max y\n");
fprintf(fileptr, "# column 19 = min z\n");
fprintf(fileptr, "# column 20 = max z\n");
fprintf(fileptr, "# column 21 = xy-plane circumference\n");
fprintf(fileptr, "# column 22 = xz-plane circumference\n");
fprintf(fileptr, "# column 23 = yz-plane circumference\n");
fprintf(fileptr, "# column 24 = ratio of xz/xy-plane circumferences\n");
fprintf(fileptr, "# column 25 = ratio of yz/xy-plane circumferences\n");
fprintf(fileptr, "# column 26 = area\n");
fprintf(fileptr, "# column 27 = m_irreducible\n");
fprintf(fileptr, "# column 28 = areal radius\n");
fprintf(fileptr, "# column 29 = expansion Theta_(l)\n");
fprintf(fileptr, "# column 30 = inner expansion Theta_(n)\n");
fprintf(fileptr, "# column 31 = product of the expansions\n");
fprintf(fileptr, "# column 32 = mean curvature\n");
fprintf(fileptr, "# column 33 = gradient of the areal radius\n");
fprintf(fileptr, "# column 34 = gradient of the expansion Theta_(l)\n");
fprintf(fileptr, "# column 35 = gradient of the inner expansion Theta_(n)\n");
fprintf(fileptr, "# column 36 = gradient of the product of the expansions\n");
fprintf(fileptr, "# column 37 = gradient of the mean curvature\n");
fprintf(fileptr, "# column 38 = minimum  of the mean curvature\n");
fprintf(fileptr, "# column 39 = maximum  of the mean curvature\n");
fprintf(fileptr, "# column 40 = integral of the mean curvature\n");
fflush(fileptr);

return fileptr;
}

//******************************************************************************

//
// This function outputs a BH-diagnostics line to a stdio stream, then
// flushes the stream (to help in examining the output while Cactus is
// still running).
//
// Arguments:
// BH_diagnostics = The BH diagnostics to be written
// fileptr = The stdio file pointer to append to
//
void BH_diagnostics::output(FILE*fileptr, const struct IO_info& IO_info)
	const
{
assert(fileptr != NULL);

fprintf(fileptr,
     //  cctk_iteration        min radius      mean radius
     //  ==  cctk_time         ======  max radius
     //  ==  ====  centroid_[xyz]      ======  ======
     //  ==  ====  ==========  ======  ======  ======
	"%d\t%.3f\t%f\t%f\t%f\t%#.10g\t%#.10g\t%#.10g\t",
	IO_info.time_iteration, double(IO_info.time),
	double(centroid_x), double(centroid_y), double(centroid_z),
	double(min_radius), double(max_radius), double(mean_radius));

fprintf(fileptr,
     //  quadrupole_{xx,xy,xz,yy,yz,zz}
     //  ================================================
	"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t",
	double(quadrupole_xx - centroid_x*centroid_x),
        double(quadrupole_xy - centroid_x*centroid_y),
	double(quadrupole_xz - centroid_x*centroid_z),
        double(quadrupole_yy - centroid_y*centroid_y),
	double(quadrupole_yz - centroid_y*centroid_z),
        double(quadrupole_zz - centroid_z*centroid_z));

fprintf(fileptr,
     //  {min,max}_x     {min,max}_y     {min,max}_z
     //  ==============  ==============  ==============
	"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t",
	double(min_x), double(max_x),
	double(min_y), double(max_y),
	double(min_z), double(max_z));

fprintf(fileptr,
     //  {xy,xz,yz}-plane         xz/xy  yz/xy
     //  circumferences          circumference
     //                          ratios
     //  ======================  ==============
	"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t",
	double(circumference_xy),
	double(circumference_xz),
	double(circumference_yz),
	double(circumference_xz / circumference_xy),
	double(circumference_yz / circumference_xy));

const fp m_irreducible = sqrt(area / (16*PI));;
const fp areal_radius = sqrt(area / (4*PI));
const fp areal_radius_gradient = sqrt(1 / (16*PI*area)) * area_gradient;
fprintf(fileptr,
     //  area    m_irre- areal   expan-  inner   prod.   mean    areal   expan-  inner   prod.   mean    mean    mean    mean
     //          ducible radius  sion    expan-  of the  curva-  radius  sion    expan-  of the  curva-	 curva-  curva-  curva-
     //                                  sion    expan-  ture    grad.   grad.   sion    exp.s   ture    ture    ture    ture
     //                                          sions                           grad.   grad.   grad.   min.    max.    integ.
     //  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======
     //                                                                                                
     //  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======
	"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\n",
	double(area), double(m_irreducible), double(areal_radius),
        double(expansion), double(inner_expansion), double(product_expansion), double(mean_curvature),
        double(areal_radius_gradient),
        double(expansion_gradient), double(inner_expansion_gradient), double(product_expansion_gradient), double(mean_curvature_gradient),
	double(mean_curvature_minimum), double(mean_curvature_maximum), double(mean_curvature_integral));

fflush(fileptr);
}

//******************************************************************************

//
// This function copies the BH-diagnostics into the SphericalSurface variables.
//
// Arguments:
// CCTK_ARGUMENTS
//
void BH_diagnostics::store(CCTK_ARGUMENTS,
                           const int horizon_number, const int surface_number)
	const
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  assert (surface_number >= 0 && surface_number < nsurfaces);
  
  assert(state.AH_data_array[horizon_number] != NULL);
  const struct AH_data& AH_data = *state.AH_data_array[horizon_number];
  
  // tell the world if we don't look for this horizon any more
  const int my_dont_find_after = dont_find_after_individual[horizon_number];
  const fp my_find_after_time = find_after_individual_time[horizon_number];
  const fp my_dont_find_after_time = dont_find_after_individual_time[horizon_number];
  if ((my_dont_find_after >= 0 and cctk_iteration > my_dont_find_after) or
      (my_dont_find_after_time > my_find_after_time and cctk_time > my_dont_find_after_time))
  {
    assert (! AH_data.search_flag);
    sf_active[surface_number] = 0;
  }
  
  // only try to copy AH info if we've found AHs at this time level
  if (! AH_data.search_flag) {
    sf_valid[surface_number] = 0;
    return;
  }
  
  // did we actually *find* this horizon?
  if (! AH_data.found_flag) {
    sf_valid[surface_number] = -1;
    return;
  }
  
  sf_active       [surface_number] = 1;
  sf_valid        [surface_number] = 1;
  sf_origin_x     [surface_number] = this->origin_x;
  sf_origin_y     [surface_number] = this->origin_y;
  sf_origin_z     [surface_number] = this->origin_z;
  sf_mean_radius  [surface_number] = mean_radius;
  sf_min_radius   [surface_number] = min_radius;
  sf_max_radius   [surface_number] = max_radius;
  sf_centroid_x   [surface_number] = centroid_x;
  sf_centroid_y   [surface_number] = centroid_y;
  sf_centroid_z   [surface_number] = centroid_z;
  sf_quadrupole_xx[surface_number] = quadrupole_xx - centroid_x*centroid_x;
  sf_quadrupole_xy[surface_number] = quadrupole_xy - centroid_x*centroid_y;
  sf_quadrupole_xz[surface_number] = quadrupole_xz - centroid_x*centroid_z;
  sf_quadrupole_yy[surface_number] = quadrupole_yy - centroid_y*centroid_y;
  sf_quadrupole_yz[surface_number] = quadrupole_yz - centroid_y*centroid_z;
  sf_quadrupole_zz[surface_number] = quadrupole_zz - centroid_z*centroid_z;
  sf_min_x        [surface_number] = min_x;
  sf_max_x        [surface_number] = max_x;
  sf_min_y        [surface_number] = min_y;
  sf_max_y        [surface_number] = max_y;
  sf_min_z        [surface_number] = min_z;
  sf_max_z        [surface_number] = max_z;
  sf_area         [surface_number] = area;
  
  std::vector<CCTK_REAL> xx, yy, zz, rr;
  xx.resize (sf_ntheta[surface_number] * sf_nphi[surface_number]);
  yy.resize (sf_ntheta[surface_number] * sf_nphi[surface_number]);
  zz.resize (sf_ntheta[surface_number] * sf_nphi[surface_number]);
  rr.resize (sf_ntheta[surface_number] * sf_nphi[surface_number]);
  size_t n;
  n = 0;
  for (int j=0; j<sf_nphi[surface_number]; ++j) {
    for (int i=0; i<sf_ntheta[surface_number]; ++i) {
      const CCTK_REAL theta
        = sf_origin_theta[surface_number] + i * sf_delta_theta[surface_number];
      const CCTK_REAL phi
        = sf_origin_phi[surface_number] + j * sf_delta_phi[surface_number];
      xx[n] = sf_origin_x[surface_number] + sin(theta) * cos(phi);
      yy[n] = sf_origin_y[surface_number] + sin(theta) * sin(phi);
      zz[n] = sf_origin_z[surface_number] + cos(theta);
      ++n;
    }
  }
  AHFinderDirect_radius_in_direction
    (horizon_number, sf_ntheta[surface_number] * sf_nphi[surface_number],
     &xx[0], &yy[0], &zz[0], &rr[0]);
  n = 0;
  for (int j=0; j<sf_nphi[surface_number]; ++j) {
    for (int i=0; i<sf_ntheta[surface_number]; ++i) {
      sf_radius[i + maxntheta * (j + maxnphi * surface_number)] = rr[n];
      ++n;
    }
  }
}

//******************************************************************************

//
// This function copies the BH-diagnostics into the Cactus variables.
//
// Arguments:
// CCTK_ARGUMENTS
//
void BH_diagnostics::save(CCTK_ARGUMENTS,
                          const int horizon_number)
	const
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  assert(state.AH_data_array[horizon_number] != NULL);
  const struct verbose_info& verbose_info = state.verbose_info;
  const struct AH_data& AH_data = *state.AH_data_array[horizon_number];
  patch_system& ps = *AH_data.ps_ptr;
  
  if (AH_data.search_flag &&
      AH_data.found_flag &&
      cctk_iteration > ah_centroid_iteration[horizon_number-1])
  {
    ah_centroid_x_p[horizon_number-1] = ah_centroid_x[horizon_number-1];
    ah_centroid_y_p[horizon_number-1] = ah_centroid_y[horizon_number-1];
    ah_centroid_z_p[horizon_number-1] = ah_centroid_z[horizon_number-1];
    ah_centroid_t_p[horizon_number-1] = ah_centroid_t[horizon_number-1];
    ah_centroid_valid_p[horizon_number-1] = ah_centroid_valid[horizon_number-1];
    ah_centroid_iteration_p[horizon_number-1] =
      ah_centroid_iteration[horizon_number-1];
  
    ah_centroid_x[horizon_number-1] = centroid_x;
    ah_centroid_y[horizon_number-1] = centroid_y;
    ah_centroid_z[horizon_number-1] = centroid_z;
    ah_centroid_t[horizon_number-1] = cctk_time;
    ah_centroid_valid[horizon_number-1] = AH_data.found_flag;
    ah_centroid_iteration[horizon_number-1] = cctk_iteration;
  }
  
  ah_origin_x[horizon_number-1] = ps.origin_x();
  ah_origin_y[horizon_number-1] = ps.origin_y();
  ah_origin_z[horizon_number-1] = ps.origin_z();
  
  for (int pn = 0; pn < ps.N_patches(); ++pn) {
    assert (pn < 6);
    patch& p = ps.ith_patch(pn);
    for (int i = 0, irho = p.min_irho(); irho <= p.max_irho(); ++i, ++irho) {
      assert (i <= max_N_zones_per_right_angle);
      for (int j = 0, isigma = p.min_isigma(); isigma <= p.max_isigma(); ++j, ++isigma) {
        assert (j <= max_N_zones_per_right_angle);
        ah_radius[i + (max_N_zones_per_right_angle+1) * (j + (max_N_zones_per_right_angle+1) * (pn + 6 * (horizon_number-1)))]
          = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
      }
    }
  }
  
  ah_initial_find_flag       [horizon_number-1] = AH_data.initial_find_flag;
  ah_really_initial_find_flag[horizon_number-1] = AH_data.really_initial_find_flag;
  ah_search_flag             [horizon_number-1] = AH_data.search_flag;
  ah_found_flag              [horizon_number-1] = AH_data.found_flag;
  if (verbose_info.print_algorithm_details) {
    printf ("AHF BH_diagnostics::save[%d] initial_find_flag=%d\n",        horizon_number, (int) AH_data.initial_find_flag);
    printf ("AHF BH_diagnostics::save[%d] really_initial_find_flag=%d\n", horizon_number, (int) AH_data.really_initial_find_flag);
    printf ("AHF BH_diagnostics::save[%d] search_flag=%d\n",              horizon_number, (int) AH_data.search_flag);
    printf ("AHF BH_diagnostics::save[%d] found_flag=%d\n",               horizon_number, (int) AH_data.found_flag);
  }
}

//******************************************************************************

//
// This function copies the BH-diagnostics from the Cactus variables.
//
// Arguments:
// CCTK_ARGUMENTS
//
void BH_diagnostics::load(CCTK_ARGUMENTS,
                          const int horizon_number)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  assert(state.AH_data_array[horizon_number] != NULL);
  const struct verbose_info& verbose_info = state.verbose_info;
  struct AH_data& AH_data = *state.AH_data_array[horizon_number];
  patch_system& ps = *AH_data.ps_ptr;
  
  // only use stored origins if horizon had not yet been found!
  if (ah_found_flag[horizon_number-1]) {
    ps.origin_x(ah_origin_x[horizon_number-1]);
    ps.origin_y(ah_origin_y[horizon_number-1]);
    ps.origin_z(ah_origin_z[horizon_number-1]);
  }
  
  for (int pn = 0; pn < ps.N_patches(); ++pn) {
    assert (pn < 6);
    patch& p = ps.ith_patch(pn);
    for (int i = 0, irho = p.min_irho(); irho <= p.max_irho(); ++i, ++irho) {
      assert (i <= max_N_zones_per_right_angle);
      for (int j = 0, isigma = p.min_isigma(); isigma <= p.max_isigma(); ++j, ++isigma) {
        assert (j <= max_N_zones_per_right_angle);
        p.ghosted_gridfn(gfns::gfn__h, irho,isigma) =
          ah_radius[i + (max_N_zones_per_right_angle+1) * (j + (max_N_zones_per_right_angle+1) * (pn + 6 * (horizon_number-1)))];
      }
    }
  }
  
  // recover the full ghosted-grid horizon shape
  // (we only save and load the nominal-grid shape)
  ps.synchronize();
  
  AH_data.initial_find_flag        = ah_initial_find_flag       [horizon_number-1];
  AH_data.really_initial_find_flag = ah_really_initial_find_flag[horizon_number-1];
  AH_data.search_flag              = ah_search_flag             [horizon_number-1];
  AH_data.found_flag               = ah_found_flag              [horizon_number-1];
  if (verbose_info.print_algorithm_details) {
    printf ("AHF BH_diagnostics::load[%d] initial_find_flag=%d\n",        horizon_number, (int) AH_data.initial_find_flag);
    printf ("AHF BH_diagnostics::load[%d] really_initial_find_flag=%d\n", horizon_number, (int) AH_data.really_initial_find_flag);
    printf ("AHF BH_diagnostics::load[%d] search_flag=%d\n",              horizon_number, (int) AH_data.search_flag);
    printf ("AHF BH_diagnostics::load[%d] found_flag=%d\n",               horizon_number, (int) AH_data.found_flag);
  }
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
