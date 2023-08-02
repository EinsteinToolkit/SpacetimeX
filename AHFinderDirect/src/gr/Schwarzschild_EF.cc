// Schwarzschild_EF.cc -- set up Schwarzschild/EF geometry variables
// $Header$
//
// <<<prototypes for functions local to this file>>>
// Schwarzschild_EF_geometry - top-level driver
/// Schwarzschild_EF_gij_xyz - (x,y,z) g_ij
/// Schwarzschild_EF_Kij_xyz - (x,y,z) g_ij
/// Schwarzschild_EF_gij_rthetaphi - (r,theta,phi) g_ij
/// Schwarzschild_EF_Kij_rthetaphi - (r,theta,phi) g_ij
/// tensor_xform_rthetaphi_to_xyz - tensor xform diag T_dd
///

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "util_Table.h"
#include "cctk.h"
#include "cctk_Arguments.h"

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

#include "gfns.hh"
#include "gr.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// ***** prototypes for functions local to this file *****
//

namespace {
void Schwarzschild_EF_gij_xyz(fp m, fp epsilon,
			      fp x, fp y, fp z,
			      fp& g_xx, fp& g_xy, fp& g_xz,
					fp& g_yy, fp& g_yz,
						  fp& g_zz);
void Schwarzschild_EF_Kij_xyz(fp m, fp epsilon,
			      fp x, fp y, fp z,
			      fp& K_xx, fp& K_xy, fp& K_xz,
					fp& K_yy, fp& K_yz,
						  fp& K_zz);

void Schwarzschild_EF_gij_rthetaphi(fp m,
				    fp r, fp theta, fp phi,
				    fp& g_rr, fp& g_theta_theta);
void Schwarzschild_EF_Kij_rthetaphi(fp m,
				    fp r, fp theta, fp phi,
				    fp& K_rr, fp& K_theta_theta);

void xform_from_rthetaphi_to_xyz(fp epsilon,
				 fp x, fp y, fp z,
				 fp T_rr, fp T_theta_theta,
				 fp& T_xx, fp& T_xy, fp& T_xz,
					   fp& T_yy, fp& T_yz,
						     fp& T_zz);
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function sets the geometry gridfns to their analytic values
// for a unit-mass Schwarzschild spacetime in Eddington-Finkelstein
// coordinates.
//
// The (r,theta,phi) $g_{ij}$ and $K_{ij}$ are as given in appendix A2
// of my (Jonathan Thornburg's) Ph.D thesis, available on my web page at
//	http://www.aei.mpg.de/~jthorn/phd/html/phd.html
// The (r,theta,phi) --> (x,y,z) tensor transformation is as worked out
// on pages 16-17 of my AHFinderDirect working notes.
// The $\partial_k g_{kl}$ is done by numerical finite differencing.
//
// Inputs (angular gridfns, on ghosted grid):
// ... defined on ghosted grid
// ... only values on nominal grid are actually used as input
//	h				# shape of trial surface
//
// Inputs (angular gridfns, all on the nominal grid):
//	global_[xyz]			# xyz positions of grid points
//
// Outputs (angular gridfns, all on the nominal grid):
//	g_dd_{11,12,13,22,23,33}			# $g_{ij}$
//	K_dd_{11,12,13,22,23,33}			# $K_{ij}$
//	partial_d_g_dd_{1,2,3}{11,12,13,22,23,33}	# $\partial_k g_{ij}$
//
void Schwarzschild_EF_geometry(patch_system& ps,
			       const struct geometry_info& gi,
			       bool msg_flag)
{
const fp mass      = gi.geometry__Schwarzschild_EF__mass;
const fp x_posn    = gi.geometry__Schwarzschild_EF__x_posn;
const fp y_posn    = gi.geometry__Schwarzschild_EF__y_posn;
const fp z_posn    = gi.geometry__Schwarzschild_EF__z_posn;
const fp epsilon   = gi.geometry__Schwarzschild_EF__epsilon;
const fp Delta_xyz = gi.geometry__Schwarzschild_EF__Delta_xyz;

if (msg_flag)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "      setting up exact Schwarzschild/EF geometry");
	CCTK_VInfo(CCTK_THORNSTRING,
		   "                 posn=(%g,%g,%g) mass=%g",
		   double(x_posn),double(y_posn),double(z_posn), double(mass));
	}

    for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
    {
    patch& p = ps.ith_patch(pn);

	for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
	{
	for (int isigma = p.min_isigma() ;
	     isigma <= p.max_isigma() ;
	     ++isigma)
	{
	const fp r = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
	const fp rho = p.rho_of_irho(irho);
	const fp sigma = p.sigma_of_isigma(isigma);
	fp local_x, local_y, local_z;
	p.xyz_of_r_rho_sigma(r,rho,sigma, local_x, local_y, local_z);
	const fp global_x = ps.origin_x() + local_x;
	const fp global_y = ps.origin_y() + local_y;
	const fp global_z = ps.origin_z() + local_z;

	const fp x_rel = global_x - x_posn;
	const fp y_rel = global_y - y_posn;
	const fp z_rel = global_z - z_posn;

	// compute g_ij and K_ij
	Schwarzschild_EF_gij_xyz(mass, epsilon,
				 x_rel, y_rel, z_rel,
				 p.gridfn(gfns::gfn__g_dd_11, irho,isigma),
				 p.gridfn(gfns::gfn__g_dd_12, irho,isigma),
				 p.gridfn(gfns::gfn__g_dd_13, irho,isigma),
				 p.gridfn(gfns::gfn__g_dd_22, irho,isigma),
				 p.gridfn(gfns::gfn__g_dd_23, irho,isigma),
				 p.gridfn(gfns::gfn__g_dd_33, irho,isigma));
	Schwarzschild_EF_Kij_xyz(mass, epsilon,
				 x_rel, y_rel, z_rel,
				 p.gridfn(gfns::gfn__K_dd_11, irho,isigma),
				 p.gridfn(gfns::gfn__K_dd_12, irho,isigma),
				 p.gridfn(gfns::gfn__K_dd_13, irho,isigma),
				 p.gridfn(gfns::gfn__K_dd_22, irho,isigma),
				 p.gridfn(gfns::gfn__K_dd_23, irho,isigma),
				 p.gridfn(gfns::gfn__K_dd_33, irho,isigma));

	fp g11p, g12p, g13p, g22p, g23p, g33p;
	fp g11m, g12m, g13m, g22m, g23m, g33m;

	// compute partial_x g_ij by finite differencing in xyz
	Schwarzschild_EF_gij_xyz(mass, epsilon,
				 x_rel+Delta_xyz, y_rel, z_rel,
				 g11p, g12p, g13p,
				       g22p, g23p,
					     g33p);
	Schwarzschild_EF_gij_xyz(mass, epsilon,
				 x_rel-Delta_xyz, y_rel, z_rel,
				 g11m, g12m, g13m,
				       g22m, g23m,
					     g33m);
	const fp fx = 1.0 / (2.0*Delta_xyz);
	p.gridfn(gfns::gfn__partial_d_g_dd_111,irho,isigma) = fx * (g11p-g11m);
	p.gridfn(gfns::gfn__partial_d_g_dd_112,irho,isigma) = fx * (g12p-g12m);
	p.gridfn(gfns::gfn__partial_d_g_dd_113,irho,isigma) = fx * (g13p-g13m);
	p.gridfn(gfns::gfn__partial_d_g_dd_122,irho,isigma) = fx * (g22p-g22m);
	p.gridfn(gfns::gfn__partial_d_g_dd_123,irho,isigma) = fx * (g23p-g23m);
	p.gridfn(gfns::gfn__partial_d_g_dd_133,irho,isigma) = fx * (g33p-g33m);

	// compute partial_y g_ij by finite differencing in xyz
	Schwarzschild_EF_gij_xyz(mass, epsilon,
				 x_rel, y_rel+Delta_xyz, z_rel,
				 g11p, g12p, g13p,
				       g22p, g23p,
					     g33p);
	Schwarzschild_EF_gij_xyz(mass, epsilon,
				 x_rel, y_rel-Delta_xyz, z_rel,
				 g11m, g12m, g13m,
				       g22m, g23m,
					     g33m);
	const fp fy = 1.0 / (2.0*Delta_xyz);
	p.gridfn(gfns::gfn__partial_d_g_dd_211,irho,isigma) = fy * (g11p-g11m);
	p.gridfn(gfns::gfn__partial_d_g_dd_212,irho,isigma) = fy * (g12p-g12m);
	p.gridfn(gfns::gfn__partial_d_g_dd_213,irho,isigma) = fy * (g13p-g13m);
	p.gridfn(gfns::gfn__partial_d_g_dd_222,irho,isigma) = fy * (g22p-g22m);
	p.gridfn(gfns::gfn__partial_d_g_dd_223,irho,isigma) = fy * (g23p-g23m);
	p.gridfn(gfns::gfn__partial_d_g_dd_233,irho,isigma) = fy * (g33p-g33m);

	// compute partial_x g_ij by finite differencing in xyz
	Schwarzschild_EF_gij_xyz(mass, epsilon,
				 x_rel, y_rel, z_rel+Delta_xyz,
				 g11p, g12p, g13p,
				       g22p, g23p,
					     g33p);
	Schwarzschild_EF_gij_xyz(mass, epsilon,
				 x_rel, y_rel, z_rel-Delta_xyz,
				 g11m, g12m, g13m,
				       g22m, g23m,
					     g33m);
	const fp fz = 1.0 / (2.0*Delta_xyz);
	p.gridfn(gfns::gfn__partial_d_g_dd_311,irho,isigma) = fz * (g11p-g11m);
	p.gridfn(gfns::gfn__partial_d_g_dd_312,irho,isigma) = fz * (g12p-g12m);
	p.gridfn(gfns::gfn__partial_d_g_dd_313,irho,isigma) = fz * (g13p-g13m);
	p.gridfn(gfns::gfn__partial_d_g_dd_322,irho,isigma) = fz * (g22p-g22m);
	p.gridfn(gfns::gfn__partial_d_g_dd_323,irho,isigma) = fz * (g23p-g23m);
	p.gridfn(gfns::gfn__partial_d_g_dd_333,irho,isigma) = fz * (g33p-g33m);
	}
	}
    }
}

//******************************************************************************

//
// This function computes the Schwarzschild/Eddington-Finkelstein
// 3-metric $g_{ij}$ in $(x,y,z)$ spatial coordinates.
//
namespace {
void Schwarzschild_EF_gij_xyz(fp m, fp epsilon,
			      fp x, fp y, fp z,
			      fp& g_xx, fp& g_xy, fp& g_xz,
					fp& g_yy, fp& g_yz,
						  fp& g_zz)
{
fp r, theta, phi;
local_coords::r_theta_phi_of_xyz(x,y,z, r,theta,phi);

fp g_rr, g_theta_theta;
Schwarzschild_EF_gij_rthetaphi(m,
			       r, theta, phi,
			       g_rr, g_theta_theta);
xform_from_rthetaphi_to_xyz(epsilon,
			    x, y, z,
			    g_rr, g_theta_theta,
			    g_xx, g_xy, g_xz,
				  g_yy, g_yz,
					g_zz);
}
	  }

//******************************************************************************

//
// This function computes the Schwarzschild/Eddington-Finkelstein
// 3-extrinsic curvature $K_{ij}$ in $(x,y,z)$ spatial coordinates.
//
namespace {
void Schwarzschild_EF_Kij_xyz(fp m, fp epsilon,
			      fp x, fp y, fp z,
			      fp& K_xx, fp& K_xy, fp& K_xz,
					fp& K_yy, fp& K_yz,
						  fp& K_zz)
{
fp r, theta, phi;
local_coords::r_theta_phi_of_xyz(x,y,z, r,theta,phi);

fp K_rr, K_theta_theta;
Schwarzschild_EF_Kij_rthetaphi(m,
			       r, theta, phi,
			       K_rr, K_theta_theta);
xform_from_rthetaphi_to_xyz(epsilon,
			    x, y, z,
			    K_rr, K_theta_theta,
			    K_xx, K_xy, K_xz,
				  K_yy, K_yz,
					K_zz);
}
	  }

//******************************************************************************

//
// This function computes the two independent components of the 3-metric
// $g_{ij}$ for Schwarzschild spacetime in Eddington-Finkelstein coordinates
// $(r,theta,phi)$, as per equation (A2.2) of my Ph.D thesis, rescaled in
// the obvious way for non-unit mass.
//
// By the spherical symmetry, $g_{\phi\phi} = g_{\theta\theta} \sin^2\theta$.
//
namespace {
void Schwarzschild_EF_gij_rthetaphi(fp m,
				    fp r, fp theta, fp phi,
				    fp& g_rr, fp& g_theta_theta)
{
g_rr          = 1.0 + 2.0*m/r;
g_theta_theta = r*r;
}
	  }

//******************************************************************************

//
// This function computes the two independent components of the 3-extrinsic
// curvature $K_{ij}$ for Schwarzschild spacetime in Eddington-Finkelstein
// coordinates $(r,theta,phi)$, as per equation (A2.2) of my Ph.D thesis,
// rescaled in the obvious way for non-unit mass.
//
// By the spherical symmetry, $K_{\phi\phi} = K_{\theta\theta} \sin^2\theta$.
//
//
namespace {
void Schwarzschild_EF_Kij_rthetaphi(fp m,
				    fp r, fp theta, fp phi,
				    fp& K_rr, fp& K_theta_theta)
{
fp temp = 1.0 / sqrt(1.0 + 2.0*m/r);
K_rr          = -(2.0/(r*r)) * (1.0 + m/r) * temp;
K_theta_theta = 2.0 * temp;
}
	  }

//******************************************************************************

//
// This function transforms a diagonal T_dd (rank 2 covariant) tensor
// in spherical symmetry,
//	T_dd = diag[ T_rr   T_theta_theta   T_theta_theta*sin(theta)**2 ]
// from (r,theta,phi) to (x,y,z) coordinates, as per pages 16-17 of my
// AHFInderDirect working notes.
//
// Arguments:
// epsilon = Tolerance parameter for deciding when to switch from
//	     generic expressions to z-axis limits: we switch iff
//	     sin^2 theta = (x^2+y^2)/r^2 is <= epsilon
//
namespace {
void xform_from_rthetaphi_to_xyz(fp epsilon,
				 fp x, fp y, fp z,
				 fp T_rr, fp T_theta_theta,
				 fp& T_xx, fp& T_xy, fp& T_xz,
					   fp& T_yy, fp& T_yz,
						     fp& T_zz)
{
const fp x2 = x*x;
const fp y2 = y*y;
const fp x2py2 = x2 + y2;
const fp z2 = z*z;
const fp r2 = x2 + y2 + z2;		const fp r4 = r2*r2;
const bool z_flag = x2py2/r2 <= epsilon;

T_xx = (x2/r2)*T_rr + (z_flag ? 1.0/r2 : (x2*z2/r2 + y2) / (r2*x2py2))
		      *T_theta_theta;
T_yy = (y2/r2)*T_rr + (z_flag ? 1.0/r2 : (y2*z2/r2 + x2) / (r2*x2py2))
		      *T_theta_theta;
T_zz = (z2/r2)*T_rr + (x2py2/r4)*T_theta_theta;
T_xy = (x*y/r2)*T_rr + (z_flag ? 0.0 : (z2/r2 - 1.0) * x*y / (r2*x2py2))
		       *T_theta_theta;
T_xz = (x*z/r2)*T_rr - (x*z/r4)*T_theta_theta;
T_yz = (y*z/r2)*T_rr - (y*z/r4)*T_theta_theta;
}
	  }

//******************************************************************************

	  }	// namespace AHFinderDirect
