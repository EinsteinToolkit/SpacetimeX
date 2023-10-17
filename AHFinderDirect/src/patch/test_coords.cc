// test_coords.cc -- test driver for coordinate systems/conversions
// $Header$

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#ifdef STANDALONE_TEST
  #include "fake_cctk.h"
#else
  #include "cctk.h"
#endif

#include "config.h"
#include "stdc.h"
#include "../jtutil/util.hh"
using jtutil::error_exit;
using jtutil::radians_of_degrees;
using jtutil::degrees_of_radians;

#include "coords.hh"

using namespace AHFinderDirect;
using namespace local_coords;

//******************************************************************************

//
// This program is a test driver for the local_coords:: coordinate-conversion
// functions.  See the help message below for details.
//
int main(int argc, const char* argv[])
{
const char* help_msg =
"Usage:\n"
"n.b. all angles are input/output in degrees!\n"
"   test_coords  [  r=number ]   [ theta=number ]\n"
"                [ mu=number ]   [    nu=number ]   [ phi=number ]\n"
"                [  x=number ]   [     y=number ]   [   z=number ]\n"
;

//
// ***** command line arguments *****
//

if (    (argc == 1)
     || ((argc == 2) && STRING_EQUAL(argv[1], "--help"))    )
   then {
	printf("%s", help_msg);
	return 0;						/*NOTREACHED*/
	}

fp r, theta;
fp mu, nu, phi;
fp x, y, z;
fp dtheta;
fp dmu, dnu, dphi;
bool r_flag = false, theta_flag = false;
bool mu_flag = false, nu_flag = false, phi_flag = false;
bool x_flag = false, y_flag = false, z_flag = false;

	for (int ap = 1 ; ap < argc ; ++ap)
	{
	if	(sscanf(argv[ap], "r=" FP_SCANF_FORMAT, &r) == 1)
	   then r_flag = true;
	else if	(sscanf(argv[ap], "theta=" FP_SCANF_FORMAT, &dtheta) == 1)
	   then {
		theta = radians_of_degrees(dtheta);
		theta_flag = true;
		}

	else if	(sscanf(argv[ap], "mu=" FP_SCANF_FORMAT, &dmu) == 1)
	   then {
		mu = radians_of_degrees(dmu);
		mu_flag = true;
		}
	else if	(sscanf(argv[ap], "nu=" FP_SCANF_FORMAT, &dnu) == 1)
	   then {
		nu = radians_of_degrees(dnu);
		nu_flag = true;
		}
	else if	(sscanf(argv[ap], "phi=" FP_SCANF_FORMAT, &dphi) == 1)
	   then {
		phi = radians_of_degrees(dphi);
		phi_flag = true;
		}

	else if	(sscanf(argv[ap], "x=" FP_SCANF_FORMAT, &x) == 1)
	   then x_flag = true;
	else if	(sscanf(argv[ap], "y=" FP_SCANF_FORMAT, &y) == 1)
	   then y_flag = true;
	else if	(sscanf(argv[ap], "z=" FP_SCANF_FORMAT, &z) == 1)
	   then z_flag = true;

	else	error_exit(ERROR_EXIT, "%s", help_msg);		/*NOTREACHED*/
	}


//
// ***** input coordinates *****
//

printf("input:");
if (r_flag) then printf(" r=%g", r);
if (theta_flag) then printf(" theta=%g", dtheta);
if (mu_flag) then printf(" mu=%g", dmu);
if (nu_flag) then printf(" nu=%g", dnu);
if (phi_flag) then printf(" phi=%g", dphi);
if (x_flag) then printf(" x=%g", x);
if (y_flag) then printf(" y=%g", y);
if (z_flag) then printf(" z=%g", z);
printf("\n");


//
// ***** coordinate conversions *****
//

// ((mu,nu,phi)) --> the 3rd
if (mu_flag && nu_flag)
   then printf("phi_of_mu_nu() ==> phi=%g\n",
	       degrees_of_radians(phi_of_mu_nu(mu,nu)));
if (mu_flag && phi_flag)
   then printf("nu_of_mu_phi() ==> nu=%g\n",
	       degrees_of_radians(nu_of_mu_phi(mu,phi)));
if (nu_flag && phi_flag)
   then printf("mu_of_nu_phi() ==> mu=%g\n",
	       degrees_of_radians(mu_of_nu_phi(nu,phi)));

// (r,(mu,nu,phi)) --> (x,y,z)
if (r_flag && mu_flag && nu_flag)
   then {
	fp x2, y2, z2;
	xyz_of_r_mu_nu(r,mu,nu, x2,y2,z2);
	printf("xyz_of_r_mu_nu() ==> x=%g y=%g z=%g\n", x2,y2,z2);
	}
if (r_flag && mu_flag && phi_flag)
   then {
	fp x2, y2, z2;
	xyz_of_r_mu_phi(r,mu,phi, x2,y2,z2);
	printf("xyz_of_r_mu_phi() ==> x=%g y=%g z=%g\n", x2,y2,z2);
	}
if (r_flag && nu_flag && phi_flag)
   then {
	fp x2, y2, z2;
	xyz_of_r_nu_phi(r,nu,phi, x2,y2,z2);
	printf("xyz_of_r_nu_phi() ==> x=%g y=%g z=%g\n", x2,y2,z2);
	}

// (x,y,z) --> (r,(mu,nu,phi))
if (x_flag && y_flag && z_flag)
   then {
	fp r2 = r_of_xyz(x,y,z);
	fp mu2 = mu_of_yz(y,z);
	fp nu2 = nu_of_xz(x,z);
	fp phi2 = phi_of_xy(x,y);
	printf("*_of_xyz() ==> r=%g mu=%g nu=%g phi=%g\n",
	       r2,
	       degrees_of_radians(mu2),
	       degrees_of_radians(nu2),
	       degrees_of_radians(phi2));
	}

// usual polar spherical (r,theta,phi) <--> (x,y,z)
if (r_flag && theta_flag && phi_flag)
   then {
	fp x2, y2, z2;
	xyz_of_r_theta_phi(r,theta,phi, x2,y2,z2);
	printf("xyz_of_r_theta_phi() ==> x=%g y=%g z=%g\n", x2,y2,z2);
	}
if (x_flag && y_flag && z_flag)
   then {
	fp r2 = r_of_xyz(x,y,z);
	fp theta2 = theta_of_xyz(x,y,z);
	fp phi2 = phi_of_xy(x,y);
	printf("*_of_xyz() ==> r=%g theta=%g phi=%g\n",
	       r2,
	       degrees_of_radians(theta2),
	       degrees_of_radians(phi2));
	}

// ((mu,nu,phi)) --> usual polar spherical (theta,phi)
// ... note phi is the same coordinate in both systems
if (mu_flag && nu_flag)
   then {
	fp ps_theta, ps_phi;
	theta_phi_of_mu_nu(mu,nu, ps_theta,ps_phi);
	printf("theta_phi_of_mu_nu() ==> ps_theta=%g ps_phi=%g\n",
	       degrees_of_radians(ps_theta),degrees_of_radians(ps_phi));
	}
if (mu_flag && phi_flag)
   then {
	fp ps_theta, ps_phi;
	theta_phi_of_mu_phi(mu,phi, ps_theta,ps_phi);
	printf("theta_phi_of_mu_phi() ==> ps_theta=%g ps_phi=%g\n",
	       degrees_of_radians(ps_theta),degrees_of_radians(ps_phi));
	}
if (nu_flag && phi_flag)
   then {
	fp ps_theta, ps_phi;
	theta_phi_of_nu_phi(nu,phi, ps_theta,ps_phi);
	printf("theta_phi_of_nu_phi() ==> ps_theta=%g ps_phi=%g\n",
	       degrees_of_radians(ps_theta),degrees_of_radians(ps_phi));
	}

// usual polar spherical (theta,phi) --> ((mu,nu,phi))
// ... note phi is the same coordinate in both systems
if (theta_flag && phi_flag)
   then {
	fp mu2, nu2;
	mu_nu_of_theta_phi(theta,phi, mu2,nu2);
	printf("mu_nu_of_theta_phi() ==> mu=%g nu=%g\n",
	       degrees_of_radians(mu2),degrees_of_radians(nu2));
	}
if (theta_flag && phi_flag)
   then {
	fp mu2, mnp_phi;
	mu_phi_of_theta_phi(theta,phi, mu2,mnp_phi);
	printf("mu_phi_of_theta_phi() ==> mu=%g mnp_phi=%g\n",
	       degrees_of_radians(mu2),degrees_of_radians(mnp_phi));
	}
if (theta_flag && phi_flag)
   then {
	fp nu2, mnp_phi;
	nu_phi_of_theta_phi(theta,phi, nu2,mnp_phi);
	printf("nu_phi_of_theta_phi() ==> nu=%g mnp_phi=%g\n",
	       degrees_of_radians(nu2),degrees_of_radians(mnp_phi));
	}

return 0;
}
