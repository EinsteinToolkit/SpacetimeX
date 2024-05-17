// test_coords2.cc -- test driver #2 for coordinate systems/conversions
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

#include "coords.hh"

using jtutil::fuzzy;
using jtutil::error_exit;
using jtutil::radians_of_degrees;
using jtutil::degrees_of_radians;

using namespace AHFinderDirect;
using namespace local_coords;

// prototypes
namespace {
void test_r_mu_nu_phi(fp x, fp y, fp z);
void test_r_theta_phi(fp x, fp y, fp z);
	  }

//******************************************************************************

//
// This program is a test driver for the local_coords:: coordinate-conversion
// functions.  It tries a large number of coordinate conversions, and checks
// them all for mutual consistency.
//
int main(int argc, const char* argv[])
{
bool verbose_flag = (argc == 2) && STRING_EQUAL(argv[1], "--verbose");
const fp xyz_minmax = 1.0;
const fp xyz_delta = 0.1;

	// these loops go *down* so we start in (+,+,+) quadrant
	// ==> get some successful tests before hard stuff
	// ==> help test the tests themselves
	for (fp x = xyz_minmax ; fuzzy<fp>::GE(x,-xyz_minmax) ; x -= xyz_delta)
	{
	for (fp y = xyz_minmax ; fuzzy<fp>::GE(y,-xyz_minmax) ; y -= xyz_delta)
	{
	for (fp z = xyz_minmax ; fuzzy<fp>::GE(z,-xyz_minmax) ; z -= xyz_delta)
	{
	// avoid places where angular coords might not be defined
	if (    fuzzy<fp>::EQ(x,0.0)
	     || fuzzy<fp>::EQ(y,0.0)
	     || fuzzy<fp>::EQ(z,0.0)    )
	   then continue;

	if (verbose_flag)
	   then printf("testing x=%g y=%g z=%g\n", x,y,z);

	test_r_mu_nu_phi(x, y, z);
	test_r_theta_phi(x, y, z);
	}
	}
	}

printf("all ok!\n");
return 0;
}

//******************************************************************************

namespace {
void test_r_mu_nu_phi(fp x, fp y, fp z)
{
const fp r = r_of_xyz(x, y, z);
const fp mu  =  mu_of_yz(y, z);
const fp nu  =  nu_of_xz(x, z);
const fp phi = phi_of_xy(x, y);

const fp  mu2 = mu_of_nu_phi(nu,phi);
const fp  nu2 = nu_of_mu_phi(mu,phi);
const fp phi2 = phi_of_mu_nu(mu,nu );
assert( fuzzy_EQ_ang(mu , mu2 ) );
assert( fuzzy_EQ_ang(nu , nu2 ) );
assert( fuzzy_EQ_ang(phi, phi2) );

fp x2, y2, z2;
xyz_of_r_mu_nu(r,mu,nu, x2,y2,z2);
assert( fuzzy<fp>::EQ(x, x2) );
assert( fuzzy<fp>::EQ(y, y2) );
assert( fuzzy<fp>::EQ(z, z2) );

fp x3, y3, z3;
xyz_of_r_mu_phi(r,mu,phi, x3,y3,z3);
assert( fuzzy<fp>::EQ(x, x3) );
assert( fuzzy<fp>::EQ(y, y3) );
assert( fuzzy<fp>::EQ(z, z3) );

fp x4, y4, z4;
xyz_of_r_nu_phi(r,nu,phi, x4,y4,z4);
assert( fuzzy<fp>::EQ(x, x4) );
assert( fuzzy<fp>::EQ(y, y4) );
assert( fuzzy<fp>::EQ(z, z4) );

}
	  }

//******************************************************************************

namespace {
void test_r_theta_phi(fp x, fp y, fp z)
{
const fp r        =     r_of_xyz(x, y, z);
const fp ps_theta = theta_of_xyz(x, y, z);
const fp ps_phi   =   phi_of_xy (x, y);

fp x2, y2, z2;
xyz_of_r_theta_phi(r,ps_theta,ps_phi, x2,y2,z2);
assert( fuzzy<fp>::EQ(x, x2) );
assert( fuzzy<fp>::EQ(y, y2) );
assert( fuzzy<fp>::EQ(z, z2) );

fp mu3, nu3;
mu_nu_of_theta_phi(ps_theta,ps_phi, mu3,nu3);
fp ps_theta3, ps_phi3;
theta_phi_of_mu_nu(mu3,nu3, ps_theta3,ps_phi3);
assert( fuzzy_EQ_ang(ps_theta3, ps_theta) );
assert( fuzzy_EQ_ang(ps_phi3  , ps_phi  ) );

fp mu4, phi4;
mu_phi_of_theta_phi(ps_theta,ps_phi, mu4,phi4);
fp ps_theta4, ps_phi4;
theta_phi_of_mu_phi(mu4,phi4, ps_theta4,ps_phi4);
assert( fuzzy_EQ_ang(ps_theta4, ps_theta) );
assert( fuzzy_EQ_ang(ps_phi4  , ps_phi  ) );

fp nu5, phi5;
nu_phi_of_theta_phi(ps_theta,ps_phi, nu5,phi5);
fp ps_theta5, ps_phi5;
theta_phi_of_nu_phi(nu5,phi5, ps_theta5,ps_phi5);
assert( fuzzy_EQ_ang(ps_theta5, ps_theta) );
assert( fuzzy_EQ_ang(ps_phi5  , ps_phi  ) );
}
	  }
