// test_cpm_map.cc -- test driver for cpm_map class
// $Header$

#include <assert.h>
#include <stdio.h>

#include "stdc.h"
#include "util.hh"
#include "cpm_map.hh"
using jtutil::cpm_map;

//******************************************************************************

// prototypes

void test_shift_map(int min_i, int sample_i, int max_i,
		    int min_j, int sample_j, int max_j);
void test_mirror_map(int min_i, int sample_i, int max_i,
		     int min_j, int sample_i, int max_j,
		     double fixed_point);

void verify_shift_map(const cpm_map<double> sm,
		      const int min_i, const int sample_i, const int max_i,
		      const int min_j, const int sample_j, const int max_j);
void verify_mirror_map(const cpm_map<double> mm,
		       const int min_i, const int sample_i, const int max_i,
		       const int min_j, const int sample_j, const int max_j,
		       const double fixed_point);

//******************************************************************************
//******************************************************************************
//******************************************************************************

int main()
{
test_shift_map(4 , 10, 14,
	       42, 48, 52);

test_mirror_map(4 , 10, 14,
		16, 20, 26,
		15.0);
test_mirror_map(4 , 10, 14,
		17, 21, 27,
		15.5);

return 0;
}

//******************************************************************************

//
// this function tests a shift map
//
void test_shift_map(const int min_i, const int sample_i, const int max_i,
		    const int min_j, const int sample_j, const int max_j)
{
printf("testing shift_map [%d,%d,%d] --> [%d,%d,%d]...\n",
       min_i, sample_i, max_i,
       min_j, sample_j, max_j);
cpm_map<double> sm(min_i, max_i,
		   sample_j - sample_i);
verify_shift_map(sm,
		 min_i, sample_i, max_i,
		 min_j, sample_j, max_j);

  {
printf("testing generic_map (+) [%d,%d,%d] --> [%d,%d,%d]...\n",
       min_i, sample_i, max_i,
       min_j, sample_j, max_j);
cpm_map<double> gm(min_i, max_i,
		   sample_i, sample_j,
		   true);
verify_shift_map(gm,
		 min_i, sample_i, max_i,
		 min_j, sample_j, max_j);
  }

  {
const double delta = 0.42;
const double d_sample_i = double(sample_i) + delta;
const double d_sample_j = double(sample_j) + delta;
printf("testing generic_map (+) [%d,%g,%d] --> [%d,%g,%d]...\n",
       min_i, d_sample_i, max_i,
       min_j, d_sample_j, max_j);
cpm_map<double> gm(min_i, max_i,
		   d_sample_i, d_sample_j,
		   true);
verify_shift_map(gm,
		 min_i, sample_i, max_i,
		 min_j, sample_j, max_j);
  }
}

//******************************************************************************

//
// this function tests a mirror map
//
void test_mirror_map(const int min_i, const int sample_i, const int max_i,
		     const int min_j, const int sample_j, const int max_j,
		     const double fixed_point)
{
  {
printf("testing mirror_map [%d,%d,%d] --> [%d,%d,%d] with fixed point %g\n",
       min_i, sample_i, max_i,
       min_j, sample_j, max_j,
       fixed_point);
cpm_map<double> mm(min_i, max_i,
		   fixed_point);
verify_mirror_map(mm,
		  min_i, sample_i, max_i,
		  min_j, sample_j, max_j,
		  fixed_point);
  }

  {
printf("testing generic_map (-) [%d,%d,%d] --> [%d,%d,%d]...\n",
       min_i, sample_i, max_i,
       min_j, sample_j, max_j);
cpm_map<double> gm(min_i, max_i,
		   sample_i, sample_j,
		   false);
verify_mirror_map(gm,
		  min_i, sample_i, max_i,
		  min_j, sample_j, max_j,
		  fixed_point);
  }

  {
const double delta = 0.5;
const double d_sample_i = double(sample_i) + delta;
const double d_sample_j = double(sample_j) - delta;
printf("testing generic_map (-) [%d,%g,%d] --> [%d,%g,%d]...\n",
       min_i, d_sample_i, max_i,
       min_j, d_sample_j, max_j);
cpm_map<double> gm(min_i, max_i,
		   d_sample_i, d_sample_j,
		   false);
verify_mirror_map(gm,
		  min_i, sample_i, max_i,
		  min_j, sample_j, max_j,
		  fixed_point);
  }
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// this function verifies that a shift map is correct
//
void verify_shift_map(const cpm_map<double> sm,
		      const int min_i, const int sample_i, const int max_i,
		      const int min_j, const int sample_j, const int max_j)
{
assert(sm.min_i() == min_i);		assert(sm.min_j() == min_j);
assert(sm.max_i() == max_i);		assert(sm.max_j() == max_j);
assert(sm.N_points() == jtutil::how_many_in_range(min_i, max_i));

assert(! sm.in_domain(min_i-1));	assert(! sm.in_range(min_j-1));
assert(  sm.in_domain(min_i  ));	assert(  sm.in_range(min_j  ));
assert(  sm.in_domain(min_i+1));	assert(  sm.in_range(min_j+1));
assert(  sm.in_domain(max_i-1));	assert(  sm.in_range(max_j-1));
assert(  sm.in_domain(max_i  ));	assert(  sm.in_range(max_j  ));
assert(! sm.in_domain(max_i+1));	assert(! sm.in_range(max_j+1));

assert(sm.is_plus());
assert(! sm.is_minus());
assert(sm.sign() == +1);
assert(sm.fp_sign() == +1.0);

assert(sm.map(min_i) == min_j);
assert(sm.map(max_i) == max_j);
assert(sm.map(sample_i) == sample_j);

assert(sm.inv_map(min_j) == min_i);
assert(sm.inv_map(max_j) == max_i);
assert(sm.inv_map(sample_j) == sample_i);

printf("==> all ok!\n");
}

//******************************************************************************

//
// this function verifies that a mirror map is correct
//
void verify_mirror_map(const cpm_map<double> mm,
		       const int min_i, const int sample_i, const int max_i,
		       const int min_j, const int sample_j, const int max_j,
		       const double fixed_point)
{
assert(mm.min_i() == min_i);		assert(mm.min_j() == min_j);
assert(mm.max_i() == max_i);		assert(mm.max_j() == max_j);
assert(mm.N_points() == jtutil::how_many_in_range(min_i, max_i));

assert(! mm.in_domain(min_i-1));	assert(! mm.in_range(min_j-1));
assert(  mm.in_domain(min_i  ));	assert(  mm.in_range(min_j  ));
assert(  mm.in_domain(min_i+1));	assert(  mm.in_range(min_j+1));
assert(  mm.in_domain(max_i-1));	assert(  mm.in_range(max_j-1));
assert(  mm.in_domain(max_i  ));	assert(  mm.in_range(max_j  ));
assert(! mm.in_domain(max_i+1));	assert(! mm.in_range(max_j+1));

assert(! mm.is_plus());
assert(mm.is_minus());
assert(mm.sign() == -1);
assert(mm.fp_sign() == -1.0);

assert(mm.map(min_i) == max_j);
assert(mm.map(max_i) == min_j);
assert(mm.map(sample_i) == sample_j);

assert(mm.inv_map(min_j) == max_i);
assert(mm.inv_map(max_j) == min_i);
assert(mm.inv_map(sample_j) == sample_i);

printf("==> all ok!\n");
}
