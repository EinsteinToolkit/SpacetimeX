// test_linear_map.cc -- test driver for linear_map class
// $Header$

#include <assert.h>
#include <stdio.h>

#include "stdc.h"
#include "util.hh"
#include "linear_map.hh"
using jtutil::fuzzy;
using jtutil::linear_map;

//******************************************************************************

//
// Usage:
//	test_linear_map
// If an argument is supplied (it doesn't matter what it is), then
// the array indices will be printed as they're tested.
//
int main(int argc, const char *const argv[])
{
const int min_int = 4;		const int max_int = 14;
const int middle_int = (min_int + max_int) / 2;
const int min_widen = 13;	const int max_widen = 9;
const double min_fp = 9.0;
const double delta_fp = 0.1;
const double max_fp = 10.0;
const double middle_fp = (min_fp + max_fp) / 2.0;

printf("constructing linear_map<double>...\n");
linear_map<double> lm_wide(min_int-min_widen, max_int+max_widen,
			   min_fp-min_widen*delta_fp,
			   delta_fp,
			   max_fp+max_widen*delta_fp);
linear_map<double> lm = *new linear_map<double>(lm_wide, min_int, max_int);

printf("checking bounds and range info...\n");
assert(lm.min_int() == min_int);
assert(lm.max_int() == max_int);
assert(lm.N_points() == jtutil::how_many_in_range(min_int,max_int));

assert(   lm.is_in_range(min_int)   );
assert(   lm.is_in_range(max_int)   );
assert( ! lm.is_in_range(min_int-1) );
assert( ! lm.is_in_range(max_int+1) );

assert( lm.clamp(min_int-100) == min_int    );
assert( lm.clamp(min_int-1  ) == min_int    );
assert( lm.clamp(min_int    ) == min_int    );
assert( lm.clamp(min_int+1  ) == min_int+1  );
assert( lm.clamp(middle_int ) == middle_int );
assert( lm.clamp(max_int-1  ) == max_int-1  );
assert( lm.clamp(max_int    ) == max_int    );
assert( lm.clamp(max_int+1  ) == max_int    );
assert( lm.clamp(max_int+100) == max_int    );

assert( fuzzy<double>::EQ(lm.origin(), lm.fp_of_int_unchecked(0)) );
assert( fuzzy<double>::EQ(lm.delta_fp(), delta_fp) );
assert( fuzzy<double>::EQ(lm.inverse_delta_fp(), 1.0/delta_fp) );
assert( fuzzy<double>::EQ(lm.min_fp(), min_fp) );
assert( fuzzy<double>::EQ(lm.max_fp(), max_fp) );

assert(   lm.is_in_range(min_fp)         );
assert(   lm.is_in_range(max_fp)         );
assert( ! lm.is_in_range(min_fp-1.0e-10) );
assert( ! lm.is_in_range(max_fp+1.0e-10) );

assert( fuzzy<double>::EQ(lm.clamp(min_fp - 100.0*delta_fp), min_fp) );
assert( fuzzy<double>::EQ(lm.clamp(min_fp -   0.4*delta_fp), min_fp) );
assert( fuzzy<double>::EQ(lm.clamp(min_fp), min_fp) );
assert( fuzzy<double>::EQ(lm.clamp(min_fp + 0.4*delta_fp),
				   min_fp + 0.4*delta_fp) );
assert( fuzzy<double>::EQ(lm.clamp(middle_fp), middle_fp) );
assert( fuzzy<double>::EQ(lm.clamp(max_fp - 0.4*delta_fp),
				   max_fp - 0.4*delta_fp) );
assert( fuzzy<double>::EQ(lm.clamp(max_fp), max_fp) );
assert( fuzzy<double>::EQ(lm.clamp(max_fp +   0.4*delta_fp), max_fp) );
assert( fuzzy<double>::EQ(lm.clamp(max_fp + 100.0*delta_fp), max_fp) );

printf("checking delta roundings...\n");
      {
    for (int i = 1 ; i <= 3 ; ++i)
    {
    printf("   testing i = %d...\n", i);
    assert( fuzzy<double>::EQ( lm.delta_fp_of_delta_int(i), i*delta_fp ) );
    assert( lm.delta_int_of_delta_fp(i*delta_fp) == i );
    double di = double(i);

    // test silent roundings
    assert( lm.delta_int_of_delta_fp((di-0.49)*delta_fp, lm.nia_round) == i );
    assert( lm.delta_int_of_delta_fp((di+0.49)*delta_fp, lm.nia_round) == i );

    // test floor/ceiling
    assert( lm.delta_int_of_delta_fp((di-0.01)*delta_fp, lm.nia_floor  ) == i-1 );
    assert( lm.delta_int_of_delta_fp((di-0.01)*delta_fp, lm.nia_ceiling) == i   );
    assert( lm.delta_int_of_delta_fp((di+0.01)*delta_fp, lm.nia_round  ) == i   );
    assert( lm.delta_int_of_delta_fp((di+0.01)*delta_fp, lm.nia_ceiling) == i+1 );
    }
      }

printf("this should produce a pair of warning messages...\n");
assert( lm.delta_int_of_delta_fp(2.99*delta_fp, lm.nia_warning) == 3 );
assert( lm.delta_int_of_delta_fp(3.01*delta_fp, lm.nia_warning) == 3 );

printf("checking mappings...\n");
  {
int i;		// must declare these outside loop
double x;	// because we have *two* loop variables!
	for (i = min_int, x = min_fp ;
	     i <= max_int ;
	     ++i, x += delta_fp)
	{
	printf("   testing i = %d of [%d,%d]...\n", i, min_int, max_int);
	// test normal int <--> fp conversions
	assert( fuzzy<double>::EQ( lm.fp_of_int(i), x) );
	assert( lm.int_of_fp(x) == i );

	// test fp --> int returning result as fp
	if (i > min_int)
		assert(
		   fuzzy<double>::EQ(
		      lm.fp_int_of_fp(x - 0.5*delta_fp), double(i) - 0.5
				    )
		      );
	if (i < max_int)
		assert(
		   fuzzy<double>::EQ(
		      lm.fp_int_of_fp(x + 0.5*delta_fp), double(i) + 0.5
				    )
		      );

	// test silent rounding
	if (i > min_int)
		assert( lm.int_of_fp(x - 0.49*delta_fp, lm.nia_round) == i );
	if (i < max_int)
		assert( lm.int_of_fp(x + 0.49*delta_fp, lm.nia_round) == i );

	// test floor/ceiling
	if (i > min_int)
		assert( lm.int_of_fp(x - 0.01*delta_fp, lm.nia_floor) == i-1 );
	if (i < max_int)
		assert( lm.int_of_fp(x + 0.01*delta_fp, lm.nia_ceiling) == i+1 );

	// test zero-origin conversions
	assert( lm.zero_origin_int(i) == i - min_int );
	assert( lm.map_int(lm.zero_origin_int(i)) == i );
	}
  }

printf("this should produce a pair of warning messages...\n");
assert( lm.int_of_fp(min_fp + 0.01*delta_fp,
		     linear_map<double>::nia_warning) == min_int );
assert( lm.int_of_fp(max_fp - 0.01*delta_fp,
		     linear_map<double>::nia_warning) == max_int );

printf("this should produce an error message...\n");
#ifdef TEST_RANGE
assert( lm.int_of_fp(min_fp - 0.01*delta_fp) == min_int );
#else
assert( lm.int_of_fp(min_fp + 0.01*delta_fp) == min_int );
#endif

return 0;
}
