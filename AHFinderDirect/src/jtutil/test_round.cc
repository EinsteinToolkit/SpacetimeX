// test_round.cc -- test driver for round<float> and round<double>
// $Header$

#include <stdio.h>
#include <assert.h>

#include "stdc.h"
#include "util.hh"
using jtutil::round;

//
// Usage:
//	test_round x
//
int main(int argc, const char* const argv[])
{
	for (int i = -10 ; i <= +10 ; ++i)
	{
	printf("testing %d +/- various things...\n", i);

	float f = float(i);
	double d = double(i);


	//
	// round<fp_t>::to_integer()
	//

	assert( round<double>::to_integer(d - 0.49) == i );
	assert( round<double>::to_integer(d - 0.01) == i );
	assert( round<double>::to_integer(d - 1e-4) == i );
	assert( round<double>::to_integer(d - 1e-10) == i );
	assert( round<double>::to_integer(d        ) == i );
	assert( round<double>::to_integer(d + 1e-10) == i );
	assert( round<double>::to_integer(d + 1e-4) == i );
	assert( round<double>::to_integer(d + 0.01) == i );
	assert( round<double>::to_integer(d + 0.49) == i );

	assert( round<float>::to_integer(f - 0.49) == i );
	assert( round<float>::to_integer(f - 0.01) == i );
	assert( round<float>::to_integer(f - 1e-4) == i );
	assert( round<float>::to_integer(f - 1e-10) == i );
	assert( round<float>::to_integer(f        ) == i );
	assert( round<float>::to_integer(f + 1e-10) == i );
	assert( round<float>::to_integer(f + 1e-4) == i );
	assert( round<float>::to_integer(f + 0.01) == i );
	assert( round<float>::to_integer(f + 0.49) == i );


	//
	// round<fp_t>::floor()
	//

	assert( round<double>::floor(d - 0.49) == i-1 );
	assert( round<double>::floor(d - 0.01) == i-1);
	assert( round<double>::floor(d - 1e-4) == i-1 );
	assert( round<double>::floor(d - 1e-10) == i-1 );
	assert( round<double>::floor(d        ) == i );
	assert( round<double>::floor(d + 1e-10) == i );
	assert( round<double>::floor(d + 1e-4) == i );
	assert( round<double>::floor(d + 0.01) == i );
	assert( round<double>::floor(d + 0.49) == i );

	assert( round<float>::floor(f - 0.49) == i-1 );
	assert( round<float>::floor(f - 0.01) == i-1 );
	assert( round<float>::floor(f - 1e-4) == i-1 );
	assert( round<float>::floor(f - 1e-10) == ((i == 0) ? i-1 : i) );
					// i != 0 ==> not enough precision
					//            to see as noninteger
	assert( round<float>::floor(f        ) == i );
	assert( round<float>::floor(f + 1e-10) == i );
	assert( round<float>::floor(f + 1e-4) == i );
	assert( round<float>::floor(f + 0.01) == i );
	assert( round<float>::floor(f + 0.49) == i );


	//
	// round<fp_t>::ceiling()
	//

	assert( round<double>::ceiling(d - 0.49) == i );
	assert( round<double>::ceiling(d - 0.01) == i);
	assert( round<double>::ceiling(d - 1e-4) == i );
	assert( round<double>::ceiling(d - 1e-10) == i );
	assert( round<double>::ceiling(d        ) == i );
	assert( round<double>::ceiling(d + 1e-10) == i+1 );
	assert( round<double>::ceiling(d + 1e-4) == i+1 );
	assert( round<double>::ceiling(d + 0.01) == i+1 );
	assert( round<double>::ceiling(d + 0.49) == i+1 );

	assert( round<float>::ceiling(f - 0.49) == i );
	assert( round<float>::ceiling(f - 0.01) == i);
	assert( round<float>::ceiling(f - 1e-4) == i );
	assert( round<float>::ceiling(f - 1e-10) == i );
	assert( round<float>::ceiling(f        ) == i );
	assert( round<float>::ceiling(f + 1e-10) == ((i == 0) ? i+1 : i) );
					// i != 0 ==> not enough precision
					//            to see as noninteger
	assert( round<float>::ceiling(f + 1e-4) == i+1 );
	assert( round<float>::ceiling(f + 0.01) == i+1 );
	assert( round<float>::ceiling(f + 0.49) == i+1 );
	}

printf("everything looks ok!\n");
return 0;
}
