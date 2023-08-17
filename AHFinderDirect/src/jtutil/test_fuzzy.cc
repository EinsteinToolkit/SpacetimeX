// test_fuzzy.cc -- test driver for fuzzy<float> and fuzzy<double>
// $Header$
//
// main
// check<fp_t>
// check_binary<fp_t>
// check_unary<fp_t>
// *** template instantiations

#include <stdio.h>
#include <assert.h>

#include "stdc.h"
#include "util.hh"
using jtutil::error_exit;
using jtutil::fuzzy;

// prototypes
template <typename fp_t>
  void check(const char* print_string,
	     fp_t x, fp_t y,
	     char xy_relation,
	     bool unary_flag,
	     bool x_is_fuzzy_integer,
	     int fuzzy_floor_of_x,
	     int fuzzy_ceiling_of_x);
template <typename fp_t>
  void check_binary(fp_t x, fp_t y,
		    char xy_relation);
template <typename fp_t>
  void check_unary(fp_t x,
		   bool x_is_fuzzy_integer,
		   int fuzzy_floor_of_x,
		   int fuzzy_ceiling_of_x);

//******************************************************************************

int main(int argc, const char* const argv[])
{
const bool print_flag = (argc > 1);
const char* double_str = print_flag ? "<double>" : NULL;
const char* float_str  = print_flag ? "<float>" : NULL;

	for (int ii = -10 ; ii <= +10 ; ++ii)
	{
	bool unary_flag = ((ii % 2) == 0);	// the "correct answers"
						// for the unary tests
						// assume f and d are integral

	float f = 0.5*float(ii);
	double d = 0.5*double(ii);
	int i = ii / 2;

	printf("testing %+-4g +/- various things...\n", d);


	//                        x        y  relation  do unary  is_integer  floor  ceiling
	//                                     (x,y)     tests?      (x)       (x)     (x)
	check<double>(double_str, d - 0.49 , d,   '<',   unary_flag,  false,     i-1,    i  );
	check<double>(double_str, d - 0.01 , d,   '<',   unary_flag,  false,     i-1,    i  );
	check<double>(double_str, d - 1e-6 , d,   '<',   unary_flag,  false,     i-1,    i  );
	check<double>(double_str, d - 1e-10, d,   '<',   unary_flag,  false,     i-1,    i  );
	// ... distinct floating point numbers, but within tolerance
	assert( d - 1e-14 != d );
	check<double>(double_str, d - 1e-14, d,   '=',   unary_flag,  true ,     i  ,    i  );
	check<double>(double_str, d        , d,   '=',   unary_flag,  true ,     i  ,    i  );
	// ... distinct floating point numbers, but within tolerance
	assert( d + 1e-14 != d );
	check<double>(double_str, d + 1e-14, d,   '=',   unary_flag,  true ,     i  ,    i  );
	check<double>(double_str, d + 1e-10, d,   '>',   unary_flag,  false,     i  ,    i+1);
	check<double>(double_str, d + 1e-6 , d,   '>',   unary_flag,  false,     i  ,    i+1);
	check<double>(double_str, d + 0.01 , d,   '>',   unary_flag,  false,     i  ,    i+1);
	check<double>(double_str, d + 0.49 , d,   '>',   unary_flag,  false,     i  ,    i+1);

	//                        x        y  relation  do unary  is_integer  floor  ceiling
	//                                     (x,y)     tests?      (x)       (x)     (x)
	check<float> (float_str , f - 0.49 , f,   '<',   unary_flag,  false,     i-1,    i  );
	check<float> (float_str , f - 0.01 , f,   '<',   unary_flag,  false,     i-1,    i  );
	check<float> (float_str , f - 1e-4 , f,   '<',   unary_flag,  false,     i-1,    i  );
	// ... distinct floating point numbers, but within tolerance
	assert( f - 1e-6 != f );
	check<float> (float_str , f - 1e-6 , f,   '=',   unary_flag,  true ,     i  ,    i  );
	// ... not enough precision to see as noninteger
	check<float> (float_str , f - 1e-10, f,   '=',   unary_flag,  true ,     i  ,    i  );
	check<float> (float_str , f        , f,   '=',   unary_flag,  true ,     i  ,    i  );
	// ... not enough precision to see as noninteger
	check<float> (float_str , f + 1e-10, f,   '=',   unary_flag,  true ,     i  ,    i  );
	// ... distinct floating point numbers, but within tolerance
	assert( f + 1e-6 != f );
	check<float> (float_str , f + 1e-6 , f,   '=',   unary_flag,  true ,     i  ,    i  );
	check<float> (float_str , f + 0.01 , f,   '>',   unary_flag,  false,     i  ,    i+1);
	check<float> (float_str , f + 1e-4 , f,   '>',   unary_flag,  false,     i  ,    i+1);
	check<float> (float_str , f + 0.49 , f,   '>',   unary_flag,  false,     i  ,    i+1);

	if (print_flag)
	   then printf("\n");
	}

printf("everything looks ok!\n");
return 0;
}

//******************************************************************************

//
// This function template is a driver to do all the tests for a single
// argument or pair of arguments.  It calls  check_binary<fp_t>()  to do
// the binary-function tests, and optionally calls  check_unary<fp_t>()
// to do the unary-function tests.
//
template <typename fp_t>
  void check(const char* print_string,
	     fp_t x, fp_t y,
	     char xy_relation,
	     bool unary_flag,
	     bool x_is_fuzzy_integer, int fuzzy_floor_of_x,
				      int fuzzy_ceiling_of_x)
{
if (print_string != NULL)
   then printf("   testing %s %s: x=%.15f y=%.15f\n",
	       print_string,
	       unary_flag ? "op(x,y) + op(x)" : "op(x,y) only   ",
	       double(x), double(y));

check_binary<fp_t>(x, y, xy_relation);

if (unary_flag)
   then check_unary(x,
		    x_is_fuzzy_integer, fuzzy_floor_of_x, fuzzy_ceiling_of_x);
}

//******************************************************************************

//
// This function template tests the binary relations
//	fuzzy<fp_t>::EQ()
//	fuzzy<fp_t>::NE()
//	fuzzy<fp_t>::LT()
//	fuzzy<fp_t>::GT()
//	fuzzy<fp_t>::LE()
//	fuzzy<fp_t>::GE()
// for a single argument.
//
// The  xy_relation  argument must be one of '<', '=', or '>', indicating
// the desired fuzzy relationship between the two arguments.  The function
// infers the correct values of the binary functions from this.
//
template <typename fp_t>
  void check_binary(fp_t x, fp_t y,
		    char xy_relation)
{
//
// basic properties which should always hold
//

// trichotomy (we have exactly one of <, =, >
assert( fuzzy<fp_t>::LT(x,y) + fuzzy<fp_t>::EQ(x,y) + fuzzy<fp_t>::GT(x,y)
	== 1 );

// consistency of EQ() with NE()
assert( fuzzy<fp_t>::NE(x,y) == ! fuzzy<fp_t>::EQ(x,y) );

// consistency of LE() with LT() and EQ(), GE() with GT() and EQ()
assert( fuzzy<fp_t>::LE(x,y)
	== (fuzzy<fp_t>::LT(x,y) || fuzzy<fp_t>::EQ(x,y)) );
assert( fuzzy<fp_t>::GE(x,y)
	== (fuzzy<fp_t>::GT(x,y) || fuzzy<fp_t>::EQ(x,y)) );


//
// now check desired ordering relation
// ... by virtue of the previous tests, we need only check LT(), EQ(), and GT()
// ... in fact, by virtue of trichotomy (already tested), we really
//     only need to assert() the true one for each case, but the slight
//     redundancy of checking all of them seems nicer in this context
//
switch	(xy_relation)
	{
case '<':
	assert( fuzzy<fp_t>::LT(x,y) == true  );
	assert( fuzzy<fp_t>::EQ(x,y) == false );
	assert( fuzzy<fp_t>::GT(x,y) == false );
	break;
case '=':
	assert( fuzzy<fp_t>::LT(x,y) == false );
	assert( fuzzy<fp_t>::EQ(x,y) == true  );
	assert( fuzzy<fp_t>::GT(x,y) == false );
	break;
case '>':
	assert( fuzzy<fp_t>::LT(x,y) == false );
	assert( fuzzy<fp_t>::EQ(x,y) == false );
	assert( fuzzy<fp_t>::GT(x,y) == true  );
	break;
default:
	error_exit(PANIC_EXIT,
"***** check_binary<fp_t>: bad xy_relation=(int)'%c'\n"
"                          (this should never happen!)\n"
,
		   int(xy_relation));				/*NOTREACHED*/
	}
}

//******************************************************************************

//
// This function template tests the unary functions
//	fuzzy<fp_t>::is_integer()
//	fuzzy<fp_t>::floor()
//	fuzzy<fp_t>::ceiling()
// for a single argument.
//
template <typename fp_t>
  void check_unary(fp_t x,
		   // remaining arguments are the correct results of...
		   bool x_is_fuzzy_integer,	// fuzzy<fp_t>::is_integer(x)
		   int fuzzy_floor_of_x,	// fuzzy<fp_t>::floor(x)
		   int fuzzy_ceiling_of_x)	// fuzzy<fp_t>::ceiling(x)
{
assert( fuzzy<fp_t>::is_integer(x) == x_is_fuzzy_integer );
assert( fuzzy<fp_t>::floor(x) == fuzzy_floor_of_x );
assert( fuzzy<fp_t>::ceiling(x) == fuzzy_ceiling_of_x );
}

//******************************************************************************

//
// template instantiations for <float> and <double>
//

template void check<float>(const char*,
			   float, float,
			   char, bool, bool, int, int);
template void check<double>(const char*,
			    double, double,
			    char, bool, bool, int, int);

template void check_binary<float>(float, float, char);
template void check_binary<double>(double, double, char);

template void check_unary<float>(float, bool, int, int);
template void check_unary<double>(double, bool, int, int);
