// test_norm -- quick-n-dirty test of norm template class
// $Header$

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "stdc.h"
#include "util.hh"
using jtutil::norm;
using jtutil::fuzzy;

int main()
{
norm<double> foo;
assert(  foo.is_empty()   );
assert(! foo.is_nonempty());
foo.data(3.0);
assert(! foo.is_empty()   );
assert(  foo.is_nonempty());
foo.data(-3.0);
foo.data(3.0);
foo.data(-3.0);
foo.data(-3.0);
assert( fuzzy<double>::EQ(foo.two_norm(), sqrt(45.0)) );
assert( fuzzy<double>::EQ(foo.rms_norm(), 3.0) );
assert( fuzzy<double>::EQ(foo.infinity_norm(), 3.0) );
assert( fuzzy<double>::EQ(foo.max_abs_value(), 3.0) );
assert( fuzzy<double>::EQ(foo.min_abs_value(), 3.0) );
assert( fuzzy<double>::EQ(foo.min_value(), -3.0) );
assert( fuzzy<double>::EQ(foo.max_value(), +3.0) );

foo.reset();
assert(  foo.is_empty() );
assert(! foo.is_nonempty());

foo.data(-3.0);
assert(! foo.is_empty() );
assert(  foo.is_nonempty());
foo.data(1.0);
foo.data(-4.0);
foo.data(1.0);
foo.data(-5.0);
assert( fuzzy<double>::EQ(foo.two_norm(), sqrt(52.0)) );
assert( fuzzy<double>::EQ(foo.rms_norm(), sqrt(10.4)) );
assert( fuzzy<double>::EQ(foo.infinity_norm(), 5.0) );
assert( fuzzy<double>::EQ(foo.max_abs_value(), 5.0) );
assert( fuzzy<double>::EQ(foo.min_abs_value(), 1.0) );
assert( fuzzy<double>::EQ(foo.min_value(), -5.0) );
assert( fuzzy<double>::EQ(foo.max_value(), +1.0) );

printf("everything ok!\n");
return 0;
}
