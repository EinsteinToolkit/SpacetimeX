// test_modulo -- quick-n-dirty test of modulo_reduce()
// $Header$

#include <stdio.h>
#include <assert.h>

#include "stdc.h"
#include "util.hh"
using jtutil::error_exit;
using jtutil::fuzzy;

// prototypes
void tryit(double x, double xmod, double xmin, double xmax, double correct);

//******************************************************************************

int main()
{
//        x    xmod    xmin   xmax  correct
tryit(   0.5,   1.0,   0.0,   1.0,   0.5);
tryit(   0.0,   1.0,   0.0,   1.0,   0.0);
tryit(   1.0,   1.0,   0.0,   1.0,   1.0);
tryit(   2.0,   1.0,   0.0,   1.0,   1.0);
tryit(   1.5,   1.0,   0.0,   1.0,   0.5);
tryit(  -0.6,   1.0,   0.0,   1.0,   0.4);
tryit(  -3.6,   1.0,   0.0,   1.0,   0.4);

tryit(-145.0, 360.0, 180.0, 270.0, 215.0);

printf("all ok!\n");
}

//******************************************************************************

void tryit(double x, double xmod, double xmin, double xmax, double correct)
{
printf("trying %g mod %g [%g,%g] ==> ",
       x, xmod, xmin, xmax);

double got = jtutil::modulo_reduce(x, xmod, xmin, xmax);

if (fuzzy<double>::EQ(got, correct))
   then printf("ok\n");
   else {
	error_exit(ERROR_EXIT,
		   "no: got %g, correct %g\n", got, correct);	/*NOTREACHED*/
	}
}
