// test_array2.cc -- test driver for array.hh classes
// $Header$

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "stdc.h"
#include "util.hh"
#include "array.hh"

void  print(      jtutil::array2d<double>& Aref);
void cprint(const jtutil::array2d<double>& Aref);

int main()
{
const int imin = 3, imax = 7, jmin = 10, jmax = 12;
jtutil::array2d<double> A(imin,imax, jmin,jmax);

	  {
	for (int i = imin ; i <= imax ; ++i)
	{
	for (int j = jmin ; j <= jmax ; ++j)
	{
	A(i,j) = i*j;
	}
	}
	  }

 print(A);
cprint(A);
}

void print(jtutil::array2d<double>& Aref)
{
printf("=== nonconst ===\n");
	  {
	for (int i = Aref.min_i() ; i <= Aref.max_i() ; ++i)
	{
	for (int j = Aref.min_j() ; j <= Aref.max_j() ; ++j)
	{
	printf("Aref(%d,%d) = %g at address %p\n",
	       i, j, Aref(i,j), &Aref(i,j));
	}
	}
	  }
}

void cprint(const jtutil::array2d<double>& Aref)
{
printf("=== const ===\n");
	  {
	for (int i = Aref.min_i() ; i <= Aref.max_i() ; ++i)
	{
	for (int j = Aref.min_j() ; j <= Aref.max_j() ; ++j)
	{
	printf("Aref(%d,%d) = %g at address %p\n",
	       i, j, Aref(i,j), &Aref(i,j));
	}
	}
	  }
}
