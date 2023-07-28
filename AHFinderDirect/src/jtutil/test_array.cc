// test_array.cc -- test driver for array.hh classes
// $Header$

#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "stdc.h"
#include "util.hh"
#include "array.hh"

using jtutil::fuzzy;
using jtutil::array1d;
using jtutil::array2d;
using jtutil::array3d;
#ifdef NOT_USED
using jtutil::array4d;
#endif

//
// prototypes for non-class functions defined in this file
//
double test_array1d_double(bool fancy_storage, bool print_flag);
double test_array2d_float(bool fancy_storage, bool print_flag);
double test_array3d_double(bool fancy_storage, bool print_flag);
#ifdef NOT_USED
double test_array4d_double(bool fancy_storage, bool print_flag);
#endif

//******************************************************************************

//
// Usage:
//	test_array  [ print ]
// If an argument is supplied (it doesn't matter what it is), then
// the array indices will be printed as they're tested.
//
int main(int argc, const char* const argv[])
{
const bool print_flag = (argc > 1);

double sum_RMS_diff = 0.0;

sum_RMS_diff += test_array1d_double(false, print_flag);
sum_RMS_diff += test_array1d_double(true,  print_flag);

sum_RMS_diff += test_array2d_float(false, print_flag);
sum_RMS_diff += test_array2d_float(true,  print_flag);

sum_RMS_diff += test_array3d_double(false, print_flag);
sum_RMS_diff += test_array3d_double(true,  print_flag);

#ifdef NOT_USED
sum_RMS_diff += test_array4d_double(false, print_flag);
sum_RMS_diff += test_array4d_double(true,  print_flag);
#endif

printf("\n");
printf("grand total RMS_diff = %g\n", sum_RMS_diff);
if (sum_RMS_diff == 0)
   then {
	printf("==> all ok!\n");
	return 0;
	}
   else {
	printf("==> ***** test(s) failed *****\n");
	return 1;
	}
}

//******************************************************************************

//
// This function tests the  array1d<double>  class.  It returns the
// "RMS difference" error measure (should be 0.0 if all is working ok).
//
double test_array1d_double(bool fancy_storage, bool print_flag)
{
const int min_i = 11;		const int max_i = 13;
const double delta = 1.0/7.0;
const double change_value = 144.025036146;

double *array;
array1d<double> *Aptr;
int stride;

// create the array
printf("\n");
printf("testing array1d<double>");
if (fancy_storage)
   then {
	printf(" with noncontiguous non-owned storage array...\n");
	const int N_i = jtutil::how_many_in_range(min_i, max_i);
	stride = 4;
	const int N = N_i*stride;
	printf(
   "allocating %d * stride=%d = %d underlying storage array...\n",
	       N_i, stride, N);
	array = new double[N];
	printf("constructing...\n");
	Aptr = new array1d<double>(min_i, max_i,
				   array,
				   stride);
	}
   else {
	printf(" with its own storage array...\n");
	Aptr = new array1d<double>(min_i, max_i);
	stride = 1;
	}
array1d<double>& A = *Aptr;

assert(A.min_i() == min_i);	assert(A.max_i() == max_i);
assert(A.N_i() == jtutil::how_many_in_range(min_i, max_i));

printf("assigning values...\n");
#define VALUE_1D(ix)	\
		double(1000.0 + ix + delta)
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	if (print_flag)
	   then printf("   assigning %d\n", i);
	A(i) = VALUE_1D(i);
	}
	  }

printf("checking N-D vs 1-D addresses...\n");
const double* array_ptr = A.data_array();
int sub = 0;
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	if (print_flag)
	   then printf("   checking %d\n", i);
	assert( &A(i) == &array_ptr[sub] );
	sub += stride;
	}
	  }
assert( A.N_array() == sub );
assert( A.subscript_stride_i() == &A(min_i+1) - &A(min_i) );

printf("checking what changes when we assign an array element...\n");
double sumsq_of_diff = 0.0;
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	if (print_flag)
	   then printf("   changing %d\n", i);
	double save_a = A(i);
	A(i) = change_value;

		for (int ii = min_i ; ii <= max_i ; ++ii)
		{
		if (print_flag)
		   then printf("      checking %d\n", ii);

		double should_be = (ii == i)
				   ? change_value : VALUE_1D(ii);
		assert( fuzzy<double>::EQ(A(ii), should_be) );

		double diff = A(ii) - should_be;
		sumsq_of_diff += diff*diff;
		}
	A(i) = save_a;
	}
	  }
double RMS_diff = sqrt(sumsq_of_diff / A.N_array());
printf("==> everything looks ok (RMS_diff=%.3g)\n", double(RMS_diff));

delete Aptr;

if (fancy_storage)
   then {
	printf("deleting underlying storage array...\n");
	delete[] array;
	}

return RMS_diff;
}

//******************************************************************************

//
// This function tests the  array2d<float>  class.  It returns the
// "RMS difference" error measure (should be 0.0 if all is working ok).
//
double test_array2d_float(bool fancy_storage, bool print_flag)
{
const int min_i = 11;		const int max_i = 13;
const int min_j = 21;		const int max_j = 24;
const float delta = 1.0/7.0;
const float change_value = 144.025036146;

float *array;
array2d<float> *Aptr;
int stride;

// create the array
printf("\n");
printf("testing array2d<float>");
if (fancy_storage)
   then {
	printf(" with noncontiguous non-owned storage array...\n");
	const int N_i = jtutil::how_many_in_range(min_i, max_i);
	const int N_j = jtutil::how_many_in_range(min_j, max_j);
	stride = 3;
	const int N = N_i*N_j*stride;
	printf(
   "allocating %d*%d * stride=%d = %d underlying storage array...\n",
	       N_i,N_j, stride, N);
	array = new float[N];
	printf("constructing...\n");
	Aptr = new array2d<float>(min_i, max_i,
				  min_j, max_j,
				  array,
				  N_j*stride, stride);
	}
   else {
	printf(" with its own storage array...\n");
	Aptr = new array2d<float>(min_i, max_i,
				  min_j, max_j);
	stride = 1;
	}
array2d<float>& A = *Aptr;

assert(A.min_i() == min_i);	assert(A.max_i() == max_i);
assert(A.min_j() == min_j);	assert(A.max_j() == max_j);
assert(A.N_i() == jtutil::how_many_in_range(min_i, max_i));
assert(A.N_j() == jtutil::how_many_in_range(min_j, max_j));

printf("assigning values...\n");
#define VALUE_2D(ix,jx)	\
		float(1000.0 + 10.0*jx + ix + delta)
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	for (int j = min_j ; j <= max_j ; ++j)
	{
	if (print_flag)
	   then printf("   assigning %d %d\n", i, j);
	A(i,j) = VALUE_2D(i,j);
	}
	}
	  }

printf("checking N-D vs 1-D addresses...\n");
const float* array_ptr = A.data_array();
int sub = 0;
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	for (int j = min_j ; j <= max_j ; ++j)
	{
	if (print_flag)
	   then printf("   checking %d %d\n", i, j);
	assert( &A(i,j) == &array_ptr[sub] );
	sub += stride;
	}
	}
	  }
assert( A.N_array() == sub );
assert( A.subscript_stride_i() == &A(min_i+1,min_j) - &A(min_i,min_j) );
assert( A.subscript_stride_j() == &A(min_i,min_j+1) - &A(min_i,min_j) );

printf("checking what changes when we assign an array element...\n");
double sumsq_of_diff = 0.0;
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	for (int j = min_j ; j <= max_j ; ++j)
	{
	if (print_flag)
	   then printf("   changing %d %d\n", i, j);
	float save_a = A(i,j);
	A(i,j) = change_value;

		for (int ii = min_i ; ii <= max_i ; ++ii)
		{
		for (int jj = min_j ; jj <= max_j ; ++jj)
		{
		if (print_flag)
		   then printf("      checking %d %d\n", ii, jj);

		float should_be = ((ii == i) && (jj == j))
				  ? change_value : VALUE_2D(ii,jj);
		assert( fuzzy<float>::EQ(A(ii,jj), should_be) );

		double diff = A(ii,jj) - should_be;
		sumsq_of_diff += diff*diff;
		}
		}
	A(i,j) = save_a;
	}
	}
	  }
double RMS_diff = sqrt(sumsq_of_diff / A.N_array());
printf("==> everything looks ok (RMS_diff=%.3g)\n", double(RMS_diff));

delete Aptr;

if (fancy_storage)
   then {
	printf("deleting underlying storage array...\n");
	delete[] array;
	}

return RMS_diff;
}

//******************************************************************************

//
// This function tests the  array3d<double>  class.
//
double test_array3d_double(bool fancy_storage, bool print_flag)
{
const int min_i = 11;		const int max_i = 13;
const int min_j = 21;		const int max_j = 24;
const int min_k = 31;		const int max_k = 35;
const double delta = 1.0/7.0;
const double change_value = 144.025036146;

double *array;
array3d<double> *Aptr;
int stride;

// create the array
printf("\n");
printf("testing array3d<double>");
if (fancy_storage)
   then {
	printf(" with noncontiguous non-owned storage array...\n");
	const int N_i = jtutil::how_many_in_range(min_i, max_i);
	const int N_j = jtutil::how_many_in_range(min_j, max_j);
	const int N_k = jtutil::how_many_in_range(min_k, max_k);
	stride = 4;
	const int N = N_i*N_j*N_k*stride;
	printf(
   "allocating %d*%d*%d * stride=%d = %d   underlying storage array...\n",
	       N_i,N_j,N_k, stride, N);
	array = new double[N];
	printf("constructing...\n");
	Aptr = new array3d<double>(min_i, max_i,
				   min_j, max_j,
				   min_k, max_k,
				   array,
				   N_j*N_k*stride,N_k*stride, stride);
	}
   else {
	printf(" with its own storage array...\n");
	Aptr = new array3d<double>(min_i, max_i,
				   min_j, max_j,
				   min_k, max_k);
	stride = 1;
	}
array3d<double>& A = *Aptr;

assert(A.min_i() == min_i);	assert(A.max_i() == max_i);
assert(A.min_j() == min_j);	assert(A.max_j() == max_j);
assert(A.min_k() == min_k);	assert(A.max_k() == max_k);
assert(A.N_i() == jtutil::how_many_in_range(min_i, max_i));
assert(A.N_j() == jtutil::how_many_in_range(min_j, max_j));
assert(A.N_k() == jtutil::how_many_in_range(min_k, max_k));

printf("assigning values...\n");
#define VALUE_3D(ix,jx,kx)	\
		double(1000.0 + 100.0*kx + 10.0*jx + ix + delta)
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	for (int j = min_j ; j <= max_j ; ++j)
	{
	for (int k = min_k ; k <= max_k ; ++k)
	{
	if (print_flag)
	   then printf("   assigning %d %d %d\n", i, j, k);
	A(i,j,k) = VALUE_3D(i,j,k);
	}
	}
	}
	  }

printf("checking N-D vs 1-D addresses...\n");
const double* array_ptr = A.data_array();
int sub = 0;
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	for (int j = min_j ; j <= max_j ; ++j)
	{
	for (int k = min_k ; k <= max_k ; ++k)
	{
	if (print_flag)
	   then printf("   checking %d %d %d\n", i, j, k);
	assert( &A(i,j,k) == &array_ptr[sub] );
	sub += stride;
	}
	}
	}
	  }
assert(sub == A.N_array());
assert( A.subscript_stride_i()
	== &A(min_i+1,min_j,min_k) - &A(min_i,min_j,min_k) );
assert( A.subscript_stride_j()
	== &A(min_i,min_j+1,min_k) - &A(min_i,min_j,min_k) );
assert( A.subscript_stride_k()
	== &A(min_i,min_j,min_k+1) - &A(min_i,min_j,min_k) );

printf("checking which things change when we assign an array element...\n");
double sumsq_of_diff = 0.0;
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	for (int j = min_j ; j <= max_j ; ++j)
	{
	for (int k = min_k ; k <= max_k ; ++k)
	{
	if (print_flag)
	   then printf("   changing %d %d %d\n", i, j, k);
	double save_a = A(i,j,k);
	A(i,j,k) = change_value;

		for (int ii = min_i ; ii <= max_i ; ++ii)
		{
		for (int jj = min_j ; jj <= max_j ; ++jj)
		{
		for (int kk = min_k ; kk <= max_k ; ++kk)
		{
		if (print_flag)
		   then printf("      checking %d %d %d\n", ii, jj, kk);

		double should_be = ((ii == i) && (jj == j) && (kk == k))
				   ? change_value : VALUE_3D(ii,jj,kk);
		assert( fuzzy<double>::EQ(A(ii,jj,kk), should_be) );

		double diff = A(ii,jj,kk) - should_be;
		sumsq_of_diff += diff*diff;
		}
		}
		}
	A(i,j,k) = save_a;
	}
	}
	}
	  }
double RMS_diff = sqrt(sumsq_of_diff / A.N_array());
printf("==> everything looks ok (RMS_diff=%.3g)\n", RMS_diff);

if (fancy_storage)
   then {
	printf("deleting underlying storage array...\n");
	delete[] array;
	}

return RMS_diff;
}

//******************************************************************************

#ifdef NOT_USED
//
// This function tests the  array4d<double>  class.
//
double test_array4d_double(bool fancy_storage, bool print_flag)
{
const int min_i = 11;		const int max_i = 13;
const int min_j = 21;		const int max_j = 24;
const int min_k = 31;		const int max_k = 35;
const int min_l = 41;		const int max_l = 46;
const double delta = 1.0/7.0;
const double change_value = 144.025036146;

double *array;
array4d<double> *Aptr;
int stride;

// create the array
printf("\n");
printf("testing array4d<double>");
if (fancy_storage)
   then {
	printf(" with noncontiguous non-owned storage array...\n");
	const int N_i = jtutil::how_many_in_range(min_i, max_i);
	const int N_j = jtutil::how_many_in_range(min_j, max_j);
	const int N_k = jtutil::how_many_in_range(min_k, max_k);
	const int N_l = jtutil::how_many_in_range(min_l, max_l);
	stride = 3;
	const int N = N_i*N_j*N_k*N_l*stride;
	printf(
   "allocating %d*%d*%d*%d * stride=%d = %d   underlying storage array...\n",
	       N_i,N_j,N_k,N_l, stride, N);
	array = new double[N];
	printf("constructing...\n");
	Aptr = new array4d<double>(min_i, max_i,
				   min_j, max_j,
				   min_k, max_k,
				   min_l, max_l,
				   array,
				   N_j*N_k*N_l*stride, N_k*N_l*stride,
				   N_l*stride, stride);
	}
   else {
	printf(" with its own storage array...\n");
	Aptr = new array4d<double>(min_i, max_i,
				   min_j, max_j,
				   min_k, max_k,
				   min_l, max_l);
	stride = 1;
	}
array4d<double>& A = *Aptr;

assert(A.min_i() == min_i);	assert(A.max_i() == max_i);
assert(A.min_j() == min_j);	assert(A.max_j() == max_j);
assert(A.min_k() == min_k);	assert(A.max_k() == max_k);
assert(A.min_l() == min_l);	assert(A.max_l() == max_l);
assert(A.N_i() == jtutil::how_many_in_range(min_i, max_i));
assert(A.N_j() == jtutil::how_many_in_range(min_j, max_j));
assert(A.N_k() == jtutil::how_many_in_range(min_k, max_k));
assert(A.N_l() == jtutil::how_many_in_range(min_l, max_l));

printf("assigning values...\n");
#define VALUE_4D(ix,jx,kx,lx)	\
		double(10000.0 + 1000.0*kx + 100.0*jx + 10.0*ix + lx + delta)
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	for (int j = min_j ; j <= max_j ; ++j)
	{
	for (int k = min_k ; k <= max_k ; ++k)
	{
	for (int l = min_l ; l <= max_l ; ++l)
	{
	if (print_flag)
	   then printf("   assigning %d %d %d %d\n", i, j, k, l);
	A(i,j,k,l) = VALUE_4D(i,j,k,l);
	}
	}
	}
	}
	  }

printf("checking N-D vs 1-D addresses...\n");
const double* array_ptr = A.data_array();
int sub = 0;
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	for (int j = min_j ; j <= max_j ; ++j)
	{
	for (int k = min_k ; k <= max_k ; ++k)
	{
	for (int l = min_l ; l <= max_l ; ++l)
	{
	if (print_flag)
	   then printf("   checking %d %d %d %d\n", i, j, k, l);
	assert( &A(i,j,k,l) == &array_ptr[sub] );
	sub += stride;
	}
	}
	}
	}
	  }
assert(sub == A.N_array());
assert( A.subscript_stride_i()
	== &A(min_i+1,min_j,min_k,min_l) - &A(min_i,min_j,min_k,min_l) );
assert( A.subscript_stride_j()
	== &A(min_i,min_j+1,min_k,min_l) - &A(min_i,min_j,min_k,min_l) );
assert( A.subscript_stride_k()
	== &A(min_i,min_j,min_k+1,min_l) - &A(min_i,min_j,min_k,min_l) );
assert( A.subscript_stride_l()
	== &A(min_i,min_j,min_k,min_l+1) - &A(min_i,min_j,min_k,min_l) );

printf("checking which things change when we assign an array element...\n");
double sumsq_of_diff = 0.0;
	  {
	for (int i = min_i ; i <= max_i ; ++i)
	{
	for (int j = min_j ; j <= max_j ; ++j)
	{
	for (int k = min_k ; k <= max_k ; ++k)
	{
	for (int l = min_l ; l <= max_l ; ++l)
	{
	if (print_flag)
	   then printf("   changing %d %d %d %d\n", i, j, k, l);
	double save_a = A(i,j,k,l);
	A(i,j,k,l) = change_value;

		for (int ii = min_i ; ii <= max_i ; ++ii)
		{
		for (int jj = min_j ; jj <= max_j ; ++jj)
		{
		for (int kk = min_k ; kk <= max_k ; ++kk)
		{
		for (int ll = min_l ; ll <= max_l ; ++ll)
		{
		if (print_flag)
		   then printf("      checking %d %d %d %d\n", ii, jj, kk, ll);

		double should_be = ((ii == i) && (jj == j)
				     && (kk == k) && (ll == l))
				   ? change_value : VALUE_4D(ii,jj,kk,ll);
		assert( fuzzy<double>::EQ(A(ii,jj,kk,ll), should_be) );

		double diff = A(ii,jj,kk,ll) - should_be;
		sumsq_of_diff += diff*diff;
		}
		}
		}
		}
	A(i,j,k,l) = save_a;
	}
	}
	}
	}
	  }
double RMS_diff = sqrt(sumsq_of_diff / A.N_array());
printf("==> everything looks ok (RMS_diff=%.3g)\n", RMS_diff);

if (fancy_storage)
   then {
	printf("deleting underlying storage array...\n");
	delete[] array;
	}

return RMS_diff;
}
#endif	/* NOT_USED */
