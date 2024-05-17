// test_horizon_sequence.cc -- test driver for horizon_sequence class
// $Header$

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "../include/stdc.h"
#include "horizon_sequence.hh"

using namespace AHFinderDirect;

// prototypes for functions local to this file
namespace {
void test0();
void test1();
void test2();
void test345();
	  }

int main(void)
{
test0();
test1();
test2();
test345();

printf("all ok!\n");
return 0;
}

//******************************************************************************

namespace {
void test0()
{
printf("capacity 0 ...\n");
horizon_sequence hs0(0);
assert( hs0.N_horizons() == 0 );
assert( hs0.my_N_horizons() == 0 );
assert( ! hs0.has_genuine_horizons() );
assert( STRING_EQUAL(hs0.sequence_string(","), "") );
assert( hs0.is_dummy() );
assert( ! hs0.is_genuine() );
assert( ! hs0.is_next_genuine() );
assert( hs0.dummy_number() == 1 );
assert( hs0.init_hn() == 0 );
	assert( hs0.get_hn() == 0 );
	assert( hs0.dummy_number() == 1 );
assert( hs0.next_hn() == 0 );
	assert( hs0.get_hn() == 0 );
	assert( hs0.dummy_number() == 2 );
assert( hs0.next_hn() == 0 );
	assert( hs0.get_hn() == 0 );
	assert( hs0.dummy_number() == 3 );
assert( hs0.next_hn() == 0 );
	assert( hs0.get_hn() == 0 );
	assert( hs0.dummy_number() == 4 );
assert( hs0.N_horizons() == 0 );
assert( hs0.my_N_horizons() == 0 );
assert( ! hs0.is_hn_genuine(-1) );
assert( ! hs0.is_hn_genuine(0) );
assert( ! hs0.is_hn_genuine(1) );
assert( ! hs0.is_hn_genuine(42) );
}
	  }

//******************************************************************************

namespace {
void test1()
{
printf("capacity 1 ...\n");
horizon_sequence hs1(1);
assert( hs1.N_horizons() == 1 );
assert( hs1.my_N_horizons() == 0 );
assert( ! hs1.has_genuine_horizons() );
assert( ! hs1.is_hn_genuine(-1) );
assert( ! hs1.is_hn_genuine(0) );
assert( ! hs1.is_hn_genuine(1) );
assert( ! hs1.is_hn_genuine(42) );
assert( hs1.append_hn(42) == 1 );
assert( hs1.N_horizons() == 1 );
assert( hs1.my_N_horizons() == 1 );
assert( hs1.has_genuine_horizons() );
assert( STRING_EQUAL(hs1.sequence_string(","), "42") );
assert( hs1.is_genuine() );
assert( ! hs1.is_next_genuine() );
assert( hs1.dummy_number() == 0 );
assert( ! hs1.is_hn_genuine(-1) );
assert( ! hs1.is_hn_genuine(0) );
assert( ! hs1.is_hn_genuine(1) );
assert(   hs1.is_hn_genuine(42) );
assert( hs1.init_hn() == 42 );
	assert( hs1.get_hn() == 42 );
	assert( hs1.is_genuine() );
	assert( ! hs1.is_next_genuine() );
	assert( hs1.dummy_number() == 0 );
assert( hs1.next_hn() == 0 );
	assert( hs1.get_hn() == 0 );
	assert( ! hs1.is_genuine() );
	assert( ! hs1.is_next_genuine() );
	assert( hs1.dummy_number() == 1 );
assert( hs1.next_hn() == 0 );
	assert( hs1.get_hn() == 0 );
	assert( ! hs1.is_genuine() );
	assert( ! hs1.is_next_genuine() );
	assert( hs1.dummy_number() == 2 );
assert( hs1.next_hn() == 0 );
	assert( hs1.get_hn() == 0 );
	assert( ! hs1.is_genuine() );
	assert( ! hs1.is_next_genuine() );
	assert( hs1.dummy_number() == 3 );
assert( hs1.next_hn() == 0 );
	assert( hs1.get_hn() == 0 );
	assert( ! hs1.is_genuine() );
	assert( ! hs1.is_next_genuine() );
	assert( hs1.dummy_number() == 4 );
assert( hs1.N_horizons() == 1 );
assert( hs1.my_N_horizons() == 1 );
}
	  }

//******************************************************************************

namespace {
void test2()
{
printf("capacity 2 ...\n");
horizon_sequence hs2(2);
assert( hs2.N_horizons() == 2 );
assert( hs2.my_N_horizons() == 0 );
assert( ! hs2.has_genuine_horizons() );
assert( ! hs2.is_hn_genuine(-1) );
assert( ! hs2.is_hn_genuine(0) );
assert( ! hs2.is_hn_genuine(1) );
assert( ! hs2.is_hn_genuine(42) );
assert( ! hs2.is_hn_genuine(69) );
assert( hs2.append_hn(42) == 1 );
assert( hs2.append_hn(69) == 2 );
assert( hs2.N_horizons() == 2 );
assert( hs2.my_N_horizons() == 2 );
assert( hs2.has_genuine_horizons() );
assert( ! hs2.is_hn_genuine(-1) );
assert( ! hs2.is_hn_genuine(0) );
assert( ! hs2.is_hn_genuine(1) );
assert(   hs2.is_hn_genuine(42) );
assert(   hs2.is_hn_genuine(69) );
assert( ! hs2.is_hn_genuine(1000) );
assert( STRING_EQUAL(hs2.sequence_string(","), "42,69") );
assert( hs2.is_genuine() );
assert( hs2.is_next_genuine() );
assert( hs2.dummy_number() == 0 );
assert( hs2.init_hn() == 42 );
	assert( hs2.get_hn() == 42 );
	assert( hs2.is_genuine() );
	assert( hs2.is_next_genuine() );
	assert( hs2.dummy_number() == 0 );
assert( ! hs2.is_hn_genuine(-1) );
assert( ! hs2.is_hn_genuine(0) );
assert( ! hs2.is_hn_genuine(1) );
assert(   hs2.is_hn_genuine(42) );
assert(   hs2.is_hn_genuine(69) );
assert( ! hs2.is_hn_genuine(1000) );
assert( hs2.next_hn() == 69 );
	assert( hs2.get_hn() == 69 );
	assert( hs2.is_genuine() );
	assert( ! hs2.is_next_genuine() );
	assert( hs2.dummy_number() == 0 );
assert( hs2.next_hn() == 0 );
	assert( hs2.get_hn() == 0 );
	assert( ! hs2.is_genuine() );
	assert( ! hs2.is_next_genuine() );
	assert( hs2.dummy_number() == 1 );
assert( hs2.next_hn() == 0 );
	assert( hs2.get_hn() == 0 );
	assert( ! hs2.is_genuine() );
	assert( ! hs2.is_next_genuine() );
	assert( hs2.dummy_number() == 2 );
assert( hs2.next_hn() == 0 );
	assert( hs2.get_hn() == 0 );
	assert( ! hs2.is_genuine() );
	assert( ! hs2.is_next_genuine() );
	assert( hs2.dummy_number() == 3 );
assert( hs2.N_horizons() == 2 );
assert( hs2.my_N_horizons() == 2 );
assert( ! hs2.is_hn_genuine(-1) );
assert( ! hs2.is_hn_genuine(0) );
assert( ! hs2.is_hn_genuine(1) );
assert(   hs2.is_hn_genuine(42) );
assert(   hs2.is_hn_genuine(69) );
assert( ! hs2.is_hn_genuine(1000) );
}
	  }

//******************************************************************************

namespace {
void test345()
{
	for (int capacity = 3 ; capacity <= 5 ; ++capacity)
	{
	printf("capacity %d ...\n", capacity);
	horizon_sequence hs(capacity);
	assert( hs.N_horizons() == capacity );
	assert( hs.my_N_horizons() == 0 );
	assert( ! hs.has_genuine_horizons() );

	assert( hs.init_hn() == 0 );
	assert( hs.dummy_number() == 1 );

	assert( hs.next_hn() == 0 );
	assert( hs.dummy_number() == 2 );

	assert( hs.next_hn() == 0 );
	assert( hs.dummy_number() == 3 );

	assert( ! hs.is_hn_genuine(-1) );
	assert( ! hs.is_hn_genuine(0) );
	assert( ! hs.is_hn_genuine(1) );
	assert( ! hs.is_hn_genuine(42) );
	assert( ! hs.is_hn_genuine(69) );
	assert( ! hs.is_hn_genuine(1000) );

	assert( hs.append_hn(42) == 1 );
	assert( hs.append_hn(69) == 2 );
	assert( hs.append_hn(105) == 3 );
	assert( hs.N_horizons() == capacity );
	assert( hs.my_N_horizons() == 3 );
	assert( hs.has_genuine_horizons() );
	assert( hs.is_genuine() );
	assert( hs.is_next_genuine() );
	assert( hs.dummy_number() == 0 );

	assert( ! hs.is_hn_genuine(-1) );
	assert( ! hs.is_hn_genuine(0) );
	assert( ! hs.is_hn_genuine(1) );
	assert(   hs.is_hn_genuine(42) );
	assert(   hs.is_hn_genuine(69) );
	assert( ! hs.is_hn_genuine(100) );
	assert(   hs.is_hn_genuine(105) );
	assert( ! hs.is_hn_genuine(1000) );

	assert( STRING_EQUAL(hs.sequence_string(","), "42,69,105") );

		for (int i = 1 ; i <= 4 ; ++i)
		{
		printf("   try %d...\n", i);
		assert( hs.N_horizons() == capacity );
		assert( hs.my_N_horizons() == 3 );
		assert( hs.init_hn() == 42 );
			assert( hs.get_hn() == 42 );
			assert( hs.is_genuine() );
			assert( hs.is_next_genuine() );
			assert( hs.dummy_number() == 0 );

		assert( ! hs.is_hn_genuine(-1) );
		assert( ! hs.is_hn_genuine(0) );
		assert( ! hs.is_hn_genuine(1) );
		assert(   hs.is_hn_genuine(42) );
		assert(   hs.is_hn_genuine(69) );
		assert( ! hs.is_hn_genuine(100) );
		assert(   hs.is_hn_genuine(105) );
		assert( ! hs.is_hn_genuine(1000) );

		assert( hs.next_hn() == 69 );
			assert( hs.get_hn() == 69 );
			assert( hs.is_genuine() );
			assert( hs.is_next_genuine() );
			assert( hs.dummy_number() == 0 );
		assert( hs.next_hn() == 105 );
			assert( hs.get_hn() == 105 );
			assert( hs.is_genuine() );
			assert( ! hs.is_next_genuine() );
			assert( hs.dummy_number() == 0 );
		assert( hs.next_hn() == 0 );
			assert( hs.get_hn() == 0 );
			assert( ! hs.is_genuine() );
			assert( ! hs.is_next_genuine() );
			assert( hs.dummy_number() == 1 );
		assert( hs.next_hn() == 0 );
			assert( hs.get_hn() == 0 );
			assert( ! hs.is_genuine() );
			assert( ! hs.is_next_genuine() );
			assert( hs.dummy_number() == 2 );
		assert( hs.next_hn() == 0 );
			assert( hs.get_hn() == 0 );
			assert( ! hs.is_genuine() );
			assert( ! hs.is_next_genuine() );
			assert( hs.dummy_number() == 3 );

		assert( ! hs.is_hn_genuine(-1) );
		assert( ! hs.is_hn_genuine(0) );
		assert( ! hs.is_hn_genuine(1) );
		assert(   hs.is_hn_genuine(42) );
		assert(   hs.is_hn_genuine(69) );
		assert( ! hs.is_hn_genuine(100) );
		assert(   hs.is_hn_genuine(105) );
		assert( ! hs.is_hn_genuine(1000) );
		}
	}
}
	  }
