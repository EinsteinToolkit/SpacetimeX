# Makefile for standalone test drivers in this directory
# $Header$
#
#	CC, CXX = C and C++ compilers.  Defaults are gcc and g++ if
#		  variables aren't set from command line or environment.
#	CC_FLAGS, CXX_FLAGS = C++ compiler flags.  Defaults are
#				$(STD_GCC_FLAGS) -I.. -g
#			      and
#				$(STD_GXX_FLAGS) -I.. -g
#			      if variables aren't set from command line or
#			      environment.
#
# Targets:
# test		==> build test programs
# clean		==> delete object files, test drivers
#
# Bugs:
# - Dependencies on *.hh are omitted.
#

CXX      := g++
CFLAGS   := $(STD_GCC_FLAGS) -I../include -DSTANDALONE_TEST -g
CXXFLAGS := $(STD_GXX_FLAGS) -I../include -DSTANDALONE_TEST -g

ALL_TESTS := test_array test_array2 \
             test_cpm_map test_linear_map \
             test_fuzzy test_round \
             test_modulo test_norm \
             test_error_exit

################################################################################

test	: $(ALL_TESTS)

test_array	: test_array.o array.o \
                  fuzzy.o round.o -lm
	$(CXX) $(CXXFLAGS) -o $@ $?
test_array2	: test_array2.o array.o
	$(CXX) $(CXXFLAGS) -o $@ $?
test_cpm_map	: test_cpm_map.o cpm_map.o \
                  fuzzy.o round.o error_exit.o -lm
	$(CXX) $(CXXFLAGS) -o $@ $?
test_linear_map	: test_linear_map.o linear_map.o \
                  fuzzy.o round.o error_exit.o -lm
	$(CXX) $(CXXFLAGS) -o $@ $?
test_fuzzy	: test_fuzzy.o fuzzy.o \
                  round.o error_exit.o
	$(CXX) $(CXXFLAGS) -o $@ $?
test_round	: test_round.o round.o
	$(CXX) $(CXXFLAGS) -o $@ $?
test_modulo	: test_modulo.o miscfp.o \
                  fuzzy.o round.o error_exit.o -lm
	$(CXX) $(CXXFLAGS) -o $@ $?
test_norm	: test_norm.o norm.o \
                  fuzzy.o round.o -lm
	$(CXX) $(CXXFLAGS) -o $@ $?
test_error_exit	: test_error_exit.o error_exit.o
	$(CXX) $(CXXFLAGS) -o $@ $?

.PHONY	: clean
clean	:
	-rm -f *.o
	-rm -f $(ALL_TESTS)
