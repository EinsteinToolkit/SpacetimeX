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
# superclean	==> delete object files, test drivers
#
# Bugs:
# - The makefile in ../jtutil/ should be executed first, since this
#   makefile assumes all the ../jtutil/*.o files already exist
# - Dependencies on *.hh are omitted.
#

CXX      := g++
CFLAGS   := $(STD_GCC_FLAGS) -I../include -DSTANDALONE_TEST -g
CXXFLAGS := $(STD_GXX_FLAGS) -I../include -DSTANDALONE_TEST -g

JTUTIL_OBJ	:= $(filter-out ../jtutil/test%, $(wildcard ../jtutil/*.o))

################################################################################

test		: test_coords test_coords2

test_coords	: test_coords.cc coords.cc $(JTUTIL_OBJ) -lm
test_coords2	: test_coords2.cc coords.cc $(JTUTIL_OBJ) -lm

.PHONY		: coord.derivs
coord.derivs	:
	maple <coord_derivs.maple >coord_derivs.out

.PHONY		: clean
clean		:
	-rm *.o test_coords test_coords2
