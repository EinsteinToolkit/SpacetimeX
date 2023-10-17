# Makefile for driver code
# $Header$

#
# Targets:
#	ellipsoid	==> compute ellipsoid equations
#	test		==> build standalone test drivers
#	clean		==> remove standalone test drivers
#

CXXFLAGS := -g $(STD_GXX_FLAGS)

###############################################################################

# we don't want a (non--no-op) default target in this makefile
.PHONY		: dummy
dummy		:
	@echo 'there is no default target for this makefile'

.PHONY		: ellipsoid
ellipsoid	:
	maple <ellipsoid.maple >ellipsoid.log

.PHONY		: test
test		: test_horizon_sequence

.PHONY			: test_horizon_sequence
test_horizon_sequence	:
	$(CC) $(CFLAGS) -c \
               -DSTANDALONE_TEST -I../../../../../src/include \
               ../../../../../src/util/String.c \
               ../../../../../src/util/snprintf.c
	$(CXX) $(CXXFLAGS) -o $@ \
               -I../../../../../src/include \
               test_horizon_sequence.cc horizon_sequence.cc \
               String.o snprintf.o
.PHONY		: clean
clean		:
	-rm *.o test_horizon_sequence
