# Makefile for maple code
# $Header: /usr/local/svn/cvs-repositories/numrelcvs/AEIThorns/AHFinderDirect/src/maple/Makefile,v 1.4.2.1 2004-04-06 21:04:22 schnetter Exp $

#
# Environment Variables:
#	MAPLE_VERSION	used via @ifdef for version control in Maple code;
#			typically set to something like  MAPLE_V_RELEASE_4
#			(not presently used, but may be needed in the future)
#
# Targets:
#	mm		==> preprocess all *.maple files to produce *.mm files
#	clean		==> delete *.mm (mpp Maple preprocessor output)
#	run		==> load all the code into Maple
#

ifneq ($(MAPLE_VERSION),)
MPP_FLAGS := -D$(MAPLE_VERSION)
endif

###############################################################################

.PHONY	: mm
mm	: $(patsubst %.maple, %.mm, $(wildcard *.maple))

%.mm	: %.maple $(wildcard *.minc $(gfa_dir)/*.minc $(gfa_dir)/*.maple)
	../misc/mpp $(MPP_FLAGS) <$< >$@

.PHONY	: clean
clean	:
	-rm *.mm

.PHONY	: run
run	: mm
	maple <setup.mm >maple.log
