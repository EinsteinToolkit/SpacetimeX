# Makefile for GR code
# $Header$

#
# Environment Variables:
#	MAPLE_VERSION	used via @ifdef for version control in Maple code;
#			typically set to something like  MAPLE_V_RELEASE_4
#			(not presently used, but may be needed in the future)
#
# Targets:
#	cg		==> run all the *.mm files to generate C code 
#	mm		==> preprocess all *.maple files to produce *.mm files
#	clean		==> delete *.mm (mpp Maple preprocessor output)
#

ifneq ($(MAPLE_VERSION),)
MPP_FLAGS := -D$(MAPLE_VERSION)
endif

###############################################################################

# we don't want a (non--no-op) default target in this makefile
.PHONY	: dummy
dummy	:
	@echo 'there is no default target for this makefile'

.PHONY	: cg
cg	: mm
	maple <doit.mm 2>&1 | tee maple.log

.PHONY	: mm
mm	: $(patsubst %.maple, %.mm, $(wildcard *.maple))

%.mm	: %.maple $(wildcard *.minc $(gfa_dir)/*.minc $(gfa_dir)/*.maple)
	../misc/mpp $(MPP_FLAGS) <$< >$@

.PHONY	: clean
clean	:
	-rm *.mm
