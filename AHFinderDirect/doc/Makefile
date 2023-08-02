# Makefile for documentation
# $Header: /usr/local/svn/cvs-repositories/numrelcvs/AEIThorns/AHFinderDirect/doc/Makefile,v 1.7 2003-07-31 17:30:20 jthorn Exp $

#
# ***** main targets *****
# default	latex & dvips everything
# bib		bibtex everything
# clean		remove all machine-generated files
#
# ***** targets for rebuilding figures in the thorn guide *****
# 3patch-cactus		run Cactus to generate the data file for the
#			3patch figure
# 3patch-gnuplot	run gnuplot to generate the postscript file
#			for the 3patch figure
#

#
# Note this Makefile always re-TeXs things, even if the .tex files
# haven't changed (and similarly for dvips); this ensures that changes
# to macros, included stuff from other directories, etc, get "used".
# To implement this, we declare the .dvi and .ps files as "PHONY".
#

###############################################################################

#
# ***** main targets *****
#

.PHONY	: default
default	: ps

.PHONY	: ps
ps	: dvi
	dvips -f <documentation.dvi >documentation.ps

.PHONY	: dvi
dvi	:
	latex documentation

.PHONY	: bib
bib	:
	bibtex documentation

.PHONY	: clean
clean	:
	-rm *.dvi *.ps *.aux *.bbl *.log *.blg

###############################################################################

#
# ***** targets for rebuilding figures in the thorn guide *****
#

.PHONY		: 3patch-cactus
3patch-cactus	:
	cactus_test-ahfinderdirect 3patch.par 2>&1 | tee 3patch.log

.PHONY		: 3patch-gnuplot
3patch-gnuplot	:
	gnuplot 3patch.gnuplot
