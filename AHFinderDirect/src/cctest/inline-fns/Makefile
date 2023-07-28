# Makefile to test various C/C++ compiler features
# $Header: /usr/local/svn/cvs-repositories/numrelcvs/AEIThorns/AHFinderDirect/src/cctest/inline-fns/Makefile,v 1.1 2003-07-31 17:29:48 jthorn Exp $
#
# Targets:
#	default		inline-fn inline-class
#	inline-fn	test non-class inline function in .cc file
#	inline-fn2	test non-class inline function in header file
#	inline-class	test class with inline function
#	clean		remove all binaries
#
# Arguments:
# CXX		how to invoke the C++ compiler, including any desired flags
#

# default setting, may be overridden from command line
CXX := g++ -W -Wall -pedantic -ansi

.PHONY	: default
default	: inline-fn inline-fn2 inline-class

.PHONY	: inline-fn
inline-fn:
	$(CXX) -o inline-fn inline-fn.cc

.PHONY	: inline-fn2
inline-fn2:
	$(CXX) -o inline-fn2 inline-fn2.cc

.PHONY	: inline-class
inline-class:
	$(CXX) -o inline-class inline-class.cc

.PHONY	: clean
clean	:
	-rm inline-fn inline-fn2 inline-class
