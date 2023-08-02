# Makefile to test various C/C++ compiler features
#
# Targets:
#	default		test namespace handling for <stdio.h> and <cstdio>
#	stdio_h		test namespace handling for <stdio.h>
#	cstdio		test namespace handling for <cstdio>
#	clean		remove all binaries
#
# Arguments:
# CXX		how to invoke the C++ compiler, including any desired flags
#

# default setting, may be overridden from command line
CXX := g++ -W -Wall -pedantic -ansi

.PHONY	: default
default	: stdio_h cstdio
	
.PHONY	: stdio_h
stdio_h	:
	-$(CXX) -o stdio_h-global \
                   stdio_h-global.cc
	        -./stdio_h-global
	-$(CXX) -o stdio_h-std--std-printf \
                   stdio_h-std--std-printf.cc
	        -./stdio_h-std--std-printf
	-$(CXX) -o stdio_h-std--using-std-printf \
                   stdio_h-std--using-std-printf.cc
	        -./stdio_h-std--using-std-printf
	-$(CXX) -o stdio_h-std--using-namespace-std \
                   stdio_h-std--using-namespace-std.cc
	        -./stdio_h-std--using-namespace-std

.PHONY	: cstdio
cstdio	:
	-$(CXX) -o cstdio-global \
                   cstdio-global.cc
	        -./cstdio-global
	-$(CXX) -o cstdio-std--std-printf \
                   cstdio-std--std-printf.cc
	        -./cstdio-std--std-printf
	-$(CXX) -o cstdio-std--using-std-printf \
                   cstdio-std--using-std-printf.cc
	        -./cstdio-std--using-std-printf
	-$(CXX) -o cstdio-std--using-namespace-std \
                   cstdio-std--using-namespace-std.cc
	        -./cstdio-std--using-namespace-std

.PHONY	: clean
clean	:
	-rm -f stdio_h-global \
               stdio_h-std--std-printf \
               stdio_h-std--using-std-printf \
               stdio_h-std--using-namespace-std
	-rm -f cstdio-global \
               cstdio-std--std-printf \
               cstdio-std--using-std-printf \
               cstdio-std--using-namespace-std
