# Makefile to test various C/C++ compiler features
#
# Targets:
#	default		test namespace handling for <vector.h> and <vector>
#	vector_h	test namespace handling for <vector.h>
#	vector		test namespace handling for <vector>
#	clean		remove all binaries
#
# Arguments:
# CXX		how to invoke the C++ compiler, including any desired flags
#

# default setting, may be overridden from command line
CXX := g++ -W -Wall -pedantic -ansi

.PHONY	: default
default	: vector_h vector
	
.PHONY	: vector_h
vector_h	:
	-$(CXX) -o vector_h-global \
                   vector_h-global.cc
	        -./vector_h-global
	-$(CXX) -o vector_h-std--std-vector \
                   vector_h-std--std-vector.cc
	        -./vector_h-std--std-vector
	-$(CXX) -o vector_h-std--using-std-vector \
                   vector_h-std--using-std-vector.cc
	        -./vector_h-std--using-std-vector
	-$(CXX) -o vector_h-std--using-namespace-std \
                   vector_h-std--using-namespace-std.cc
	        -./vector_h-std--using-namespace-std

.PHONY	: vector
vector	:
	-$(CXX) -o vector-global \
                   vector-global.cc
	        -./vector-global
	-$(CXX) -o vector-std--std-vector \
                   vector-std--std-vector.cc
	        -./vector-std--std-vector
	-$(CXX) -o vector-std--using-std-vector \
                   vector-std--using-std-vector.cc
	        -./vector-std--using-std-vector
	-$(CXX) -o vector-std--using-namespace-std \
                   vector-std--using-namespace-std.cc
	        -./vector-std--using-namespace-std

.PHONY	: clean
clean	:
	-rm -f vector_h-global \
               vector_h-std--std-vector \
               vector_h-std--using-std-vector \
               vector_h-std--using-namespace-std
	-rm -f vector-global \
               vector-std--std-vector \
               vector-std--using-std-vector \
               vector-std--using-namespace-std
