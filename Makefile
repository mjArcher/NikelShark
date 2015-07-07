#Matt Archer. Makefile for Solid solved with riemann method v0.1 Feb 2014
#   2 # Can use $(shell ...) which expands to the output of the argument(s))
CUDA:=
CXX=g++
CPPFLAGS= -g -fopenmp -Wall -std=c++0x #- $(shell root-config --cflags)
LDFLAGS=-g #$(shell root-config --ldflags)
LDLIBS=-I ~/install/Libraries/ -I ~/install/Libraries/eigen/ -I ~/install/Libraries/eigen-dev/unsupported/ -I ~/install/Libraries/ -I ~/install/Libraries/libconfig/lib/ -lconfig++
RM=rm -f
# LDLIBS=-I ~/Libraries/eigen/ -I ~/Libraries/ -lconfig++ 

SRCS=ElasticState.cpp ElasticPrimState.cpp SolidSystem.cpp ElasticEOS.cpp Utils.cpp SquareTensor3.cpp
OBJS=$(subst .cpp,.o,$(SRCS))#substitutes all .cpp in SRCS variable with .o (objects)
LDADD=#~/Libraries/lib/libconfig++.la
LIBTOOL=~/Libraries/libconfig/libtool
CXXLINK=$(CXX)# $(LIBTOOL)

Elastic1D: Elastic1D.o $(OBJS) $(LDADD)
	rm -f Elastic1D
	$(CXXLINK) $(CPPFLAGS) $(LDFLAGS) -o Elastic1D Elastic1D.o $(OBJS) $(LDLIBS) $(LDADD)

EOSUnitTests: EOSUnitTests.o $(OBJS) $(LDADD)
	rm -f EOSUnitTests
	$(CXXLINK) $(CPPFLAGS) $(LDFLAGS) -o EOSUnitTests EOSUnitTests.o $(OBJS) $(LDLIBS) $(LDADD)

Elastic1D.o: Elastic1D.cpp ElasticState.h ElasticPrimState.h SolidSystem.h ElasticEOS.h Utils.h SquareTensor3.h 
	$(CXX) $(CPPFLAGS) -c Elastic1D.cpp -o Elastic1D.o $(LDLIBS)

EOSUnitTests.o: EOSUnitTests.cpp ElasticState.h ElasticPrimState.h SolidSystem.h ElasticEOS.h SquareTensor3.h
	$(CXX) $(CPPFLAGS) -c EOSUnitTests.cpp -o EOSUnitTests.o $(LDLIBS)

ElasticState.o: ElasticState.h ElasticState.cpp
	$(CXX) $(CPPFLAGS) -c ElasticState.cpp -o ElasticState.o $(LDLIBS)

ElasticPrimState.o: ElasticPrimState.h ElasticPrimState.cpp
	$(CXX) $(CPPFLAGS) -c ElasticPrimState.cpp -o ElasticPrimState.o $(LDLIBS)

SolidSystem.o: SolidSystem.h SolidSystem.cpp ElasticPrimState.h ElasticState.h ElasticEOS.h Utils.h
	$(CXX) $(CPPFLAGS) -c SolidSystem.cpp -o SolidSystem.o $(LDLIBS)

ElasticEOS.o: ElasticEOS.h ElasticEOS.cpp
	$(CXX) $(CPPFLAGS) -c ElasticEOS.cpp -o ElasticEOS.o $(LDLIBS)

Utils.o: Utils.cpp Utils.h 
	$(CXX) $(CPPFLAGS) -c Utils.cpp -o Utils.o $(LDLIBS)

SquareTensor3.o: SquareTensor3.cpp SquareTensor3.h 
	$(CXX) $(CPPFLAGS) -c SquareTensor3.cpp -o SquareTensor3.o $(LDLIBS)


# Tensor3.o: Tensor3.cpp Tensor3.h
# 	$(CXX) $(CPPFLAGS) -c Tensor3.cpp -o Tensor3.o 

# Tensor4.o: Tensor4.cpp Tensor4.h
# 	$(CXX) $(CPPFLAGS) -c Tensor4.cpp -o Tensor4.o 

# SymmetricMatrix.o: SymmetricMatrix.cpp SymmetrixMatrix.h
# 	$(CXX) $(CPPFLAGS) -c SymmetricMatrix.cpp -o SymmetricMatrix.o

# SquareMatrix.o: SquareMatrix.cpp SquareMatrix.h
# 	$(CXX) $(CPPFLAGS) -c SquareMatrix.cpp -o SquareMatrix.o $(LDLIBS)

clean:
	$(RM) Elastic1D
	$(RM) EOSUnitTests
	$(RM) *.o

dist-clean: clean
	$(RM) tool
# Sort.o: Sort.cpp Sort.h
#     g++ $(CPPFLAGS) -c Sort.cpp

