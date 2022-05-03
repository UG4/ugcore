	
TESTS = \
	boost_test0 \
	boost_test1 \
	boost_test2

# sparsematrixgraph_test
#	sm_test0

all: ${TESTS}

MPI_INCLUDE=-I/usr/lib/x86_64-linux-gnu/openmpi/include

boost_test0: CXXFLAGS=-std=c++11 -g -O0
boost_test0: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}

boost_test1 boost_test2: CXXFLAGS=-std=c++11 -g -O0
boost_test1 boost_test2: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}

sm_test0: CXXFLAGS=-std=c++11 -g -O0
sm_test0: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}

sparsematrixgraph_test: CXXFLAGS=-std=c++11 -g -O0
sparsematrixgraph_test: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}
