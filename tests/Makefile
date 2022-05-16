
TESTS = \
	boost_ptest0 \
	boost_ptest1 \
	boost_test0 \
	boost_test1 \
	boost_test2 \
	boost_test3 \
	boost_test4

# sparsematrixgraph_test
#	sm_test0

all: ${TESTS}

MPI_INCLUDE=-I/usr/lib/x86_64-linux-gnu/openmpi/include

${TESTS}: CXXFLAGS=-std=c++11 -g -O0 -Wall
${TESTS}: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}

sm_test0: CXXFLAGS=-std=c++11 -g -O0 -Wall
sm_test0: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}

sparsematrixgraph_test: CXXFLAGS=-std=c++11 -g -O0 -Wall
sparsematrixgraph_test: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}

LIBS =
boost_ptest0 boost_ptest1: LIBS = -lpcl_common -lmpi_cxx -lmpi
boost_ptest0 boost_ptest1: CXX = mpiCC
# boost_ptest0: CPPFLAGS=-I../ugbase ${MPI_INCLUDE} -DUG_PARALLEL
boost_ptest0 boost_ptest1: %: %.o
	${CXX} -o $@ $< ${LIBS}
