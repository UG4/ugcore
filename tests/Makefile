
PTESTS = \
	boost_ptest0 \
	boost_ptest1 \
	boost_ptest3

TESTS = \
	${PTESTS} \
	sm_transpose \
	boost_test0 \
	boost_test1 \
	boost_test3 \
	boost_test4

TEST_OUT = ${TESTS:%=out/%.out}

# sparsematrixgraph_test
#	sm_test0

all: ${TESTS}

out: ${TEST_OUT}

out/%.out: %
	mkdir -p out
	./$* | grep -v "\ refresh\ " > $@

MPI_INCLUDE=-I/usr/lib/x86_64-linux-gnu/openmpi/include

${TESTS}: CXXFLAGS=-std=c++11 -g -O0 -Wall
${TESTS}: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}

sm_test0: CXXFLAGS=-std=c++11 -g -O0 -Wall
sm_test0: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}

sparsematrixgraph_test: CXXFLAGS=-std=c++11 -g -O0 -Wall
sparsematrixgraph_test: CPPFLAGS=-I../ugbase ${MPI_INCLUDE}

LIBS =
${PTESTS}: LIBS = -lpcl_common -lmpi_cxx -lmpi
${PTESTS}: CXX = mpiCC
# boost_ptest0: CPPFLAGS=-I../ugbase ${MPI_INCLUDE} -DUG_PARALLEL
${PTESTS}: %: %.o
	${CXX} -o $@ $< ${LIBS}

clean:
	rm -rf *~ ${TESTS} out *.vtu
