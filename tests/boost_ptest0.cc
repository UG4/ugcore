#define UG_PARALLEL

#include "lib_algebra/graph_interface/sparsematrix_boost.h"
#include "lib_algebra/cpu_algebra/sparsematrix_impl.h"

#include "lib_algebra/graph_interface/parallel_matrix_boost.h"
#include "pcl/pcl_base.cpp"
#include "pcl/pcl_util.cpp"

// library?
#include "common/log.cpp"
#include "common/util/file_util.cpp"
#include "common/debug_id.cpp"
#include "common/assert.cpp"
#include "common/util/crc32.cpp"
#include "common/util/ostream_buffer_splitter.cpp"
#include "common/util/string_util.cpp"
#include "common/util/binary_buffer.cpp"
#include "common/util/os_dependent_impl/file_util_posix.cpp"
#include "common/util/os_dependent_impl/os_info_linux.cpp"

// #include "lib_algebra/parallelization/parallel_nodes.cpp" // ?

#include "pcl/pcl_process_communicator.cpp"
#include "pcl/pcl_comm_world.cpp"

#include "boost/graph/graph_utility.hpp"

#include "lib_algebra/parallelization/algebra_layouts.cpp" // operator<<
#include "lib_algebra/parallelization/parallel_index_layout.cpp"

// #include "ug.cpp"

// typedef ug::SparseMatrix<ug::DenseMatrix<ug::FixedArray2<double, 3, 3> > > T; test1?
typedef ug::SparseMatrix<double> T;
typedef ug::ParallelMatrix<T> M;

static const int N=4;
static const int threads=2;

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	M A;
	A.resize_and_clear(N, N);

	// if (process_id(g.process_group()) == 0)
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		// Only process 0 loads the graph, which is distributed automatically
		A(0, 1) = 1.;
		//for (int i=1; i<n; ++i){
		//	std::cout << i << " g " << im[boost::vertex(i,g)] << " owner: " << owner(boost::vertex(i,g)) << "\n";
		//	add_edge(boost::vertex(i, g), boost::vertex(i-1, g), g);
		//	add_edge(boost::vertex(i-1, g), boost::vertex(i, g), g);
		//}
	}else{
	}
	auto L = A.layouts();
	assert(L);

	std::cout << "rank: " << rank << " " << *L << "\n";

	boost::graph_traits<M>::adjacency_iterator i;

	auto p = boost::adjacent_vertices(size_t(1), A);
	for(; p.first!=p.second; ++p.first){
		std::cout << " " << *p.first;
	}

	std::cout << "pg\n";
	boost::print_graph(A);

	{
		auto e = boost::vertices(A);
		for(; e.first!=e.second; ++e.first){
			auto v = *e.first;
			std::cout << v << " deg: " << boost::out_degree(v, A) << "\n";
		}
	}
	MPI_Finalize();
	return 0;
}

BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<ug::ParallelMatrix<ug::SparseMatrix<double>>> ));
BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<ug::ParallelMatrix<ug::SparseMatrix<double>>> ));
BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<ug::ParallelMatrix<ug::SparseMatrix<double>>> ));
