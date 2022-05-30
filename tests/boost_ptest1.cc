
#define UG_PARALLEL

#include "lib_algebra/graph_interface/sparsematrix_boost.h"
#include "lib_algebra/cpu_algebra/sparsematrix_impl.h"

#include "lib_algebra/graph_interface/parallel_matrix_boost.h"
#include "lib_algebra/graph_interface/bidirectional_boost.h"
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

static const int N=7;

typedef ug::SparseMatrix<double> T;
typedef ug::ParallelMatrix<T> PM;
typedef ug::BidirectionalMatrix<PM> BM;

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	PM M;

	typedef typename boost::graph_traits<BM>::in_edge_iterator in_edge_iterator;
	typedef typename boost::graph_traits<BM>::out_edge_iterator out_edge_iterator;

	M.resize_and_clear(N, N+2);

	M(0, 0) = 2.;
	M(1, N-1) = -1.;
	for(unsigned i=1; i<N; ++i){
		M(i, i) = 2.;
		M(i, i-1) = 1.;
		M(i-1, i) = 1.;
	}
	M(N-1, N+1) = -1.;

	BM b(&M);

	{
		auto e = boost::vertices(b);
		for(; e.first!=e.second; ++e.first){
			auto v = *e.first;
			std::cout << v << " " << boost::out_degree(v, b) << " " << boost::in_degree(v, b) << "\n";
		}
	}

	auto wtmap = boost::get(boost::edge_weight, b);

	for(int w=0; w<N+2; ++w){
		std::cout << "=== out edges ===\n";
		std::pair<out_edge_iterator, out_edge_iterator> e = boost::out_edges(w, b);
		for(; e.first!=e.second; ++e.first){
			boost::graph_traits<BM>::edge_descriptor const& edg = *e.first;
			std::cout << boost::source(edg, b) << ":" << boost::target(edg, b)
			          << " -- " << boost::get(wtmap, edg) << "\n";
		}
	}

	for(int w=0; w<N+2; ++w){
		std::cout << "=== in edges ===\n";
		std::pair<in_edge_iterator, in_edge_iterator> e = boost::in_edges(w, b);
		for(; e.first!=e.second; ++e.first){
			boost::graph_traits<BM>::edge_descriptor const& edg = *e.first;
			std::cout << boost::source(edg, b) << ":" << boost::target(edg, b)
			          << " -- " << boost::get(wtmap, edg) << "\n";
		}
	}

	auto n = boost::num_vertices(b);

	std::cout << "pg " << n << "\n";
	boost::print_graph(b);
}

BOOST_CONCEPT_ASSERT(( boost::BidirectionalGraphConcept<BM> ));
BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<BM> ));
BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<BM> ));
