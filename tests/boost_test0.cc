// #define UG_PARALLEL

#include "lib_algebra/cpu_algebra/sparsematrix_boost.h"
#include "lib_algebra/cpu_algebra/sparsematrix_impl.h"

#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallel_matrix_boost.h"
#include "pcl/pcl_base.cpp"
#endif

#include "common/log.cpp" // ?
#include "common/debug_id.cpp" // ?
#include "common/assert.cpp" // ?
#include "common/util/crc32.cpp" // ?
#include "common/util/ostream_buffer_splitter.cpp" // ?
#include "common/util/string_util.cpp" // ?

#include "boost/graph/graph_utility.hpp"

#ifdef UG_PARALLEL
typedef ug::SparseMatrix<ug::DenseMatrix<ug::FixedArray2<double, 3, 3> > > T;
typedef ug::ParallelMatrix<T> M;
#else
typedef ug::SparseMatrix<double> M;
#endif

static const int N=7;

int main()
{
	M A;
	A.resize_and_clear(N, N);

	A(0 , N-1) = -1.;
	A(1 , N-1) = 0.;
	A(0, 0) = 2.;
	for(unsigned i=1; i<N; ++i){
		A(i, i) = 2.;
		A(i, i-1) = 1.;
		A(i-1, i) = 1.;
	}

	M B;
	B.set_as_transpose_of(A);
	boost::graph_traits<M>::adjacency_iterator i;

	auto p = boost::adjacent_vertices(size_t(1), B);
	for(; p.first!=p.second; ++p.first){
		std::cout << " " << *p.first;
	}
	std::cout << "\nedges\n";
	auto e = boost::out_edges(1, B);
	auto wtmap = boost::get(boost::edge_weight, B);
	for(; e.first!=e.second; ++e.first){
		boost::graph_traits<M>::edge_descriptor const& edg = *e.first;
		std::cout << boost::source(edg, B) << ":" << boost::target(edg, B) <<
		" -- " << boost::get(wtmap, edg) << "\n";
	}
	std::cout << "pg\n";
	boost::print_graph(A);

	{
		auto e = boost::vertices(B);
		for(; e.first!=e.second; ++e.first){
			auto v = *e.first;
			std::cout << v << " " << boost::out_degree(v, B) << "\n";
		}
	}
}

BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<ug::SparseMatrix<double>> ));
BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<ug::SparseMatrix<double>> ));
BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<ug::SparseMatrix<double>> ));
