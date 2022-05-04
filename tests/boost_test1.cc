
#include "lib_algebra/cpu_algebra/sparsematrix_boost.h"
#include "lib_algebra/cpu_algebra/sparsematrix_impl.h"
#include "lib_algebra/graph_interface/bidirectional_boost.h"
#include "lib_algebra/graph_interface/boost_util.h"

#include "common/log.cpp" // ?
#include "common/debug_id.cpp" // ?
#include "common/assert.cpp" // ?
#include "common/util/crc32.cpp" // ?
#include "common/util/ostream_buffer_splitter.cpp" // ?
#include "common/util/string_util.cpp" // ?

#include "boost/graph/graph_utility.hpp"

static const int N=7;

int main()
{
	typedef ug::SparseMatrix<double> T;
	typedef ug::BidirectionalMatrix<T> BM;
	ug::SparseMatrix<double> M;

	typedef typename boost::graph_traits<BM>::in_edge_iterator in_edge_iterator;
	typedef typename boost::graph_traits<BM>::out_edge_iterator out_edge_iterator;

	M.resize_and_clear(N, N);

	M(0, 0) = 2.;
	M(0, 6) = -1.;
	for(unsigned i=1; i<N; ++i){
		M(i, i) = 2.;
		M(i, i-1) = 1.;
		M(i-1, i) = 1.;
	}

	BM b(&M);

	{
		auto e = boost::vertices(b);
		for(; e.first!=e.second; ++e.first){
			auto v = *e.first;
			std::cout << v << " " << boost::out_degree(v, b) << " " << boost::in_degree(v, b) << "\n";
		}
	}

	auto wtmap = boost::get(boost::edge_weight, b);

	{
		std::cout << "=== out edges ===\n";
		std::pair<out_edge_iterator, out_edge_iterator> e_ = boost::out_edges(1, b);
		auto e = ug::util::omit_loops(e_, b);
		for(; e.first!=e.second; ++e.first){
			boost::graph_traits<BM>::edge_descriptor const& edg = *e.first;
			std::cout << boost::source(edg, b) << ":" << boost::target(edg, b)
			          << " -- " << boost::get(wtmap, edg) << "\n";
		}
	}

	{
		std::cout << "=== in edges ===\n";
		std::pair<in_edge_iterator, in_edge_iterator> e = boost::in_edges(1, b);
		for(; e.first!=e.second; ++e.first){
			boost::graph_traits<BM>::edge_descriptor const& edg = *e.first;
			std::cout << boost::source(edg, b) << ":" << boost::target(edg, b)
			          << " -- " << boost::get(wtmap, edg) << "\n";
		}
	}
	
}

typedef ug::DenseMatrix<ug::VariableArray2<double> > T;
BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<ug::SparseMatrix<T>> ));
BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<ug::SparseMatrix<T>> ));
BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<ug::SparseMatrix<T>> ));
