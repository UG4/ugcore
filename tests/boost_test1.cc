
#include "lib_algebra/graph_interface/sparsematrix_boost.h"
#include "lib_algebra/cpu_algebra/sparsematrix_impl.h"
#include "lib_algebra/graph_interface/bidirectional_boost.h"
#include "lib_algebra/graph_interface/boost_util.h"

#include "common/log.cpp" // ?
#include "common/debug_id.cpp" // ?
#include "common/assert.cpp" // ?
#include "common/util/crc32.cpp" // ?
#include "common/util/ostream_buffer_splitter.cpp" // ?
#include "common/util/string_util.cpp" // ?

#include <boost/graph/graph_utility.hpp>

// bidirectional tests

static const int N=7;

void test1()
{
	typedef ug::SparseMatrix<double> T;
	typedef ug::BidirectionalMatrix<T> BM;
	ug::SparseMatrix<double> M;

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

	BM b_(&M);
	BM b;
	b = b_;

	{
		auto e = boost::vertices(b);
		for(; e.first!=e.second; ++e.first){
			auto v = *e.first;
			std::cout << v << " outdeg: " << boost::out_degree(v, b) << " indeg: " << boost::in_degree(v, b) << "\n";
		}
	}

	auto wtmap = boost::get(boost::edge_weight, b);

	for(int w=0; w<N+2; ++w){
		std::cout << "=== out edges " << w << " ===\n";
		std::pair<out_edge_iterator, out_edge_iterator> e = boost::out_edges(w, b);
		for(; e.first!=e.second; ++e.first){
			boost::graph_traits<BM>::edge_descriptor const& edg = *e.first;
			int s = boost::source(edg, b);
			assert(w == s);
			std::cout << s << "->" << boost::target(edg, b)
			          << " -- " << boost::get(wtmap, edg) << "\n";
		}
	}

	for(int w=0; w<N+2; ++w){
		std::cout << "=== in edges === " << w << "\n";
		std::pair<in_edge_iterator, in_edge_iterator> e = boost::in_edges(w, b);
		for(; e.first!=e.second; ++e.first){
			boost::graph_traits<BM>::edge_descriptor const& edg = *e.first;
			int t = boost::target(edg, b);
			std::cout << boost::source(edg, b) << "->" << t
			          << " -- " << boost::get(wtmap, edg) << "\n";
			assert(w == t);
		}
	}

	auto n = boost::num_vertices(b);

	std::cout << "pg " << n << "\n";
	boost::print_graph(b);
}

void test0()
{
	typedef ug::SparseMatrix<double> T;
	typedef ug::BidirectionalMatrix<T> BM;
	ug::SparseMatrix<double> M;

	typedef typename boost::graph_traits<BM>::in_edge_iterator in_edge_iterator;
	typedef typename boost::graph_traits<BM>::out_edge_iterator out_edge_iterator;

	M.resize_and_clear(N, N);

	M(0, 0) = 2.;
	M(0, 6) = -1.;
	M(5, 1) = 0.;
	for(unsigned i=1; i<N; ++i){
		M(i, i) = 2.;
		M(i, i-1) = 1.;
		M(i-1, i) = 1.;
	}

	BM b_(&M);
	BM b(b_);

	{
		auto e = boost::vertices(b);
		for(; e.first!=e.second; ++e.first){
			auto v = *e.first;
			std::cout << v << " " << boost::out_degree(v, b) << " " << boost::in_degree(v, b) << "\n";
		}
	}

	auto wtmap = boost::get(boost::edge_weight, b);

	{
		for(int k=0; k<N; ++k){
			std::cout << "=== out edges " << k << " === "
			          << boost::out_degree(k, b) << "\n";
			std::pair<out_edge_iterator, out_edge_iterator> e = boost::out_edges(k, b);
			// auto e = ug::util::omit_loops(e_, b);
			for(; e.first!=e.second; ++e.first){
				boost::graph_traits<BM>::edge_descriptor const& edg = *e.first;
				std::cout << boost::source(edg, b) << ":" << boost::target(edg, b)
							 << " -- " << boost::get(wtmap, edg) << "\n";
			}
		}
	}

	if(0){
		std::cout << "=== out edges ===\n";
		std::pair<out_edge_iterator, out_edge_iterator> e_ = boost::out_edges(5, b);
		auto e = ug::util::omit_loops(e_, b);
		for(; e.first!=e.second; ++e.first){
			boost::graph_traits<BM>::edge_descriptor const& edg = *e.first;
			auto wt = boost::get(wtmap, edg);
			std::cout << boost::source(edg, b) << ":" << boost::target(edg, b)
			          << " -- " << wt << "\n";
			assert(wt);
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

int main()
{
	std::cout << "== test0\n";
	test0();
	std::cout << "== test1\n";
	test1();
}

typedef ug::DenseMatrix<ug::VariableArray2<double> > T;
BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<ug::SparseMatrix<T>> ));
BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<ug::SparseMatrix<T>> ));
BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<ug::SparseMatrix<T>> ));

BOOST_CONCEPT_ASSERT(( boost::BidirectionalGraphConcept<ug::BidirectionalMatrix<ug::SparseMatrix<double>>> ));
BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<ug::BidirectionalMatrix<ug::SparseMatrix<double>>> ));
BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<ug::BidirectionalMatrix<ug::SparseMatrix<double>>> ));

