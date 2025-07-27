
#include "lib_algebra/cpu_algebra/sparsematrix_impl.h"
#include "lib_algebra/graph_interface/sparsematrix_boost.h"

#include "lib_algebra/graph_interface/undirected_boost.h"

// shared lib?
#include "common/log.cpp"
#include "common/debug_id.cpp"
#include "common/assert.cpp"
#include "common/util/crc32.cpp"
#include "common/util/ostream_buffer_splitter.cpp"
#include "common/util/string_util.cpp"

#include <boost/graph/graph_utility.hpp>

static const int N=7;

int main()
{
	typedef ug::SparseMatrix<double> T;
	typedef ug::UndirectedMatrix<T> UM;
	T M;

	M.resize_and_clear(N, N);

	M(0, 0) = 2.;
	M(0, 6) = -1.;
	for(unsigned i=1; i<N; ++i){
		M(i-1, i) = 1.;
	}
	M(5, 5) = 2.;
	M(4, 3) = 4.;

	/* 2  1  0  0  0  0 -1
	 * e  0  1
	 *    e  0  1
	 *       e  0  1
	 *          4  0  1
	 *             e  2  1
	 * e              e  0
	 */

	// lower fill1 0 1
	// lower fill0 0 6
	// lower fill1 1 2
	// lower fill1 2 3
	// lower fill1 4 5
	// lower fill0 5 6

	UM b(&M);

	boost::print_graph(b);

	{
		auto e = boost::vertices(b);
		for(; e.first!=e.second; ++e.first){
			auto v = *e.first;
			std::cout << v << " " << boost::out_degree(v, b) << " " << boost::in_degree(v, b) << "\n";
		}
	}

	{
		std::cout << "=== adj nodes ===\n";
		auto e = boost::adjacent_vertices(1, b);
		for(; e.first!=e.second; ++e.first){
			std::cout << *e.first << ":" <<
	//		" -- " << boost::get(wtmap, edg) <<
			"\n";
		}
	}

#if 0
	{
		std::cout << "=== out edges ===\n";
		auto e = boost::out_edges(1, b);
		for(; e.first!=e.second; ++e.first){
			boost::graph_traits<UM>::edge_descriptor const& edg = *e.first;
			std::cout << boost::source(edg, b) << ":" << boost::target(edg, b) <<
	//		" -- " << boost::get(wtmap, edg) << 
			"\n";
		}
	}
#endif
}
