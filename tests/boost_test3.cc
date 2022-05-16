
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
#include <boost/graph/cuthill_mckee_ordering.hpp>

static const int N=4;
typedef ug::SparseMatrix<double> T;
typedef ug::UndirectedMatrix<T> UM;

void test0()
{
	T M;

	M.resize_and_clear(N, N);

	M(0, 0) = 1.;
	M(0, 1) = 1.;
	M(0, 2) = 1.;
	M(0, 3) = 1.;

	/* 1  1  1  1
	 * 0  0  0  0
	 * 0  0  0  0
	 * 0  0  0  0
	 */


	UM b(&M);

	boost::print_graph(b);
//	boost::print_graph(b._extra_fill);
	auto dm = boost::make_degree_map(b);
	auto im = get(boost::vertex_index, b);

	{
		auto e = boost::vertices(b);
		for(; e.first!=e.second; ++e.first){
			auto v = *e.first;
			auto d = boost::get(dm, v);
			auto ind = boost::get(im, v);
			std::cout << "vertex " << v << " deg: " << d << "\n";

			assert(v == ind);
			assert(d==boost::out_degree(v, b));
			assert(d==boost::in_degree(v, b));
			assert(d==boost::degree(v, b));
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

#if 1
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

void test1()
{
	typedef ug::SparseMatrix<double> T;
	typedef ug::UndirectedMatrix<T> UM;
	T M;

	M.resize_and_clear(N, N);

	M(0, 0) = 1.;
	M(0, 2) = 2.;
	M(1, 0) = 3.;
	M(2, 3) = 4.;
	M(3, 2) = 4.;

	/* 1  0  2  0
	 * 3  0  0  0
	 * 0  0  0  4
	 * 0  0  4  0
	 */

	// lower fill1 2 0
	// upper fill  0 1
	//
	// 0 <--> 0 2 1
	// 1 <--> 0
	// 2 <--> 3 0
	// 3 <--> 2

	UM b(&M);

	boost::print_graph(b);
//	boost::print_graph(b._extra_fill);

	{
		auto e = boost::vertices(b);
		for(; e.first!=e.second; ++e.first){
			auto v = *e.first;
			std::cout << v << " " << boost::out_degree(v, b) << " " << boost::in_degree(v, b) << "\n";
		}
	}

	{
		std::cout << "=== adj nodes 2 ===\n";
		auto e = boost::adjacent_vertices(2, b);
		for(; e.first!=e.second; ++e.first){
			std::cout << " " << *e.first <<
	//		" -- " << boost::get(wtmap, edg) <<
			"\n";
		}
	}

	{
		std::cout << "=== out edges 2 ===\n";
		auto e = boost::out_edges(2, b);
		for(; e.first!=e.second; ++e.first){
			boost::graph_traits<UM>::edge_descriptor const& edg = *e.first;
			std::cout << boost::source(edg, b) << ":" << boost::target(edg, b) <<
	//		" -- " << boost::get(wtmap, edg) << 
			"\n";
		}
	}
}

void test2()
{
	typedef ug::SparseMatrix<double> T;
	typedef ug::UndirectedMatrix<T> UM;
	T M;

	M.resize_and_clear(N, N);

	M(0, 0) = 1.;
	M(0, 2) = 2.;
	M(1, 0) = 3.;
	M(2, 3) = 4.;
	M(3, 2) = 4.;

	/* 1  0  2  0
	 * 3  0  0  0
	 * 0  0  0  4
	 * 0  0  4  0
	 */

	// lower fill1 2 0
	// upper fill  0 1
	//
	// 0 <--> 0 2 1
	// 1 <--> 0
	// 2 <--> 3 0
	// 3 <--> 2

	UM b(&M);

	size_t N = boost::num_vertices(b);
	assert(N);
	std::vector<int> inv_perm(N);

	auto dm = boost::make_degree_map(b);

	typedef boost::iterator_property_map<unsigned*,
		 boost::identity_property_map, unsigned, unsigned&> map_type;

	std::vector<unsigned> V(N);
	map_type vc(&V[0], boost::identity_property_map());

	boost::cuthill_mckee_ordering(b, inv_perm.rbegin(), vc, dm);

	for(unsigned i = 0; i != inv_perm.size(); ++i){
		std::cout << " " << inv_perm[i];
	}
	std::cout << "\n";

} // test2: cmk

int main()
{
	std::cout << "== test0 ==\n";
	test0();
	std::cout << "== test1 ==\n";
	test1();
	std::cout << "== test2 ==\n";
	test2();
}

BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<UM> ));
BOOST_CONCEPT_ASSERT(( boost::GraphConcept<UM> ));
BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<UM> ));
BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<UM> ));
