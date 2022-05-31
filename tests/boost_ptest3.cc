
#define UG_PARALLEL

#include "lib_algebra/graph_interface/sparsematrix_boost.h"
#include "lib_algebra/cpu_algebra/sparsematrix_impl.h"

#include "lib_algebra/graph_interface/parallel_matrix.h"
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

// parallel_index_layout.h
// typedef pcl::SingleLevelLayout<pcl::OrderedInterface<size_t, std::vector> > IndexLayout;


typedef ug::SparseMatrix<double> T;
typedef ug::ParallelMatrix<T> PM;
typedef ug::BGLParallelMatrix<PM> BM;

template<class G>
class Map{
public:
	Map(G const& g):_g(g){}

	template<class T>
	std::string operator()(T& t) const{
//		using boost::detail::parallel::owner;
//		using boost::detail::parallel::local;
		return "(" + std::to_string(boost::owner(t)) + ":" + std::to_string(boost::local(t)) + ")";
	}
private:
	G const& _g;
};

template<class T, class G>
std::string get(Map<G> const& m, T t)
{
	return m(t);
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	PM M;

	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	sleep(rank);

	static const int N=7 + 2*size;
	M.resize_and_clear(N, N);

	M(0, 0) = 2.;
	M(1, N-1) = -1.;
	for(int i=1; i<N; ++i){
		M(i, i) = 2.;
		M(i, i-1) = 1.;
		M(i-1, i) = 1.;
	}

	ug::AlgebraLayouts* L = new ug::AlgebraLayouts();
	pcl::ProcessCommunicator pc = L->proc_comm();
	IndexLayout& idx_master = L->master();
	IndexLayout& idx_slave = L->slave();

	// typedef pcl::OrderedInterface IndexLayout::Interface
	// pcl_communication_structs.h

	// build some arbitrary horizontal interfaces
	if(rank){ // "slave"
		IndexLayout::Interface& s = idx_slave.interface(0, 0);
		s.push_back(2*rank+1);
		s.push_back(2*rank+2);

	}else{ // master
		for(int p=1; p<size; ++p){
			IndexLayout::Interface& m = idx_master.interface(p, 0);
			m.push_back(2*p);
			m.push_back(2*p+1);
		}
	}

	{
		std::cout << "underlying...\n";
		auto e = boost::vertices(M);
		for(; e.first!=e.second; ++e.first){
			int v = *e.first;
			std::cout << " " << v;
		}
		std::cout << "\n";
	}

	M.set_layouts(make_sp(L));
	std::cout << *L <<"\n";

	BM b(&M);
	typedef typename boost::graph_traits<BM>::vertex_descriptor vertex_descriptor;
	typedef typename boost::graph_traits<BM>::out_edge_iterator out_edge_iterator;

	if(0)
	{
		auto e = boost::vertices(b);
		for(; e.first!=e.second; ++e.first){
			vertex_descriptor v = *e.first;
			std::cout << owner(v) << ":" << local(v) << " " << boost::out_degree(v, b) << "\n";
		}

		e = boost::vertices(b);
		for(; e.first!=e.second; ++e.first){
			vertex_descriptor v = *e.first;
			std::cout << owner(v) << ":" << local(v) << " " << boost::out_degree(v, b) << "\n";

//	auto wtmap = boost::get(boost::edge_weight, b);

#if 1 // need global vertex map or so
			std::cout << "=== out edges ===\n";
			std::pair<out_edge_iterator, out_edge_iterator> e = boost::out_edges(v, b);
			for(; e.first!=e.second; ++e.first){
				boost::graph_traits<BM>::edge_descriptor const& edg = *e.first;
				auto t = boost::target(edg, b);
				std::cout << local(boost::source(edg, b)) << " -- " << owner(t)
					<< ":" << local(t) << "\n";
			}
#endif
		}
	}

	std::cout << "pg\n";
	Map<BM> m(b);
	boost::print_graph(b, m);

	MPI_Finalize();
	
}

// BOOST_CONCEPT_ASSERT(( boost::DirectedGraphConcept<BM> ));
// BOOST_CONCEPT_ASSERT(( boost::DistributedGraphConcept<BM> )); // checker not available.
BOOST_CONCEPT_ASSERT(( boost::GraphConcept<BM> ));
BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<BM> ));
BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<BM> ));

typedef ug::BGLParallelMatrix<ug::ParallelMatrix<ug::SparseMatrix<ug::DenseMatrix<ug::FixedArray2<double, 3, 3> > > > > BMDF;
BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<BMDF> ));
