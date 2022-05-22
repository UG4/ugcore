/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Author: Lukas Larisch
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_SPARSEMATRIX_GRAPH_TEST__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_SPARSEMATRIX_GRAPH_TEST__

#include "lib_algebra/graph_interface/sparsematrix_boost.h"
#include "lib_algebra/graph_interface/parallel_matrix.h"
#include "lib_algebra/graph_interface/parallel_matrix_boost.h"
#include "lib_algebra/graph_interface/bidirectional.h"
#include "lib_algebra/graph_interface/bidirectional_boost.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>

#include <vector>
#include <utility> //for pair
#include <climits> //INT_MAX

#include "IOrderingAlgorithm.h"
#include "util.cpp"

//debug
#include "common/error.h"
#include "common/log.h"


namespace ug{

namespace {

template<class G>
class Map{
public:
	Map(G const& g):_g(g){}

	template<class T>
	std::string operator()(T& t) const{
//		using boost::detail::parallel::owner;
//		using boost::detail::parallel::local;
		return "(" + std::to_string(owner(t)) + ":" + std::to_string(local(t)) + ")";
	}
private:
	G const& _g;
};

template<class T, class G>
std::string get(Map<G> const& m, T t)
{
	return m(t);
}

}

//for cycle-free matrices only
template <typename TAlgebra, typename O_t>
class SparseMatrixGraphTest : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> G_t;
	typedef boost::graph_traits<G_t>::vertex_descriptor vd_t;
	typedef boost::graph_traits<G_t>::vertex_iterator vIt_t;
	typedef boost::graph_traits<G_t>::adjacency_iterator nIt_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	SparseMatrixGraphTest(){}

	/// clone constructor
	SparseMatrixGraphTest( const SparseMatrixGraphTest<TAlgebra, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new SparseMatrixGraphTest<TAlgebra, O_t>(*this));
	}

	//TODO: this is a very inefficient implementation, rewrite if necessary
	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){

			return;
			UG_THROW(name() << "::compute: Graph is empty!");
		}

		o.resize(n);

		for(size_t i = 0; i < n; ++i){
			o[i] = i;
		}

		g = G_t(0);

		#ifdef UG_DEBUG
		check();
		#endif
	}

	void check(){
		UG_COND_THROW(!is_permutation(o), name() << "::check: Not a permutation!");
	}

	O_t& ordering(){
		return o;
	}

	void init(M_t* A, const V_t&){
		init(A);
	}

	void init(M_t* A){
		if(!A){
			std::cerr<<"no A?!\n";
			return;
		}else{
		}

		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << "\n");
		#endif

		m_A = A;

		unsigned rows = A->num_rows();

		g = G_t(rows);

		for(unsigned i = 0; i < rows; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){ //TODO: think about this!!
					boost::add_edge(i, conn.index(), g);
				}
			}
		}

		typedef ug::SparseMatrix<double> T;

		boost::graph_traits<T>::adjacency_iterator i;

		if(boost::num_vertices(*A)>2){
			auto p = boost::adjacent_vertices(1, *A);
			for(; p.first!=p.second; ++p.first){
				std::cout << " " << *p.first;
			}
			std::cout << "\n";
			auto e = boost::out_edges(1, *A);
			for(; e.first!=e.second; ++e.first){
				auto edg = *e.first;
				std::cout << boost::source(edg, *A) << ":" << boost::target(edg, *A) << "\n";
			}
		}

		std::cout << "pg " << boost::num_vertices(*A) << "\n";

		// need to print to cerr, because cout from rank!=0 gets lost
		boost::print_graph(*A, std::cerr);

		std::cout << "bidirectional" << std::endl;

		BidirectionalMatrix<M_t> bim(A);

		for(int i = 0; i < boost::num_vertices(bim); ++i){
			auto inedge = boost::in_edges(i, bim);
			for(; inedge.first != inedge.second; ++inedge.first){
				if(boost::source(*inedge.first, bim) != boost::target(*inedge.first, bim)){
				}
			}
		}
#ifdef UG_PARALLEL

		std::cout << "parallel pg \n";
		typedef BGLParallelMatrix<M_t> BPM;
		BPM bpm(A);
		Map<BPM> map(bpm);
		// need to print to cerr, because cout from rank!=0 gets lost
		boost::print_graph(bpm, map, std::cerr);

#endif

#if 0
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		sleep(rank);

		std::cerr << " process " << rank << " num vert " << boost::num_vertices(*A) << "\n";
		boost::print_graph(*A, std::cerr);

		// ../../ugcore/ugbase/pcl/pcl_communication_structs.h
		ConstSmartPtr<ug::AlgebraLayouts> L = A->layouts();

		IndexLayout const& idx_master = L->master();
		IndexLayout const& idx_slave = L->slave();

		// typedef pcl::SingleLevelLayout<pcl::OrderedInterface<size_t, std::vector> >
		// 		IndexLayout;
		// 		../../ugcore/ugbase/pcl/pcl_communication_structs.h
		// typedef pcl::OrderedInterface<size_t, std::vector> Interface
		//


		for(auto k = idx_master.begin(); k!=idx_master.end(); ++k){
			std::cerr << "rank " << rank << " master " << idx_master.proc_id(k) << "\n";
			IndexLayout::Interface const& I = idx_master.interface(k);
			for(auto i : I){
				std::cerr << " [ elem " << i.elem;
				std::cerr << " id " << i.localID << "]";
			}
			std::cerr << "\n";
		}
		for(auto k = idx_slave.begin(); k!=idx_slave.end(); ++k){
			std::cerr << "rank " << rank << " slave " << idx_slave.proc_id(k) << "\n";
			IndexLayout::Interface const& I = idx_slave.interface(k);
			for(auto i : I){
				std::cerr << " [ elem " << i.elem;
				std::cerr << " id " << i.localID << "]";
			}
			std::cerr << "\n";
		}
		assert(L);
		std::cerr << "rank: " << rank << "\n"
			       << "overlap: " << L->overlap_enabled() << "\n"
		          << "Layouts\n " << *L << "\n";
#endif
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	virtual const char* name() const {return "SparseMatrixGraphTest";}

private:
	G_t g;
	O_t o;
	M_t* m_A;
};

} //namespace

#endif //guard
