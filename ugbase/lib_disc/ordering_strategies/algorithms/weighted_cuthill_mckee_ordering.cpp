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

#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_WEIGHTED_CUTHILL_MCKEE_ORDERING__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_WEIGHTED_CUTHILL_MCKEE_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <set>
#include <algorithm> //reverse
#include <utility> //pair

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"

#include <assert.h>
#include "common/error.h"


namespace ug{

template <typename S_t>
void print(S_t &s){
	for(auto sIt = s.begin(); sIt != s.end(); ++sIt){
		std::cout << sIt->source << " -> " << sIt->target << " ( " << sIt->w << ")" << std::endl;
	} std::cout << std::endl;
}

#ifndef PRINT_GRAPH
#define PRINT_GRAPH
template <typename G_t>
void print_graph(G_t& g){
	typedef typename boost::graph_traits<G_t>::edge_descriptor Edge;

	typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
	for(boost::tie(eIt, eEnd) = boost::edges(g); eIt != eEnd; ++eIt){
		std::pair<Edge, bool> e = boost::edge(boost::source(*eIt, g), boost::target(*eIt, g), g);
		double w = boost::get(boost::edge_weight_t(), g, e.first);
		std::cout << boost::source(*eIt, g) << " -> " << boost::target(*eIt, g) << " ( " << w << " )" << std::endl;
	}
}
#endif

template <typename T>
void print(std::set<T> &s){
	for(auto sIt = s.begin(); sIt != s.end(); ++sIt){
		std::cout << *sIt << " ";
	} std::cout << std::endl;
}

template <typename T>
void print(std::vector<T> &s){
	for(auto sIt = s.begin(); sIt != s.end(); ++sIt){
		std::cout << *sIt << " ";
	} std::cout << std::endl;
}

#define DUMP

/* TODO

	use iterators

*/


template <typename TAlgebra, typename TDomain, typename O_t>
class WeightedCuthillMcKeeOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	/// Grid function type for the solution
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;

	typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeWeightProperty> G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;
	typedef typename boost::graph_traits<G_t>::out_edge_iterator oute_iter;

	WeightedCuthillMcKeeOrdering() : m_bReverse(false), m_ssDirichletIdx(-1){}

	/// clone constructor
	WeightedCuthillMcKeeOrdering( const WeightedCuthillMcKeeOrdering<TAlgebra, TDomain, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse), m_ssDirichletIdx(parent.m_ssDirichletIdx){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new WeightedCuthillMcKeeOrdering<TAlgebra, TDomain, O_t>(*this));
	}

	double compute_inflow(vd v, std::vector<BOOL>& visited, bool ignore_visited=true){
		double w = .0f;
		typename boost::graph_traits<G_t>::in_edge_iterator in_nIt, in_nEnd;
		for(boost::tie(in_nIt, in_nEnd) = boost::in_edges(v, g); in_nIt != in_nEnd; ++in_nIt){
			if(ignore_visited && visited[v]){ continue; }
#ifdef UG_CPU_1
			w += abs(boost::get(boost::edge_weight_t(), g, *in_nIt)); //TODO: think about this!
#endif
		}
		return w;
	}

	//TODO: do not recompute
	vd min_inflow_vertex(std::set<vd>& front, std::vector<BOOL>& visited){
		double min_w = -1u;
		vd min_vertex;

		assert(front.size());

		typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
		for(auto sIt = front.begin(); sIt != front.end(); ++sIt){
			double w = compute_inflow(*sIt, visited);
#ifdef DUMP
			std::cout << *sIt << ": " << w << std::endl;
#endif
			if(w < min_w){
				min_w = w;
				min_vertex = *sIt;
			}
		}
		return min_vertex;
	}


	//overload
	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}

		o.resize(n);

		UG_LOG(name() << "::compute: n = " << n << ", m = " << boost::num_edges(g) << std::endl);
#ifdef DUMP
		std::cout << "graph: " << std::endl; print_graph(g); std::cout << "end graph" << std::endl;
#endif


		for(unsigned i = 0; i < n; ++i){
			o[i] = i;
		}

#ifdef REWORK_DONE
		vd cur;
		std::vector<BOOL> visited(n, false);

		for(unsigned k = 0; k < n; ++k){
			/*
				the algorithm can be tuned here: start vertex and next vertex
			*/

			visited[cur] = true;

			//insert downstream neighbors to front
			oute_iter out_nIt, out_nEnd;
			for(boost::tie(out_nIt, out_nEnd) = boost::out_edges(cur, g); out_nIt != out_nEnd; ++out_nIt){
				if(!visited[boost::target(*out_nIt, g)]){
					//front.insert(boost::target(*out_nIt, g));
				}
			}

			o[k] = cur;

#ifdef DUMP
			std::cout << "ordering: ";
			for(unsigned i = 0; i < k; ++i){
				std::cout << o[i] << " ";
			}std::cout << std::endl;
#endif
		}

		if(m_bReverse){
			std::reverse(o.begin(), o.end());
		}

#endif
		g = G_t(0);

		std::cout << "[CuthillMcKeeOrderingWeighted::compute] done. " << std::endl;
	}

	void init(M_t* A, const V_t& V){
		if(m_ssDirichletIdx < 0)
			UG_THROW(name() << "::init: No dirichlet subset selected! Call 'select_dirichlet_subset(int)'.");

		const GridFunc_t* pGridF;
		size_t numDirichlet = 0;
		std::string ssDirichletName;
		unsigned n;

		try{
			if((pGridF = dynamic_cast<const GridFunc_t*>(&V)) == 0){
				UG_THROW(name() << "::init: No DoFDistribution specified.");
			}

			//throws an exception if m_ssIdx is invalid
			ssDirichletName = pGridF->domain()->subset_handler()->get_subset_name(m_ssDirichletIdx);

			UG_LOG("ssDirichletName = " << ssDirichletName << "\n");

			n = A->num_rows();

			g = G_t(n);

			for(unsigned i = 0; i < n; i++){
				for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
					if(conn.value() != 0.0 && conn.index() != i){ //TODO: think about this!!
						//double w;
						//w = abs(conn.value()); //TODO: think about this
						//boost::add_edge(i, conn.index(), w, g);
					}
				}
			}
		}
		catch(...){
			throw;
		}

		o.resize(n);
		m_dirichlet = std::vector<BOOL>(n, false);

		//select dirichlet vertices according to m_ssDirichletIdx
		typedef typename GridFunc_t::template traits<ug::Vertex>::const_iterator ugVertIt_t;

		ugVertIt_t ugVertIt = pGridF->template begin<ug::Vertex>();
		ugVertIt_t ugVertEnd = pGridF->template end<ug::Vertex>();
		size_t k = 0;
		for(; ugVertIt != ugVertEnd; ++ugVertIt, ++k){
			ug::Vertex* v = *ugVertIt;

			int si = pGridF->domain()->subset_handler()->get_subset_index(v);

			if(si == m_ssDirichletIdx){
				m_dirichlet[k] = true;
				++numDirichlet;
			}
		}

		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " (subset " << m_ssDirichletIdx << ", " << ssDirichletName
				<< ", n=" << boost::num_vertices(g) << ", m=2*" << boost::num_edges(g)/2
				<< ", d=" <<  numDirichlet << ")\n");
		#endif
	}

	void init(M_t*){
		UG_THROW(name() << "::init: Cannot initialize smoother without a geometry. Specify the 2nd argument for init!");
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void check(){
		if(!is_permutation(o)){
			print(o);
			UG_THROW(name() << "::check: Not a permutation!");
		}
	}

	O_t& ordering(){
		return o;
	}

	SmartPtr<LuaOrdering> get_lua_ordering(){
		SmartPtr<LuaOrdering> lua_ord = SmartPtr<LuaOrdering>(new LuaOrdering());
		lua_ord->ordering = o;
		return lua_ord;
	}

	void set_reverse(bool b){
		m_bReverse = b;
	}

	void select_dirichlet_subset(int ssDirichletIdx){
		m_ssDirichletIdx = ssDirichletIdx;
	}

	virtual const char* name() const {return "WeightedCuthillMcKeeOrdering";}

private:
	G_t g;
	O_t o;

	bool m_bReverse;
	int m_ssDirichletIdx;

	std::vector<BOOL> m_dirichlet;
};

} //namespace


#endif //guard
