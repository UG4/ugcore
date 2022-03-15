/*
 * Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_SCC_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_SCC_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>

#include "IOrderingAlgorithm.h"
#include "util.cpp"

//debug
#include "common/error.h"
#include "common/log.h"

#include <boost/graph/strong_components.hpp>

#include "lib_algebra/algebra_common/permutation_util.h"


namespace ug{
template <typename TAlgebra, typename O_t>
class SCCOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> G_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	typedef std::vector<size_t> ordering_container_type;
	typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> ordering_algo_type;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;
	typedef typename boost::graph_traits<G_t>::in_edge_iterator inedge_iter;

	SCCOrdering(){}

	/// clone constructor
	SCCOrdering( const SCCOrdering<TAlgebra, O_t> &parent )
			: baseclass(), m_spOrderingSubAlgo(parent.m_spOrderingSubAlgo){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new SCCOrdering<TAlgebra, O_t>(*this));
	}

	/// sets an ordering algorithm
	void set_ordering_subalgorithm(SmartPtr<ordering_algo_type> ordering_subalgo){
		m_spOrderingSubAlgo = ordering_subalgo;
	}

	vd get_source_vertex(std::vector<BOOL>& visited, G_t& g){
		for(unsigned i = 0; i < boost::num_vertices(g); ++i){
			if(!visited[i] && boost::in_degree(i, g) == 0){
				visited[i] = true;
				boost::clear_vertex(i, g);
				return i;
			}
		}

		return -1u;
	}

	void topological_ordering(O_t& o, G_t& g){
		size_t n = boost::num_vertices(g);
		std::vector<BOOL> visited(n, false);
		for(unsigned i = 0; i < n; ++i){
			o[i] = get_source_vertex(visited, g);
		}
	}

	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}

		std::vector<int> component(n), discover_time(n);
		std::vector<boost::default_color_type> color(n);
		std::vector<vd> root(n);
		size_t num_components = boost::strong_components(g,
			boost::make_iterator_property_map(component.begin(), boost::get(boost::vertex_index, g)),
			boost::root_map(boost::make_iterator_property_map(root.begin(), boost::get(boost::vertex_index, g)))
			    .color_map(
				boost::make_iterator_property_map(color.begin(), boost::get(boost::vertex_index, g)))
			    .discover_time_map(boost::make_iterator_property_map(
				discover_time.begin(), boost::get(boost::vertex_index, g))));


		std::vector<std::vector<vd> > comp_members(num_components);
		for(unsigned i = 0; i < n; ++i){
			comp_members[component[i]].push_back(i);
		}

		//create scc meta graph
		scc_g = G_t(num_components);

		adj_iter nIt, nEnd;
		size_t i_comp, n_comp;
		for(unsigned i = 0; i < n; ++i){
			i_comp = component[i];
			for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(i, g); nIt != nEnd; ++nIt){
				n_comp = component[*nIt];
				if(i_comp != n_comp && !boost::edge(i_comp, n_comp, scc_g).second){
					boost::add_edge(i_comp, n_comp, scc_g);
				}
			}
		}

		O_t scc_topo_ordering(num_components);
		topological_ordering(scc_topo_ordering, scc_g);

		//scc_topo_ordering is now a topological ordering of scc_g

		o.resize(n);

		if(m_spOrderingSubAlgo.invalid()){
			UG_LOG(name() << "::compute: not using ordering subalgo\n");
			//enumerate members of components where components are iterated using the
			//topological ordering of the scc meta graph

			size_t k = 0;
			for(unsigned i = 0; i < num_components; ++i){
				for(unsigned j = 0; j < comp_members[scc_topo_ordering[i]].size(); ++j){
					o[k++] = comp_members[scc_topo_ordering[i]][j];
				}
			}
		}
		else{
			UG_LOG(name() << ":compute: using " << m_spOrderingSubAlgo->name() << " as subalgo\n");

			size_t k = 0;
			for(unsigned i = 0; i < num_components; ++i){
				size_t c_size = comp_members[scc_topo_ordering[i]].size();

				if(c_size == 1){
					o[comp_members[scc_topo_ordering[i]][0]] = k++;
				}
				else{
					m_spOrderingSubAlgo->init(m, comp_members[scc_topo_ordering[i]]);
					m_spOrderingSubAlgo->compute();
					std::vector<size_t>& sub_o = m_spOrderingSubAlgo->ordering();
					std::vector<size_t> inv_sub_o;
					GetInversePermutation(sub_o, inv_sub_o);

					for(unsigned j = 0; j < c_size; ++j){
						o[comp_members[scc_topo_ordering[i]][inv_sub_o[j]]] = k++;
					}
				}
			}

			UG_COND_THROW(k != n, "k!=n, k=" << k << ", n=" << n);
		}

		//reset
		scc_g = G_t(0);
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
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << "\n");
		#endif

		m = A;

		unsigned rows = A->num_rows();

		g = G_t(rows);

		for(unsigned i = 0; i < rows; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){ //TODO: think about this!!
					//double w;
					//w = abs(conn.value()); //TODO: think about this
					//boost::add_edge(i, conn.index(), w, g);
					boost::add_edge(i, conn.index(), g);
				}
			}
		}
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	virtual const char* name() const {return "SCCOrdering";}

private:
	G_t g, scc_g;
	O_t o;

	M_t* m;

	SmartPtr<ordering_algo_type> m_spOrderingSubAlgo;
};


} //namespace


#endif //guard
