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
//#include "common/debug_id.h"
#include "common/log.h"

#include <boost/graph/strong_components.hpp>

namespace ug{
template <typename TAlgebra, typename O_t>
class SCCOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> G_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;
	typedef typename boost::graph_traits<G_t>::in_edge_iterator inedge_iter;

	enum options{NO_COMPONENT_ORDERING, COMPONENT_ORDERING};

	SCCOrdering(){}

	/// clone constructor
	SCCOrdering( const SCCOrdering<TAlgebra, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new SCCOrdering<TAlgebra, O_t>(*this));
	}

	void preorder(vd v, O_t& o, size_t& i, G_t& g){
		o[i] = v;
		adj_iter nIt, nEnd;
		for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, g); nIt != nEnd; ++nIt){
			preorder(*nIt, o, ++i, g);
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

		std::cout << "n=" << n << std::endl;
		std::cout << "Total number of components: " << num_components << std::endl;

		for(unsigned i = 0; i < n; ++i){
			std::cout << "v=" << i << ", c=" << component[i] << std::endl;
		}

		std::cout << "-----------------" << std::endl;

		adj_iter nIt, nEnd;
		for(unsigned i = 0; i < n; ++i){
			std::cout << "v=" << i << "(" << component[i] << ")" << ", n={";
			for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(i, g); nIt != nEnd; ++nIt){
				std::cout << *nIt << "(" << component[*nIt] << ");";
			} std::cout << "}" << std::endl;
		}

		std::vector<std::vector<vd> > comp_members(num_components);
		for(unsigned i = 0; i < n; ++i){
			comp_members[component[i]].push_back(i);
		}

		for(unsigned i = 0; i < num_components; ++i){
			std::cout << "comp_" << i << ":" << std::endl;
			for(unsigned j = 0; j < comp_members[i].size(); ++j){
				std::cout << comp_members[i][j] << " ";
			} std::cout << std::endl;
		}

		//create scc meta graph
		scc_g = G_t(num_components); //use directed graph instead of G_t?

		size_t i_comp, n_comp;
		for(unsigned i = 0; i < n; ++i){
			i_comp = component[i];
			for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(i, g); nIt != nEnd; ++nIt){
				n_comp = component[*nIt];
				if(i_comp != n_comp){
					boost::add_edge(i_comp, n_comp, scc_g);
				}
			}
		}

		//determine root of scc_g
		vd scc_root = *boost::vertices(scc_g).first;
		inedge_iter ineIt, ineEnd;
		boost::tie(ineIt, ineEnd) = boost::in_edges(scc_root, scc_g);
		while(ineIt != ineEnd){
			scc_root = boost::source(*ineIt, scc_g);
			boost::tie(ineIt, ineEnd) = boost::in_edges(scc_root, scc_g);
		}

		std::cout << "scc root: " << scc_root << std::endl;

		O_t scc_topo_ordering(num_components);
		size_t idx = 0;
		preorder(scc_root, scc_topo_ordering, idx, scc_g);

		//scc_topo_ordering is now a topological ordering of the scc_g

		std::cout << "scc ordering: " << std::endl;
		for(unsigned i = 0; i < num_components; ++i){
			std::cout << scc_topo_ordering[i] << " ";
		} std::cout << std::endl;

		//fix this option for now
		static options option = NO_COMPONENT_ORDERING;

		o.resize(n);

		if(option == NO_COMPONENT_ORDERING){
			//enumerate members of components where components are iterated using the
			//topological ordering of the scc meta graph

			size_t k = 0;
			for(unsigned i = 0; i < num_components; ++i){
				for(unsigned j = 0; j < comp_members[scc_topo_ordering[i]].size(); ++j){
					o[k++] = comp_members[scc_topo_ordering[i]][j];
				}
			}
		}
		else if(option == COMPONENT_ORDERING){
			UG_THROW(name() << "::compute: option COMPONENT_ORDERING not implemented yet!");
		}
		else{
			UG_THROW(name() << "::compute: option ELSE not implemented yet!");
		}

		std::cout << "total ordering:" << std::endl;
		for(unsigned i = 0; i < n; ++i){
			std::cout << o[i] << " ";
		} std::cout << std::endl;

		//reset
		scc_g = G_t(0);
		g = G_t(0);

		#ifdef UG_DEBUG
		check();
		#endif

		std::cout << name() << ": done. " << std::endl;
	}

	void check(){
		if(!is_permutation(o)){
			UG_THROW(name() << "::check: Not a permutation!");
		}
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

	virtual const char* name() const {return "SCCOrdering";}

private:
	G_t g, scc_g;
	O_t o;
};


} //namespace


#endif //guard
