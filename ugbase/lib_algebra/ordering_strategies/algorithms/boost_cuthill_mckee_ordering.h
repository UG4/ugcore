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
 
#ifndef UG_BASE_LIB_ALGEBRA_ORDERING_STRATEGIES_ALGORITHMS_BOOST_CUTHILL_MCKEE_ORDERING_H
#define UG_BASE_LIB_ALGEBRA_ORDERING_STRATEGIES_ALGORITHMS_BOOST_CUTHILL_MCKEE_ORDERING_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/cuthill_mckee_ordering.hpp>

#include "IOrderingAlgorithm.h"
#include "util.h"

//debug
#include "common/error.h"
#include "common/log.h"

namespace ug{


template <typename G_t, typename M_t>
void induced_subgraph(G_t& ind_g, M_t* A, const std::vector<size_t>& inv_map){
	size_t n = A->num_rows();
	size_t k = inv_map.size();
	ind_g = G_t(k);

	std::vector<int> ind_map(n, -1);
	for(unsigned i = 0; i < k; ++i){
		ind_map[inv_map[i]] = i;
	}

	typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
	for(unsigned i = 0; i < inv_map.size(); ++i){
		for(typename M_t::row_iterator conn = A->begin_row(inv_map[i]); conn != A->end_row(inv_map[i]); ++conn){
			if(conn.value() != 0.0 && conn.index() != i){
				int idx = ind_map[conn.index()];
				if(idx >= 0){
					boost::add_edge(i, idx, ind_g);
				}
			}
		}
	}
}



#ifndef MCKEE_GRAPH_T
#define MCKEE_GRAPH_T
/* boost graph type for Cuthill-McKee */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	boost::property<boost::vertex_color_t,
			 boost::default_color_type,
			 boost::property<boost::vertex_degree_t, int> > >
				Graph_t;
#endif


template <typename TAlgebra, typename O_t=std::vector<size_t> >
class BoostCuthillMcKeeOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef Graph_t G_t;
	typedef boost::graph_traits<G_t>::vertex_descriptor Vertex_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	BoostCuthillMcKeeOrdering() : m_bReverse(false){}

	/// clone constructor
	BoostCuthillMcKeeOrdering( const BoostCuthillMcKeeOrdering<TAlgebra, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new BoostCuthillMcKeeOrdering<TAlgebra, O_t>(*this));
	}

	void compute(){
		UG_COND_THROW(boost::num_vertices(g) == 0, name() << "::compute: Graph is empty!");

		boost::property_map<G_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, g);

		std::vector<Vertex_t> inv_perm(boost::num_vertices(g));

		if(m_bReverse){
			boost::cuthill_mckee_ordering(g, inv_perm.rbegin(), get(boost::vertex_color, g), boost::make_degree_map(g));
		}
		else{
			boost::cuthill_mckee_ordering(g, inv_perm.begin(), get(boost::vertex_color, g), boost::make_degree_map(g));
		}

		o.resize(boost::num_vertices(g));

		for(unsigned i = 0; i != inv_perm.size(); ++i){
			o[index_map[inv_perm[i]]] = i;
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
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << "\n");
		#endif

		unsigned rows = A->num_rows();

		g = G_t(rows);

		for(unsigned i = 0; i < rows; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){ //TODO: think about this!!
					if(!boost::edge(i, conn.index(), g).second){
						boost::add_edge(i, conn.index(), g);
					}
				}
			}
		}
	}

	void init(M_t* A, const V_t&, const O_t& inv_map){
		init(A, inv_map);
	}

	void init(M_t* A, const O_t& inv_map){
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " on induced matrix of size " << inv_map.size() << "\n");
		#endif

		induced_subgraph<G_t, M_t>(g, A, inv_map);
	}

	void set_reverse(bool b){
		m_bReverse = b;
	}

	virtual const char* name() const {
		if(m_bReverse){
			return "ReverseBoostCuthillMcKeeOrdering";
		}
		else{
			return "BoostCuthillMcKeeOrdering";
		}
	}

private:
	G_t g;
	O_t o;

	bool m_bReverse;
};

}

#endif
