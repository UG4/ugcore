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
 
#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_CUTHILL_MCKEE_NEW_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_CUTHILL_MCKEE_NEW_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include "lib_algebra/graph_interface/undirected.h"
#include "lib_algebra/graph_interface/undirected_boost.h"

#include <boost/graph/cuthill_mckee_ordering.hpp>

#include "IOrderingAlgorithm.h"
#include "util.cpp"

//debug
#include "common/error.h"
#include "common/log.h"

namespace ug{


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
class BoostCuthillMcKeeNewOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef Graph_t G_t;
	typedef boost::graph_traits<G_t>::vertex_descriptor Vertex_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	BoostCuthillMcKeeNewOrdering() : m_bReverse(false){}

	/// clone constructor
	BoostCuthillMcKeeNewOrdering( const BoostCuthillMcKeeOrdering<TAlgebra, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new BoostCuthillMcKeeNewOrdering<TAlgebra, O_t>(*this));
	}

	void compute(){ untested();
		UG_COND_THROW(boost::num_vertices(undir) == 0, name() << "::compute: Graph is empty!");

		typename boost::property_map<M_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, undir);

		size_t N = boost::num_vertices(undir);
		std::vector<Vertex_t> inv_perm(N);

		auto dm = boost::make_degree_map(undir);

		untested();

#if 0
		auto vc = get(boost::vertex_color, undir);
#else
		typedef boost::iterator_property_map<unsigned*,
		    boost::identity_property_map, unsigned, unsigned&> map_type;
		untested();

		std::vector<unsigned> V(N);
		map_type vc(&V[0], boost::identity_property_map());
#endif
		untested();

		if(m_bReverse){ untested();
			boost::cuthill_mckee_ordering(undir, inv_perm.rbegin(), vc, dm);
		}else{ untested();
			boost::cuthill_mckee_ordering(undir, inv_perm.begin(), vc, dm);
		}

		o.resize(boost::num_vertices(undir));

		for(unsigned i = 0; i != inv_perm.size(); ++i){
			o[index_map[inv_perm[i]]] = i;
		}

		#ifdef UG_DEBUG
		untested();
		check();
		#endif
	}

	void check(){ untested();
		UG_COND_THROW(!is_permutation(o), name() << "::check: Not a permutation!");
	}

	O_t& ordering(){ untested();
		return o;
	}

	void init(M_t* A, const V_t&){
		init(A);
	}

	void init(M_t* A){ untested();
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << "\n");
		#endif

		unsigned rows = A->num_rows();

		undir = UndirectedMatrix<M_t>(A);
	}

	void init(M_t* A, const V_t&, const O_t& inv_map){
		init(A, inv_map);
	}

	void init(M_t* A, const O_t& inv_map){ untested();
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " on induced matrix of size " << inv_map.size() << "\n");
		#endif

//		?
// 	induced_subgraph<G_t, M_t>(g, A, inv_map);
	}

	void set_reverse(bool b){
		m_bReverse = b;
	}

	virtual const char* name() const {
		if(m_bReverse){
			return "ReverseBoostCuthillMcKeeNewOrdering";
		}
		else{
			return "BoostCuthillMcKeeNewOrdering";
		}
	}

private:
	O_t o;
	UndirectedMatrix<M_t> undir;

	bool m_bReverse;
};

} //namespace

#endif //guard
