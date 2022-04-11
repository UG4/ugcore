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
 
#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_CUTHILL_MCKEE_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_CUTHILL_MCKEE_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/cuthill_mckee_ordering.hpp>

#include "IOrderingPreprocessor.h"
#include "IOrderingAlgorithm.h"
#include "util.cpp"

//debug
#include "common/error.h"
#include "common/log.h"

namespace ug{

#ifndef GRAPH_T_FOR_CUTHILL_MCKEE
#define GRAPH_T_FOR_CUTHILL_MCKEE
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
	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;
	typedef IOrderingPreprocessor<TAlgebra, G_t> Preprocessor_t;

	BoostCuthillMcKeeOrdering() : m_bReverse(false), m_bUseStartVertex(false),
			m_vdStartVertex(0){}

	/// clone constructor
	BoostCuthillMcKeeOrdering( const BoostCuthillMcKeeOrdering<TAlgebra, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse), m_bUseStartVertex(parent.m_bUseStartVertex),
			m_vdStartVertex(parent.m_vdStartVertex){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new BoostCuthillMcKeeOrdering<TAlgebra, O_t>(*this));
	}

	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}

		if(m_spPreprocessor.valid()){
			UG_LOG("call preprocessor\n");
			m_spPreprocessor->compute();
		}

		boost::property_map<G_t, boost::vertex_degree_t>::type deg = get(boost::vertex_degree, g);
		boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
			deg[*vIt] = boost::degree(*vIt, g);
		}

		boost::property_map<G_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, g);

		typedef boost::graph_traits<G_t>::vertex_descriptor Vertex_t;
		std::vector<Vertex_t> inv_perm(boost::num_vertices(g));

		if(m_bUseStartVertex){
			if(m_bReverse){
				boost::cuthill_mckee_ordering(g, m_vdStartVertex, inv_perm.rbegin(), get(boost::vertex_color, g), boost::make_degree_map(g));
			}
			else{
				boost::cuthill_mckee_ordering(g, m_vdStartVertex, inv_perm.begin(), get(boost::vertex_color, g), boost::make_degree_map(g));
			}
		}
		else{
			if(m_bReverse){
				boost::cuthill_mckee_ordering(g, inv_perm.rbegin(), get(boost::vertex_color, g), boost::make_degree_map(g));
			}
			else{
				boost::cuthill_mckee_ordering(g, inv_perm.begin(), get(boost::vertex_color, g), boost::make_degree_map(g));
			}
		}

		size_t offset = 0;

		if(m_bUseStartVertex){
			offset = 1;
		}

		if(m_bReverse){
			//remove start vertex (at end of o)
			o.resize(boost::num_vertices(g)-offset);

			for(unsigned i = 0; i != inv_perm.size()-offset; ++i){
				o[index_map[inv_perm[i]]] = i;
			}
		}
		else{
			//remove start vertex (at begin of o)
			o.resize(boost::num_vertices(g)-offset);

			for(unsigned i = 0+offset; i != inv_perm.size(); ++i){
				o[index_map[inv_perm[i]]] = i-offset;
			}
		}
		g = G_t(0);

		#ifdef UG_DEBUG
		check();
		#endif
	}

	void init(M_t* A, const V_t& V){
		init(A);

		if(m_spPreprocessor.valid()){
			m_spPreprocessor->init(V, g);
		}
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
				if(conn.value() != 0.0 && conn.index() != i){
					if(!boost::edge(conn.index(), i, g).second){
						boost::add_edge(conn.index(), i, g); //convection is reverse to disc. direction
					}
				}
			}
		}
	}

	void init(M_t* A, const V_t&, const O_t& inv_map){
		init(A, inv_map);
	}

	//TODO: avoid copy
	void init(M_t* A, const O_t& inv_map){
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " on induced matrix of size " << inv_map.size() << "\n");
		#endif

		induced_subgraph<G_t, M_t>(g, A, inv_map);
	}

	//called with an assembled graph by lib_disc version
	//TODO: avoid copy
	void init(G_t &g_in){
		g = g_in;
	}

	void check(){
		UG_COND_THROW(!is_permutation(o), name() << "::check: Not a permutation!");
	}

	O_t& ordering(){
		return o;
	}

	virtual const char* name() const {
		if(m_bReverse){
			return "ReverseBoostCuthillMcKeeOrdering";
		}
		else{
			return "BoostCuthillMcKeeOrdering";
		}
	}

	void set_reverse(bool b){
		m_bReverse = b;
	}

	void set_start_vertex(vd s){
		m_vdStartVertex = s;
		m_bUseStartVertex = true;
	}

	void set_preprocessor(SmartPtr<Preprocessor_t> spPreprocessor){
		m_spPreprocessor = spPreprocessor;
	}

private:
	G_t g;
	O_t o;

	bool m_bReverse;

	bool m_bUseStartVertex;
	vd m_vdStartVertex;

	SmartPtr<Preprocessor_t> m_spPreprocessor;
};

} //namespace

#endif //guard
