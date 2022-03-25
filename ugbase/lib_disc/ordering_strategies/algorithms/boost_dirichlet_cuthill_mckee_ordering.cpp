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

#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_BOOST_DIRICHLET_CUTHILL_MCKEE_ORDERING__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_BOOST_DIRICHLET_CUTHILL_MCKEE_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/cuthill_mckee_ordering.hpp>

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

#ifndef GRAPH_T_FOR_CUTHILL_MCKEE
#define GRAPH_T_FOR_CUTHILL_MCKEE
/* boost graph type for Cuthill-McKee */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	boost::property<boost::vertex_color_t,
			 boost::default_color_type,
			 boost::property<boost::vertex_degree_t, int> > >
				Graph_t;
#endif

template <typename TAlgebra, typename TDomain, typename O_t>
class BoostDirichletCuthillMcKeeOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	/// Grid function type for the solution
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;

	typedef Graph_t G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;
	typedef typename boost::graph_traits<G_t>::out_edge_iterator oute_iter;

	BoostDirichletCuthillMcKeeOrdering() : m_bReverse(false), m_ssDirichletIdx(-1){}

	/// clone constructor
	BoostDirichletCuthillMcKeeOrdering( const BoostDirichletCuthillMcKeeOrdering<TAlgebra, TDomain, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse), m_ssDirichletIdx(parent.m_ssDirichletIdx){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new BoostDirichletCuthillMcKeeOrdering<TAlgebra, TDomain, O_t>(*this));
	}

	//overload
	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}

		boost::property_map<G_t, boost::vertex_degree_t>::type deg = get(boost::vertex_degree, g);
		boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
			deg[*vIt] = boost::degree(*vIt, g);
		}

		boost::property_map<G_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, g);

		typedef boost::graph_traits<G_t>::vertex_descriptor Vertex_t;
		std::vector<Vertex_t> inv_perm(boost::num_vertices(g));

		if(m_bReverse){
			boost::cuthill_mckee_ordering(g, s, inv_perm.rbegin(), get(boost::vertex_color, g), boost::make_degree_map(g));
		}
		else{
			boost::cuthill_mckee_ordering(g, s, inv_perm.begin(), get(boost::vertex_color, g), boost::make_degree_map(g));
		}

		//skip s
		o.resize(boost::num_vertices(g)-1);
		for(unsigned i = 1; i != inv_perm.size(); ++i){
			o[index_map[inv_perm[i]]] = i-1;
		}

		g = G_t(0);

		#ifdef UG_DEBUG
		check();
		#endif


	}

	void init(M_t* A, const V_t& V){
		const GridFunc_t* pGridF;
		size_t numDirichlet = 0;
		unsigned n;

		try{
			if((pGridF = dynamic_cast<const GridFunc_t*>(&V)) == 0){
				UG_THROW(name() << "::init: No DoFDistribution specified.");
			}

			m_ssDirichletIdx = pGridF->domain()->subset_handler()->get_subset_index(m_ssDirichletName);

			if(m_ssDirichletIdx < 0)
				UG_THROW(name() << "::init: Invalid subset for sources selected! Call 'select_sources(const char*)'.");

			n = A->num_rows();


/* REMOVE THIS */
	typedef typename std::pair<MathVector<TDomain::dim>, size_t> Position_t;
	std::vector<Position_t> vPositions;
	SmartPtr<DoFDistribution> dd = ((GridFunc_t*) pGridF)->dof_distribution();
	ExtractPositions(pGridF->domain(), dd, vPositions);
/* REMOVE THIS*/

			g = G_t(n);

			for(unsigned i = 0; i < n; i++){
/* REMOVE THIS */
						UG_LOG(i << ": " << vPositions[i].first << "\n");
/* REMOVE THIS*/

				for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
					if(conn.value() != 0.0 && conn.index() != i){ //TODO: think about this!!
						boost::add_edge(i, conn.index(), g);
/* REMOVE THIS */
						UG_LOG(i << " -> " << conn.index() << ", " << vPositions[i].first << " -> " << vPositions[conn.index()].first << "\n");
/* REMOVE THIS*/
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

		if(numDirichlet == 0){
			UG_THROW(name() << "::init: #dirichlet_nodes = 0!");
		}

		//create a vertex and connect it to dirichlet nodes
		s = boost::add_vertex(g);
		for(unsigned i = 0; i < n; ++i){
			if(m_dirichlet[i]){
				boost::add_edge(s, i, g);
			}
		}

		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " (subset " << m_ssDirichletIdx << ", " << m_ssDirichletName
				<< ", n=" << boost::num_vertices(g) << ", m=2*" << boost::num_edges(g)/2
				<< ", d=" <<  numDirichlet << ")\n");
		#endif
	}

	void init(M_t* A){
		UG_THROW(name() << "::init: Cannot initialize smoother without a geometry. Specify the 2nd argument for init!");
	}

	void init(M_t* A, const V_t&, const O_t& inv_map){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void init(M_t* A, const O_t& inv_map){
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

	void set_reverse(bool b){
		m_bReverse = b;
	}

	void select_dirichlet_subset(const char* ssDirichletName){
		m_ssDirichletName = ssDirichletName;
	}

	virtual const char* name() const {return "BoostDirichletCuthillMcKeeOrdering";}

private:
	G_t g;
	O_t o;

	bool m_bReverse;
	int m_ssDirichletIdx;
	const char* m_ssDirichletName;

	std::vector<BOOL> m_dirichlet;
	vd s;
};

} //namespace


#endif //guard
