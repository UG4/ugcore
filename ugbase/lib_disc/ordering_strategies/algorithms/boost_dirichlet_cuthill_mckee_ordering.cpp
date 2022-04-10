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

#include "lib_algebra/ordering_strategies/algorithms/boost_cuthill_mckee_ordering.cpp"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"

#include <assert.h>
#include "common/error.h"


namespace ug{

template <typename TAlgebra, typename TDomain, typename O_t>
class BoostDirichletCuthillMcKeeOrdering : public BoostCuthillMcKeeOrdering<TAlgebra, O_t>//public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef BoostCuthillMcKeeOrdering<TAlgebra, O_t> baseclass;

	/// Grid function type for the solution
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;

	typedef Graph_t G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

	BoostDirichletCuthillMcKeeOrdering() : m_ssDirichletIdx(-1){}

	/// clone constructor
	BoostDirichletCuthillMcKeeOrdering( const BoostDirichletCuthillMcKeeOrdering<TAlgebra, TDomain, O_t> &parent )
			: baseclass(), m_ssDirichletIdx(parent.m_ssDirichletIdx){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new BoostDirichletCuthillMcKeeOrdering<TAlgebra, TDomain, O_t>(*this));
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

			g = G_t(n);

			for(unsigned i = 0; i < n; i++){
				for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
					if(conn.value() != 0.0 && conn.index() != i){
						if(!boost::edge(conn.index(), i, g).second){
							boost::add_edge(conn.index(), i, g);
						}
					}
				}
			}
		}
		catch(...){
			throw;
		}

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
		vd s = boost::add_vertex(g);
		for(unsigned i = 0; i < n; ++i){
			if(m_dirichlet[i]){
				boost::add_edge(s, i, g);
			}
		}

		baseclass::init(g);
		baseclass::set_start_vertex(s);

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

	O_t& ordering(){
		return baseclass::ordering();
	}

	void select_dirichlet_subset(const char* ssDirichletName){
		m_ssDirichletName = ssDirichletName;
	}

	virtual const char* name() const {return "BoostDirichletCuthillMcKeeOrdering";}

private:
	G_t g;

	int m_ssDirichletIdx;
	const char* m_ssDirichletName;

	std::vector<BOOL> m_dirichlet;
};

} //namespace


#endif //guard
