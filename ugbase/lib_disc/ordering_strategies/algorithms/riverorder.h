/*
 * Copyright (c) 2011-2022:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS__RIVERORDER__
#define __H__UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS__RIVERORDER__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
//#include <boost/graph/properties.hpp>

#include <vector>
#include <utility> //for pair
#include <climits> //INT_MAX

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/util.h"

#include "common/error.h"

namespace ug{



template <typename TAlgebra, typename TDomain, typename O_t>
class RiverOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> G_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	/// Grid function type for the solution
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;

	RiverOrdering() : m_ssIdx(-1){}

	/// clone constructor
	RiverOrdering( const LexOrdering<TAlgebra, TDomain, O_t> &parent )
			: baseclass(), m_ssIdx(parent.m_ssIdx){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new RiverOrdering<TAlgebra, TDomain, O_t>(*this));
	}


	vd get_source_vertex(std::vector<BOOL>& visited, G_t& g){
		for(unsigned i = 0; i < boost::num_vertices(g); ++i){
			if(!visited[i] && m_sources[i]){
				if(boost::in_degree(i, g) == 1){
					visited[i] = true;
					adj_iter nIt, nEnd;
					boost::tie(nIt, nEnd) = boost::adjacent_vertices(i, g);
					m_sources[*nIt] = true;
					boost::clear_vertex(i, g);
					return i;
				}
				else if(boost::in_degree(i, g) == 0){
					visited[i] = true;
					return i;
				}
				else{}
			}
		}

		return -1u;
	}

	void topological_ordering(O_t& o, G_t& g){
		size_t n = boost::num_vertices(g);
		std::vector<BOOL> visited(n, false);
		for(unsigned i = 0; i < n; ++i){
			auto v = get_source_vertex(visited, g);
			o[v] = i;
		}
	}

	void compute(){
		topological_ordering(o, g);

		g = G_t(0);

		#ifdef UG_DEBUG
		check();
		#endif
	}

	void check(){
		if(!is_permutation(o)){
			UG_THROW(name() << "::check: Not a permutation!");
		}
	}

	O_t& ordering(){
		return o;
	}

	void init(M_t* A, const V_t& V){
		const GridFunc_t* pGridF;
		size_t numSources = 0;
		std::string ssName;
		size_t n = 0;

		//graph construction
		try{
			if((pGridF = dynamic_cast<const GridFunc_t*>(&V)) == 0){
				UG_THROW(name() << "::init: No DoFDistribution specified.");
			}

			SmartPtr<DoFDistribution> dd = ((GridFunc_t*) pGridF)->dof_distribution();

			n = dd->num_indices();

			if(n != A->num_rows ()){
				UG_THROW(name() << "::init: #indices != #rows");
			}

			m_ssIdx = pGridF->domain()->subset_handler()->get_subset_index(m_ssName);

			if(m_ssIdx < 0)
				UG_THROW(name() << "::init: Invalid subset for sources selected! Call 'select_sources(const char*)'.");

			std::vector<std::vector<size_t> > vvConnection(n);
			try{
				dd->get_connections<ug::Edge>(vvConnection);
			}
			UG_CATCH_THROW(name() << "::init: No adjacency graph available!");

			g = G_t(n);

			for(unsigned i = 0; i < n; ++i){
				for(unsigned j = 0; j < vvConnection[i].size(); ++j){
					if(i != vvConnection[i][j]){ //no self loops!
						boost::add_edge(i, vvConnection[i][j], g);
					}
				}
			}
		}
		catch(...){
			throw;
		}

		o.resize(n);
		m_sources = std::vector<BOOL>(n, false);

		//select source vertices according to m_ssIdx
		typedef typename GridFunc_t::template traits<ug::Vertex>::const_iterator ugVertIt_t;

		ugVertIt_t ugVertIt = pGridF->template begin<ug::Vertex>();
		ugVertIt_t ugVertEnd = pGridF->template end<ug::Vertex>();
		size_t k = 0;
		for(; ugVertIt != ugVertEnd; ++ugVertIt, ++k){
			ug::Vertex* v = *ugVertIt;

			int si = pGridF->domain()->subset_handler()->get_subset_index(v);

			if(si == m_ssIdx){
				m_sources[k] = true;
				++numSources;
			}
		}

		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " (subset " << m_ssIdx << ", " << m_ssName
				<< ", n=" << boost::num_vertices(g) << ", m=2*" << boost::num_edges(g)/2
				<< ", s=" <<  numSources << ")\n");
		#endif
	}

	void init(M_t*){
		UG_THROW(name() << "::init: Cannot initialize smoother without a geometry. Specify the 2nd argument for init!");
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	virtual const char* name() const {return "RiverOrdering";}

	void select_sources(const char* ssName){
		//m_ssIdx = ssIdx;
		m_ssName = ssName;
	}

private:
	G_t g;
	O_t o;

	int m_ssIdx;
	const char* m_ssName;
	std::vector<BOOL> m_sources;
};

} // end namespace ug

#endif
