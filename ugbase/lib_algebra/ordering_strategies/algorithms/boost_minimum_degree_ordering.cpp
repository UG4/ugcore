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
 
#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_MINIMUM_DEGREE_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_MINIMUM_DEGREE_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/minimum_degree_ordering.hpp>

#include "IOrderingAlgorithm.h"
#include "util.cpp"

//debug
#include "common/error.h"
#include "common/log.h"

namespace ug{

//Important Note: This implementation requires the BGL graph to be
//directed.  Therefore, nonzero entry (i, j) in a symmetrical matrix
//A coresponds to two directed edges (i->j and j->i).
template <typename TAlgebra, typename O_t>
class BoostMinimumDegreeOrdering final : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> G_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	BoostMinimumDegreeOrdering(){}

	/// clone constructor
	BoostMinimumDegreeOrdering( const BoostMinimumDegreeOrdering<TAlgebra, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new BoostMinimumDegreeOrdering<TAlgebra, O_t>(*this));
	}

	void compute(){
		unsigned n = boost::num_vertices(g);
		unsigned e = boost::num_edges(g);

		O_t io(n, 0);

		o.resize(n);
		unsigned i = 0;

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}
		else if(n*(n-1u)==e || e==0){
			UG_LOG(name() << "::compute: Graph is complete or edgeless!");
			typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
			for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; vIt++){
				o[i++] = *vIt;
			}
		}

		std::vector<int> inverse_perm(n, 0);
		std::vector<int> supernode_sizes(n, 1);
		auto id = boost::get(boost::vertex_index, g);
		std::vector<int> degree(n, 0);

		/*
		 * (Graph& G,
		 *  DegreeMap degree,
		 *  InversePermutationMap inverse_perm,
		 *  PermutationMap perm,
		 *  SuperNodeMap supernode_size,
		 *  int delta,
		 *  VertexIndexMap vertex_index_map)
		 */

		boost::minimum_degree_ordering
		  (g,
		   boost::make_iterator_property_map(&degree[0], id, degree[0]),
		   &io[0],
		   &(o[0]),
		   boost::make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
		   0,
		   id
		   );

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
					boost::add_edge(i, conn.index(), g);
				}
			}
		}
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	virtual const char* name() const {return "BoostMinimumDegreeOrdering";}

private:
	G_t g;
	O_t o;
};

} //namespace

#endif //guard
