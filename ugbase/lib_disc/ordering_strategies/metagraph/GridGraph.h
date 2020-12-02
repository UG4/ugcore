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


#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_METAGRAPH_GRID_GRAPH__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_METAGRAPH_GRID_GRAPH__

#include "../../../lib_algebra/ordering_strategies/metagraph/IMetaGraph.h"

namespace ug{

/*
	Creates a graph g of type G_t based on a grid grid (accessed via ApproximationSpace->dof_distribution->get_connections). g is accessible via graph().
	The graph g is equivalent to the grid graph.
*/

//TODO: get the edges directly via approx->domain->grid()
//TODO: avoid copy via get_connections(), rather use an iterator, write one of not existent.
template <typename TDomain, typename G_t>
class GridGraph : public IMetaGraph<G_t>{
public:
	typedef IMetaGraph<G_t> baseclass;

	ConnectionGraph(SmartPtr<ug::ApproximationSpace<TDomain> > approx, unsigned grid_level){
		std::vector<SmartPtr<ug::DoFDistribution> > vDD = approxSpace->dof_distributions();
		auto dofDistr = vDD[grid_level];
		unsigned n = dofDistr->num_indices();

		std::vector<std::vector<size_t> > vvConnection;
		try{
			dofDistr->get_connections<ug::Edge>(vvConnection); //this the element-graph
		}
		UG_CATCH_THROW("[GridGraph] No adjacency graph available.");

		baseclass::g = G_t(n);

		for(unsigned i = 0; i < n; ++i){
			for(unsigned j = 1; j < vvConnection[i].size(); ++j){ //first element is the index of the vertex itself
				if(i != vvConnection[i][j]){ //no self loops!
					boost::add_edge(i, vvConnection[i][j], baseclass::g);
				}
			}
		}
	}
};

} //namespace

#endif //guard

