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
 
#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_TYPEDEFS__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_TYPEDEFS__

#include "../../common/code_marker.h" //error()

namespace ug{

#ifndef HAVE_BOOL
#define HAVE_BOOL

class BOOL{
public:
	BOOL() : value_(bool()){}
	/* explicit */ BOOL(bool const& t): value_(t) {}
	// /* explicit */ operator bool&() { return value_; }
	/* explicit */ operator bool() const { return value_; }
private:
	char value_;
};

#endif

}

//------------------

#include <boost/graph/adjacency_list.hpp>

namespace ug{

// TODO: think about edge set type
/* boost graph type for minimum degree ordering */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> graph_bidirectional_t;

typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeWeightProperty> graph_bidirectional_weighted_t;


/* boost graph type for boost Cuthill-McKee */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
	boost::property<boost::vertex_color_t, 
			 boost::default_color_type, 
			 boost::property<boost::vertex_degree_t, int> > > 
				Graph_t;

}

//------------------

//#include "metagraph/MatrixGraph.h"
#include "metagraph/WeightedMatrixGraph.h"
#include "metagraph/WeightedCliqueMatrixGraph.h"

namespace ug{

typedef graph_bidirectional_t base_graph_type;
typedef graph_bidirectional_weighted_t weighted_base_graph_type;
typedef std::vector<size_t> ordering_container_type;


template <typename TAlgebra>
using weighted_graph_type = WeightedMatrixGraph<weighted_base_graph_type, typename TAlgebra::matrix_type>;

template <typename TAlgebra>
using weighted_clique_graph_type = WeightedCliqueMatrixGraph<weighted_base_graph_type, typename TAlgebra::matrix_type>;

}

//---------------------

#include "algorithms/weighted_cuthill_mckee_ordering.cpp"
//#include "algorithms/boost_shortest_paths_ordering.cpp"

namespace ug{

template <typename TAlgebra>
using WeightedCuthillMcKeeOrdering_type = WeightedCuthillMcKeeOrdering<typename TAlgebra::matrix_type, weighted_base_graph_type, ordering_container_type>;

#if 0
template <typename TAlgebra>
using BoostShortestPathsOrdering_type = BoostShortestPathsOrdering<typename TAlgebra::matrix_type, weighted_base_graph_type, ordering_container_type>;
#endif

}

#endif //guard
