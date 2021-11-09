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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_SHORTEST_PATHS_ORDERING__BLA
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_SHORTEST_PATHS_ORDERING__BLA

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "IOrderingAlgorithm.h"
#include "util.cpp"
#include "common/error.h"

namespace ug{

//for sorting
struct Blo{
	size_t v;
	double w;
};

bool compBlo(Blo a, Blo b){
	return a.w < b.w;
}


template <typename TAlgebra, typename O_t>
class BoostShortestPathsOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeWeightProperty> G_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	BoostShortestPathsOrdering(){}

	/// clone constructor
	BoostShortestPathsOrdering( const BoostShortestPathsOrdering<TAlgebra, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new BoostShortestPathsOrdering<TAlgebra, O_t>(*this));
	}

	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}

		typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
		std::vector<vd> p(n); //parents
		std::vector<int> d(n); //distances

		vd s = *boost::vertices(g).first; //start vertex	//TODO: choose a vertex strategically

		boost::dijkstra_shortest_paths(g, s, boost::predecessor_map(&p[0]).distance_map(&d[0]));

		std::vector<Blo> blo(n);
		for(unsigned i = 0; i < n; ++i){
			blo[i].v = i;
			blo[i].w = d[i];
		}

		//sort o according to d
		std::sort(blo.begin(), blo.end(), compBlo);

		o.resize(n);
		for(unsigned i = 0; i < n; ++i){
			o[i] = blo[i].v;
		}

		g = G_t(0);
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
		UG_LOG("Using " << name() << "\n");
		unsigned rows = A->num_rows();

		g = G_t(rows);

		for(unsigned i = 0; i < rows; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){ //TODO: think about this!!
					double w;
					w = abs(conn.value()); //TODO: think about this
					boost::add_edge(i, conn.index(), w, g);
				}
			}
		}
	}

	virtual const char* name() const {return "BoostShortestPathsOrdering";}

private:
	G_t g;
	O_t o;
};

} //namespace

#endif //guard
