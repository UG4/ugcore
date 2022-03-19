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

#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_FOLLOW_CONVECTION_ORDERING__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_FOLLOW_CONVECTION_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp> //TODO: remove if not used

#include <set>
#include <algorithm> //reverse
#include <utility> //pair

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "bindings/lua/lua_user_data.h"

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"

#include <assert.h>
#include "common/error.h"


namespace ug{

bool CompareScalar(const std::pair<number, size_t> &p1,
		const std::pair<number, size_t> &p2)
{
	return p1.first<p2.first;
}

//for sorting
struct Blo{
	size_t v;
	number w;
};

bool compBlo(Blo a, Blo b){
	return a.w < b.w;
}

template <typename TAlgebra, typename TDomain, typename O_t>
class FollowConvectionOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	typedef MathVector<TDomain::dim> small_vec_t;
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;

	typedef typename std::pair<MathVector<TDomain::dim>, size_t> Position_t;
	typedef typename std::pair<number, size_t> Scalar_t;
	typedef UserData<MathVector<TDomain::dim>, TDomain::dim> Velocity_t;

	typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty; //TODO: number
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeWeightProperty> G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd_t;
	typedef typename boost::graph_traits<G_t>::edge_descriptor ed_t;
	typedef typename std::pair<ed_t, bool> edge_t;
	typedef typename boost::graph_traits<G_t>::vertex_iterator vIt_t;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator nIt_t;
	typedef typename boost::graph_traits<G_t>::in_edge_iterator in_eIt_t;
	typedef typename boost::graph_traits<G_t>::out_edge_iterator out_eIt_t;

	//TODO: use other container
	typedef std::set<vd_t> front_type;

	FollowConvectionOrdering() : m_ssDirichletIdx(-1){}

	/// clone constructor
	FollowConvectionOrdering( const FollowConvectionOrdering<TAlgebra, TDomain, O_t> &parent )
			: baseclass(), m_ssDirichletIdx(parent.m_ssDirichletIdx){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new FollowConvectionOrdering<TAlgebra, TDomain, O_t>(*this));
	}

	/* TODO: idea, is dijkstras algorithm on graph with angle weights beginning at s, s connected to dirichlet nodes solving this as well? */

	//overload
	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}

		//TODO: do this in init at graph creation?
		//compute angle between velocity and dist(v, w)
		vIt_t vIt, vEnd;
		nIt_t nIt, nEnd;
		vd_t s, t;
		small_vec_t pos_s, pos_t, vel_s, dir_st;
		edge_t edge;
		number angle;
		number prod;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
			s = *vIt;
			pos_s = m_vPositions.at(s).first;

			//call lua function, store velocity in s in vel_s
			(*m_spVelocity)(vel_s, pos_s, 0.0f, 0);

			if(VecLengthSq(vel_s) > 0){
				for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, g); nIt != nEnd; ++nIt){
					t = *nIt;
					pos_t = m_vPositions.at(t).first;
					VecSubtract(dir_st, pos_t, pos_s);
					angle = VecAngle(dir_st, vel_s);
					//prod = dir_st*vel_s;

					UG_LOG("prod: " << prod << "\n");

					edge = boost::edge(s, t, g);
					//boost::put(boost::edge_weight_t(), g, edge.first, prod);
					//boost::put(boost::edge_weight_t(), g, edge.first, angle);
					boost::put(boost::edge_weight_t(), g, edge.first, VecLength(dir_st));

					//UG_LOG("edge " << pos_s << " -> " << pos_t << ", w= " << angle << "\n");
				}
			}
			else{
				UG_LOG(name() << "::compute:: velocity at " << pos_s << " is 0!\n");
				for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, g); nIt != nEnd; ++nIt){
					t = *nIt;
					edge = boost::edge(s, t, g);
					boost::put(boost::edge_weight_t(), g, edge.first, 4*PI);
				}
			}
		}


		//find start vertices
		m_vScalars.resize(m_vPositions.size());
		small_vec_t pos;
		small_vec_t vel;
		(*m_spVelocity)(vel, m_vPositions[0].first, 0.0f, 0); //TODO: this assumes vel is constant!

		for(size_t i = 0; i < m_vPositions.size(); ++i){
			pos = m_vPositions[i].first;
			m_vScalars[i].first = pos*vel; //scalar product
			m_vScalars[i].second = m_vPositions[i].second;
		}

		std::sort(m_vScalars.begin(), m_vScalars.end(), CompareScalar);

		number min_w = m_vScalars[0].first;
		size_t min_v;

		s = boost::add_vertex(g);

		for(size_t i = 0; i < m_vScalars.size(); ++i){
			if(m_vScalars[i].first > min_w){
				break;
			}

			min_v = m_vScalars[i].second;

			boost::add_edge(s, min_v, 0.0f, g);

			UG_LOG(name() << "::compute: min_v = " << min_v << " (" << m_vPositions[min_v].first << "), min_w = " << min_w << "\n");
		}

		std::vector<vd_t> p(n+1); //parents
		std::vector<number> d(n+1); //distances

		//start vertex
		//s = *boost::vertices(g).first;

		UG_LOG("distance to s: " << d[s] << "\n");

		boost::dijkstra_shortest_paths(g, s, boost::predecessor_map(&p[0]).distance_map(&d[0]));

		UG_LOG(name() << "::compute: Done!\n");

		std::vector<Blo> blo(n);
		for(unsigned i = 0; i < n; ++i){
			UG_LOG("d: " << d[i] << "\n");
			blo[i].v = i;
			blo[i].w = d[i];
		}

		//sort o according to d
		std::sort(blo.begin(), blo.end(), compBlo);

		o.resize(n);
		for(unsigned i = 0; i < n; ++i){
			o[blo[i].v] = i;
			//o[i] = blo[i].v;
			UG_LOG("o[" << i << "]: " << m_vPositions[blo[i].v].first << ", cost: " << d[blo[i].v] << "\n");
		}

/*
		for(size_t i = 0; i < m_vScalars.size(); ++i){
			o[m_vScalars[i].second] = i;
		}
*/
/*
		//find start vertex
		s = *(boost::vertices(g).first);
		pos_s = m_vPositions.at(s).first;
		UG_LOG("start vertex " << pos_s << "\n");
		bool found;
		while(true){
			found = false;
			for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(s, g); nIt != nEnd; ++nIt){
				edge = boost::edge(s, *nIt, g);
				double w = boost::get(boost::edge_weight_t(), g, edge.first);
				if(PI-w < PI/2){
					s = *nIt;
					pos_s = m_vPositions.at(s).first;
					UG_LOG("new start vertex " << pos_s << "\n");
					found = true;
					break;
				}
			}
			if(!found){
				break;
			}
		}

		m_sFront.insert(s);
		m_dirichlet[s] = true;
*/

/*
		//choose edge (v, w) in cut(m_sFront, V(g)/m_sFront) with minimum weight
		//add w to m_sFront
		//ensure that are no edges between vertices in m_sFront
		UG_LOG(name() << "::compute: front.size(): " << m_sFront.size() << "\n");

		double w;
		vd_t min_source; //TODO: just for sanity check
		vd_t min_target;
		size_t min_idx;
		size_t k = m_sFront.size();
		//TODO: if same value, take earlier numbered vertex
		while(m_sFront.size() < n){
			min_w = 2*PI; //PI should be sufficient
			min_idx = n;
			for(front_type::iterator sIt = m_sFront.begin(); sIt != m_sFront.end(); ++sIt){
				for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*sIt, g); nIt != nEnd; ++nIt){
					if(m_dirichlet[*nIt]){
						continue;
						//UG_THROW(name() << "::compute: iterated over a dirichlet neighbor!\n");
					}

					edge = boost::edge(*sIt, *nIt, g);
					w = boost::get(boost::edge_weight_t(), g, edge.first);
					if(w < min_w){
						min_w = w;
						min_source = *sIt;
						min_target = *nIt;
						min_idx = o[min_source];
					}
					else if(w == min_w){
						if(o[*sIt] < min_idx){
							min_idx = o[*sIt];
							min_source = *sIt;
							min_target = *nIt;
						}
					}
					else{}
				}
			}
			UG_LOG("add vertex " << min_target << " via edge " << min_source << " -> " << min_target << ", w = " << min_w << "\n");

#ifdef WHATEVER
			for(front_type::iterator sIt = m_sFront.begin(); sIt != m_sFront.end(); ++sIt){
				for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*sIt, g); nIt != nEnd; ++nIt){
					if(m_dirichlet[*nIt]){
						continue;
					}
					edge = boost::edge(*sIt, *nIt, g);
					w = boost::get(boost::edge_weight_t(), g, edge.first);
					UG_LOG("w: " << w << "\n");
				}
			}
#endif
			m_sFront.insert(min_target);
			//o[k] = min_target;
			o[min_target] = k;
			m_dirichlet[min_target] = true; //use this container as visited marker
			++k;
		}
*/
		g = G_t(0);

		#ifdef UG_DEBUG
		check();
		#endif
	}

	void init(M_t* A, const V_t& V){
		const GridFunc_t* pGridF;
		unsigned n;

		try{
			if((pGridF = dynamic_cast<const GridFunc_t*>(&V)) == 0){
				UG_THROW(name() << "::init: No DoFDistribution specified.");
			}

			SmartPtr<DoFDistribution> dd = ((GridFunc_t*) pGridF)->dof_distribution();

			m_ssDirichletIdx = pGridF->domain()->subset_handler()->get_subset_index(m_ssDirichletName);

			if(m_ssDirichletIdx < 0)
				UG_THROW(name() << "::init: Invalid subset for sources selected! Call 'select_sources(const char*)'.");

			n = dd->num_indices();

			if(n != A->num_rows ()){
				UG_THROW(name() << "::init: #indices != #rows");
			}

			m_vvConnection.resize(n);
			try{
				dd->get_connections<ug::Edge>(m_vvConnection);
			}
			UG_CATCH_THROW(name() << "::init: No adjacency graph available!");

			ExtractPositions(pGridF->domain(), dd, m_vPositions);
		}
		catch(...){
			throw;
		}

		o.resize(n);

		//select dirichlet vertices according to m_ssDirichletIdx
		typedef typename GridFunc_t::template traits<ug::Vertex>::const_iterator ugVertIt_t;
		small_vec_t pos_s, pos_t;
		size_t numDirichlet = 0;
		m_dirichlet = std::vector<BOOL>(n, false);
		ugVertIt_t ugVertIt = pGridF->template begin<ug::Vertex>();
		ugVertIt_t ugVertEnd = pGridF->template end<ug::Vertex>();
		size_t k = 0;
		for(; ugVertIt != ugVertEnd; ++ugVertIt, ++k){
			ug::Vertex* v = *ugVertIt;

			int si = pGridF->domain()->subset_handler()->get_subset_index(v);

			if(si == m_ssDirichletIdx){
				//o[numDirichlet] = k;
				//o[k] = numDirichlet;
				//m_sFront.insert(k);
				//m_dirichlet[k] = true;
				++numDirichlet;
			}
		}

		//create graph, do not add edges (v, w) with w being a dirichlet vertex
		g = G_t(n);

/*
		for(unsigned i = 0; i < n; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){
					if(!m_dirichlet[conn.index()]){
						boost::add_edge(i, conn.index(), g);
						//boost::add_edge(i, conn.index(), w, g); //TODO: compute angle here?
					}
				}
			}
		}
*/
		for(unsigned i = 0; i < n; ++i){
			for(unsigned j = 0; j < m_vvConnection[i].size(); ++j){
				if(i != m_vvConnection[i][j]){ //no self loops!
					boost::add_edge(i, m_vvConnection[i][j], g);
				}
			}
		}

/*
		//create grid edges {v, w} for dirichlet vertices v, w must not be a dirichlet vertex
		for(front_type::iterator sIt = m_sFront.begin(); sIt != m_sFront.end(); ++sIt){
			for(size_t i = 0; i < m_vvConnection[*sIt].size(); ++i){
				if(m_vvConnection[*sIt][i] == *sIt){
					continue;
				}
				if(!m_dirichlet[m_vvConnection[*sIt][i]]){
					boost::add_edge(*sIt, m_vvConnection[*sIt][i], g);
					UG_LOG("create dirichlet edge: " << *sIt << " -> " << m_vvConnection[*sIt][i] << "\n");
				}
			}
		}
*/

		if(numDirichlet == 0){
			UG_THROW(name() << "::init: #dirichlet_nodes = 0!");
		}

		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " (subset " << m_ssDirichletIdx << ", " << m_ssDirichletName
				<< ", n=" << boost::num_vertices(g) << ", m=" << boost::num_edges(g)
				<< ", d=" << numDirichlet << ")\n");
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

	void select_dirichlet_subset(const char* ssDirichletName){
		m_ssDirichletName = ssDirichletName;
	}

	void set_velocity(const char* strVelocity){
		m_spVelocity = make_sp(new LuaUserData<MathVector<TDomain::dim>, TDomain::dim>(strVelocity));
	}

	virtual const char* name() const {return "FollowConvectionOrdering";}

private:
	G_t g;
	O_t o;

	SmartPtr<Velocity_t> m_spVelocity;

	int m_ssDirichletIdx;
	const char* m_ssDirichletName;

	front_type m_sFront;
	std::vector<BOOL> m_dirichlet;

	std::vector<Position_t> m_vPositions;
	std::vector<Scalar_t> m_vScalars;
	std::vector<std::vector<size_t> > m_vvConnection;
};

} //namespace


#endif //guard
