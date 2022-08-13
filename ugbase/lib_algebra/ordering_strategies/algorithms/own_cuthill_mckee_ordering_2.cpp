/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Author: Lukas Larisch, Felix Salfelder
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

#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_OWN_CUTHILL_MCKEE_ORDERING2__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_OWN_CUTHILL_MCKEE_ORDERING2__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <set>
#include <algorithm> //reverse
#include <utility> //pair
#include <deque>
#include <list>

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"

#include <assert.h>
#include "common/error.h"
#include "common/util/bucket_sorter.hpp"

#define TRACE_UNTESTED
#include "common/util/trace.h"

// is this already implemented somewhere else?
#ifndef HAVE_MYABS
#define HAVE_MYABS
namespace{
template<class T>
double my_abs(T){return 0;}

template<>
double my_abs(double v){return abs(v);}
}
#endif

namespace ug{

static unsigned INVALID = -1u;

template <typename TAlgebra, typename O_t>
class OwnCuthillMcKeeOrdering2 : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef unsigned deg_t;
	typedef typename boost::graph_traits<G_t>::edge_descriptor ed_t;
	typedef typename boost::graph_traits<G_t>::vertex_iterator vIt_t;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;
	typedef typename boost::graph_traits<G_t>::out_edge_iterator oute_iter;
	typedef boost::iterator_property_map<unsigned*,
			boost::identity_property_map, unsigned, unsigned&> map_type;
	typedef boost::bucket_sorter<unsigned, unsigned,
			map_type, boost::identity_property_map > bucket_sorter;

	std::vector<unsigned> _degrees;

public:
	OwnCuthillMcKeeOrdering2() : m_bReverse(false), m_look_for_sources(true){}

	/// clone constructor
	OwnCuthillMcKeeOrdering2( const OwnCuthillMcKeeOrdering2<TAlgebra, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse), m_look_for_sources(parent.m_look_for_sources){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone() {
		return make_sp(new OwnCuthillMcKeeOrdering2<TAlgebra, O_t>(*this));
	}

	inline void unregister_indegree(size_t v, std::vector<size_t>& indegs){
		std::pair<oute_iter, oute_iter> e;
		e = boost::out_edges(v, g);

		for(; e.first != e.second; ++e.first){
			ed_t const& edg = *e.first;
			auto t = boost::target(edg, g);
			--indegs[t];
		 }
	}

	void compute() override { untested();
		unsigned n = boost::num_vertices(g);
		unsigned min_bucket = -1u;

		if(n == 0){ untested();
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}else{ untested();
		}

		o.resize(0);
		o.resize(n, INVALID);

		vd s = *boost::vertices(g).first;
		o[s] = 0;
		_degrees_[s] = INVALID;
		_tags[s] = 0;

		auto i = boost::adjacent_vertices(s, g);
		for(; i.first != i.second; ++i.first){ untested();
			vd n = *i.first;
			assert(n!=s);
			deg_t d = boost::in_degree(n, g);

			_degrees_[n] = d;
			_bs.push_front(n);

			assert(_tags[n] == INVALID);

			_tags[n] = s;
			min_bucket = std::min(_degrees_[n], min_bucket);
		}

		size_t k = 1;

		std::pair<oute_iter, oute_iter> e;

		//main loop
		for(; k < n;){ untested();
			unsigned min_tag = -1u;
			vd next = INVALID;
			for(unsigned j=0; j<min_bucket; ++j){ untested();
				assert(_bs[j].empty());
			}
			assert(!_bs[min_bucket].empty());

			for(vd cand : _bs[min_bucket]){ untested();
				unsigned tag = _tags[cand];
				assert(tag!=INVALID);
				if(tag < min_tag){ untested();
					next = cand;
					min_tag = tag;
				}else{ untested();
				}
			}

			assert(_degrees_[next] == min_bucket);
			assert(_tags[next] != INVALID);
			assert(next!=INVALID);

			auto i = boost::adjacent_vertices(next, g);
			for(; i.first != i.second; ++i.first){ untested();
				vd n = *i.first;
				assert(n!=next);

				if(_tags[n] == INVALID){
					// not a candidate, add it, initialise degree.

					deg_t d = boost::in_degree(n, g);
					assert(d);
					_tags[n] = k;
					
					_degrees_[n] = d;
					_bs.push_front(n);

					min_bucket = std::min(_degrees_[n], min_bucket);
				}else if(_degrees_[n] == INVALID){ untested();
					// n already numbered.
				}else{ untested();
					// already reachable from S. Degree drop.
					assert(_tags[n]<k);
					assert(_degrees_[n]);


					if(min_bucket == _degrees_[n]){ untested();
						--min_bucket;
					}else{ untested();
						assert(min_bucket < _degrees_[n]);
					}

					--_degrees_[n];
					_bs.update(n);
				}
			}

			_bs.remove(next);
			_degrees_[next] = INVALID;
			assert(o[next] == INVALID);
			o[next] = k;
			++k;

			if(k==o.size()){ untested();
			}else{ untested();
				for(unsigned i=0; i<min_bucket; ++i){
					assert(_bs[i].empty());
				}
				while(_bs[min_bucket].empty()){ untested();
					++min_bucket;
					assert(min_bucket < o.size()); // not connected?
				}
			}
		} // main loop

		g = G_t(0);
	}

	void init(M_t* A, const V_t&){ untested();
		init(A);
	}

	void init(M_t* A){ untested();
		unsigned n = A->num_rows();

		g = G_t(n);

		for(unsigned i = 0; i < n; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){
					boost::add_edge(conn.index(), i, g);
				}
			}
		}

		init_bs();
	}

	void init(M_t* A, const V_t&, const O_t& inv_map, const O_t& start){ untested();
		init(A, inv_map, start);
		init_bs();
	}

	void init(M_t* A, const O_t& inv_map, const O_t& start){ untested();
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " on induced matrix of size " << inv_map.size() << "\n");
		#endif

		own_cmk_induced_subgraph<G_t, M_t>(g, A, inv_map);
		m_look_for_sources = false;

		size_t n = boost::num_vertices(g);
		UG_LOG("n: " << n << ", e: " << boost::num_edges(g) << "\n");

		init_bs();
	}

	void init_bs(){ untested();
		size_t n = boost::num_vertices(g);
		assert(n);
		_degrees.resize(0);
		_degrees.resize(n, INVALID);
		_tags.resize(0);
		_tags.resize(n, INVALID);
		_degrees_ = map_type(&_degrees[0], boost::identity_property_map());
		_bs = bucket_sorter(n+1, n+1, _degrees_);
	}

	void check(){
		if(!is_permutation(o)){
			UG_THROW(name() << "::check: Not a permutation!");
		}
	}

	O_t& ordering(){ untested();
#ifndef NDEBUG
		check();
#endif
		return o;
	}

	void set_reverse(bool b){
		m_bReverse = b;
	}

	virtual const char* name() const {return "OwnCuthillMcKeeOrdering2";}

private:
	G_t g;
	O_t o;

	std::list<size_t> front;
	bool m_bReverse;
	bool m_look_for_sources;

	std::vector<unsigned> _tags;
	map_type _degrees_;
	bucket_sorter _bs; // (10, 10, M);
}; // OwnCuthillMcKeeOrdering2

} //namespace


#endif //guard
