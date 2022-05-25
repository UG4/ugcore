/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Author: Felix Salfelder, 2022
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

#ifndef UG_GRAPH_INTERFACE_UNDIRECTED_H
#define UG_GRAPH_INTERFACE_UNDIRECTED_H

//#include "trace.h"
#include "sparsematrix_boost.h"
#include <boost/geometry/iterators/concatenate_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "common/util/bucket_sorter.hpp"
#include "common/stopwatch.h"

namespace ug{

// Symmetrify matrix stencil
template<class T>
class UndirectedMatrix {
	typedef typename T::const_row_iterator const_row_iterator;
public: // types
	typedef typename T::value_type value_type;
	typedef boost::adjacency_list<boost::vecS, boost::vecS> G;
	typedef typename boost::graph_traits<T>::adjacency_iterator T_adj_it;
	typedef typename boost::graph_traits<G>::adjacency_iterator G_adj_it;
	class edge : public std::pair<int, int>{
	public:
		explicit edge() : std::pair<int, int>(){ untested();
		}
		edge(int a, int b) : std::pair<int, int>(a, b) {
		}
		edge(boost::SM_edge<value_type>){ untested();
			incomplete();
		}
		edge(boost::graph_traits<G>::edge_descriptor){ untested();
			incomplete();
		}
	};
	typedef boost::geometry::concatenate_iterator<
	                             boost::SM_adjacency_iterator<value_type>,
	                             G_adj_it,
										  int, int> adjacency_iterator;

public:
	explicit UndirectedMatrix(T const* m=nullptr)
	    : _matrix(m) {
		if(m){
			refresh();
		}else{ untested();
		}
	}
	UndirectedMatrix(UndirectedMatrix && o) = delete;
	explicit UndirectedMatrix(UndirectedMatrix const& o)
	    : _matrix(o._matrix), _extra_fill(o._extra_fill) { untested();
	}
	UndirectedMatrix& operator=(UndirectedMatrix const& o) { untested();
		_matrix = o._matrix;
		_extra_fill = o._extra_fill;
		return *this;
	}

public: // interface
	void refresh();

	int num_rows() const {
		return std::max(_matrix->num_rows(), _matrix->num_cols());
	}
	int num_cols() const { untested();
		return num_rows();
	}
	int out_degree(int v) const {
		// could call in_degree?
		return in_degree(v);
	}
	int in_degree(int v) const {
		return boost::out_degree(v, *_matrix)
		     + boost::out_degree(v, _extra_fill);
	}
	int degree(int v) const {
		return out_degree(v);
	}

#if 0
	edge_iterator begin_row(int row) const { untested();
		auto r1 = boost::out_edges(row, *_matrix);
		auto r2 = boost::out_edges(row, _extra_fill);

		return edge_iterator(r1.first, r1.second, r2.first, r2.second);
	}
	edge_iterator end_row(int row) const { untested();
		auto r1 = boost::out_edges(row, *_matrix);
		auto r2 = boost::out_edges(row, _extra_fill);

		return edge_iterator(r1.second, r1.second, r2.second, r2.second);
	}
#endif
	adjacency_iterator begin_row(int row) const {
		auto r1 = boost::adjacent_vertices(row, *_matrix);
		auto r2 = boost::adjacent_vertices(row, _extra_fill);

		return adjacency_iterator(r1.first, r1.second, r2.first, r2.first);
	}
	adjacency_iterator end_row(int row) const {
		auto r1 = boost::adjacent_vertices(row, *_matrix);
		auto r2 = boost::adjacent_vertices(row, _extra_fill);

		return adjacency_iterator(r1.second, r2.first, r2.second);
	}

private:
	T const* _matrix;
	G _extra_fill;

private:
	class map_type{
	public:
		explicit map_type(T_adj_it const* i, T_adj_it const* e) : _it(i), _end(e){}

		int operator[](int i) const{
			assert(_it[i] != _end[i]);
			return *_it[i];
		}

	private:
		T_adj_it const* _it;
		T_adj_it const* _end;
	};
}; // UndirectedMatrix

// refresh missing fill. This is O(nonzeroes)
template<class T>
void UndirectedMatrix<T>::refresh()
{
	double t = get_clock_s();

	assert(_matrix);
	size_t N = _matrix->num_rows();
	assert(N == _matrix->num_cols());
	assert(!_matrix->iters()); // for now.

	_extra_fill = G(N);

	typedef boost::bucket_sorter<unsigned, unsigned,
			  map_type, boost::identity_property_map > container_type;

	std::vector<T_adj_it> c(N);
	std::vector<T_adj_it> e(N);
	map_type M(&c[0], &e[0]);

	container_type bs(N, N, M);

	for(size_t k=0; k<N; ++k){
		auto rr = boost::adjacent_vertices(k, *_matrix);
		c[k] = rr.first;
		e[k] = rr.second;

		if(c[k] != e[k]){
			bs.push(k);
		}else{
		}
	}

#if 0
	for(size_t j=0; j<N; ++j){
		std::cerr << "bs init " << j << "... ";
		for(auto t : bs[j]){
			std::cerr << " " << t;
		}
		std::cerr << "\n";
	}
#endif

	for(size_t j=0; j<N; ++j){
		// std::cerr << "====== " << j << " ====\n";
		if(c[j] == e[j]){
		}else if(*c[j]==j){
			++c[j];

			if(c[j] == e[j]){ itested();
				bs.remove(j);
			}else{
				bs.update(j);
			}

		}else{
			assert(*c[j]>j);
		}

		// upper triangle, visit nonzeroes.
		for(; c[j] != e[j];){
			int k = *c[j];
			if(c[k] == e[k]){
				assert(size_t(k)<N);
				assert(j<N);
				// std::cerr << "lower fill0 " << k << " " << j << "\n";
				boost::add_edge(k, j, _extra_fill);
			}else if(*c[k] > j){
				assert(size_t(k)<N);
				assert(j<N);
				// std::cerr << "lower fill1 " << k << " " << j << "\n";
				boost::add_edge(k, j, _extra_fill);
			}else if(*c[k] == j){
				++c[k];
				if(c[k] == e[k]){
					bs.remove(k);
				}else{
					bs.update(k);
				}
			}
			++c[j];

			if(c[j] == e[j]){
				bs.remove(j);
			}else{
				bs.update(j);
			}
		}

		// lower triangle, visit subset of nonzeroes.
		// std::cerr << "====== lower " << j << " ====\n";
		assert(j<N);
		auto bj = bs[j];

		for(auto t=bj.begin(); t!=bj.end(); ){
			int i = *t;
			++t;

			// std::cerr << "upper fill " << j << " " << i << "\n";
			boost::add_edge(j, i, _extra_fill);

			assert(*c[i]==j);

			++c[i];
			if(c[i] == e[i]){
				bs.remove(i);
			}else if(int(*c[i]) < i){ itested();
				// std::cerr << "bs push " << i << " " << j << " " << *c[i] << "\n";
				assert(*c[i] > j);
				bs.update(i);
			}else{ itested();
				assert(*c[i] > j);
				// bs.remove(i);
			}

		}
	}

	c.clear();
	e.clear();
	assert(!_matrix->iters());

	UG_LOG("Undirected refresh " << get_clock_s()-t << "\n");
}

} // ug

#endif // guard
