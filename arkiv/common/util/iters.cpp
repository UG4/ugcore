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

#include <tuple>
#include <vector>

namespace ug {

#ifndef HAVE_BOOL
#define HAVE_BOOL

class BOOL{
public:
	BOOL() : value_(bool()){}
	/* explicit */ BOOL(bool const& t): value_(t) {}
	// /* explicit */ operator bool &() { return value_; }
	/* explicit */ operator bool () const { return value_; }
private:
	char value_;
};

#endif

class unvisited_iterator{
public:
	unvisited_iterator(std::vector<BOOL> &vec_bool, unsigned i) : n(vec_bool.size()), cur(i), visited(vec_bool){
	}

	unvisited_iterator(std::vector<BOOL> &vec_bool) : n(vec_bool.size()), cur(0), visited(vec_bool){
		if(visited[cur]){
			operator ++ ();
		}
	}

	bool operator == (const unvisited_iterator& o) const{
		return cur==o.cur;
	}

	bool operator != (const unvisited_iterator& o) const{
		return !operator == (o);
	}

	void operator ++ (){
		do{
			++cur;
		}
		while(cur < n && visited[cur]);
	}

	unsigned operator * () const{
		return cur;
        }

public:
	unsigned n;
	unsigned cur;
	std::vector<BOOL> &visited;
};

#if 0
std::tuple<unvisited_iterator, unvisited_iterator> make_unvisited_iterator(std::vector<BOOL> &visited)
{
	auto a = unvisited_iterator(visited, 0);
	auto b = unvisited_iterator(visited, visited.size());

	return std::make_tuple(a, b);
}
#endif


/* Redundant in case that ug::Selector can do the job... */

template <typename TAlgebra>
class dirichlet_iterator{
public:
	using matrix_type = typename TAlgebra::matrix_type;

	dirichlet_iterator(matrix_type &mat, unsigned i, unsigned num) : A(mat), cur(i), n(num){
		operator ++ ();
	}

	bool operator == (const dirichlet_iterator& o) const{
		return cur==o.cur;
	}

	bool operator != (const dirichlet_iterator& o) const{
		return !operator == (o);
	}

	void operator ++ (){
		do{
			++cur;
		}
		while(cur < n && !A.is_isolated(cur));
	}

	unsigned operator * () const{
		return cur;
        }

protected:
	matrix_type &A;
	unsigned cur;
	unsigned n;
};


template <typename TAlgebra>
class non_dirichlet_iterator{
public:
	using matrix_type = typename TAlgebra::matrix_type;

	non_dirichlet_iterator(matrix_type &mat, unsigned i, unsigned num) : A(mat), cur(i), n(num){
		operator ++ ();
	}

	bool operator == (const non_dirichlet_iterator& o) const{
		return cur==o.cur;
	}

	bool operator != (const non_dirichlet_iterator& o) const{
		return !operator == (o);
	}

	void operator ++ (){
		do{
			++cur;
		}
		while(cur < n && A.is_isolated(cur));
	}

	unsigned operator * () const{
		return cur;
        }

protected:
	matrix_type &A;
	unsigned cur;
	unsigned n;
};



template <typename TAlgebra>
std::tuple<dirichlet_iterator<TAlgebra>, dirichlet_iterator<TAlgebra> > make_dirichlet_iterator(typename TAlgebra::matrix_type &mat)
{
	auto a = dirichlet_iterator<TAlgebra>(mat, 0, mat.num_rows());
	auto b = dirichlet_iterator<TAlgebra>(mat, mat.num_rows(), mat.num_rows());

	return std::make_tuple(a, b);
}

template <typename TAlgebra>
std::tuple<non_dirichlet_iterator<TAlgebra>, non_dirichlet_iterator<TAlgebra> > make_non_dirichlet_iterator(typename TAlgebra::matrix_type &mat)
{
	auto a = non_dirichlet_iterator<TAlgebra>(mat, 0, mat.num_rows());
	auto b = non_dirichlet_iterator<TAlgebra>(mat, mat.num_rows(), mat.num_rows());

	return std::make_tuple(a, b);
}

} //namespace

