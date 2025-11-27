/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef UNSORTED_VECTOR_H_
#define UNSORTED_VECTOR_H_

#include <vector>
#include "connection.h"

namespace ug{

/**
 * This is in most cases faster than the std::map-based SparseVector in sparse_vector.h
 * because it uses and "posInConnection" array which reduces all operations to O(1),
 * instead of O(log n) for std::map when number of non-zeroes in the vector.
 * The additional array is as long as the (non-sparse) size of the vector, so this makes only sense if you
 * REUSE the UnsortedSparseVector, like
 * \code
 * N = 10000; // N is the non-sparse total size of the vector
 * UnsortedSparseVector<int> vec(N); // expensive.
 * for(size_t i=0; i < N; i++)
 * {
 * 		vec.clear();
 * 		for(size_t j=0; j<50; j++)
 * 			vec(rand()%N)++;	// O(1)
 * }
 * \endcode
 */
template<typename TValue>
class UnsortedSparseVector
{
public:
	using value_type = TValue;
	using connection = AlgebraicConnection<TValue>;

	using iterator = typename std::vector<connection>::iterator;
	using const_iterator = typename std::vector<connection>::const_iterator;

private:
	std::vector<int> posInConnections;
	std::vector<connection > con;
	size_t m_size;
public:

	explicit UnsortedSparseVector(size_t s) : posInConnections(s, -1), m_size(s)
	{
		con.reserve(32);
	}
	iterator begin()
	{
		return con.begin();
	}
	iterator end()
	{
		return con.end();
	}
	const_iterator begin() const
	{
		return con.begin();
	}
	const_iterator end() const
	{
		return con.end();
	}

	[[nodiscard]] size_t num_connections() const
	{
		return con.size();
	}

	[[nodiscard]] size_t size() const
	{
		return m_size;
	}

	connection *unsorted_raw_ptr()
	{
		return &con[0];
	}


	void clear()
	{
		for(size_t i=0; i<con.size(); i++)
			posInConnections[con[i].iIndex] = -1;
		con.clear();
	}

	const value_type &operator () (size_t c) const
	{
		assert(c < m_size);
		int p = posInConnections[c];
		if(p != -1)
			return con[p].value();
		else
			assert(0 && "const_and_not_available");
	}
	value_type &operator () (size_t c)
	{
		assert(c < m_size);
		int p = posInConnections[c];
		if(p != -1)
			return con[p].value();
		else
		{
			p = posInConnections[c] = con.size();
			con.push_back(connection(c, TValue(0)));
			return con[p].value();
		}
	}
	[[nodiscard]] bool has_connection(size_t c) const
	{
		assert(c < m_size);
		return posInConnections[c] != -1;
	}
};


}
#endif