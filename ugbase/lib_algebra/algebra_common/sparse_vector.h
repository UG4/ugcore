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

#ifndef SPARSE_VECTOR_H_
#define SPARSE_VECTOR_H_

#include <map>
namespace ug{

template<typename T>
class SparseVector
{
	size_t m_size;
	using container = std::map<size_t, T>;
	container data;
public:
	using value_type = T;

	class const_iterator : public container::const_iterator
	{
		using container::const_iterator::operator *;
	public:
		const_iterator(typename container::const_iterator it) : container::const_iterator(it) {}
		const T &value() const { return (operator *()).second; }
		size_t index() const { return (operator *()).first; }
	};

	class iterator : public container::iterator
	{
		using container::iterator::operator *;
	public:
		iterator(typename container::iterator it) : container::iterator(it) {}
		const T &value() const { return (operator *()).second; }
		T &value() { return (operator *())->second; }
		size_t index() const { return (operator *()).first; }
	};

	SparseVector(size_t s) : m_size(s)
	{
	}
	const_iterator begin() const
	{
		const_iterator c(data.begin());
		return c;
	}
	const_iterator end() const
	{
		return const_iterator(data.end());
	}

	const T &operator()(size_t c) const
	{
		assert(c < m_size);
		return data[c];
	}
	T &operator()(size_t c)
	{
		assert(c < m_size);
		return data[c];
	}
	bool has_connection(size_t c) const
	{
		return data.find(c) != data.end();
	}

	size_t size()
	{
		return m_size;
	}
	size_t num_connections() const
	{
		return data.size();
	}

	void print() const
	{
		for(const_iterator it=begin(); it != end(); ++it)
		{
			//if(BlockNorm(it.value()) == 0.0) continue;
			UG_LOG("(" << it.index() << " -> " << it.value() << ")");
		}

		UG_LOG("\n");
	}

};

}
#endif