/*
 * unsorted_vector.h
 *
 *  Created on: 13.05.2013
 *      Author: mrupp
 */

#ifndef UNSORTED_VECTOR_H_
#define UNSORTED_VECTOR_H_

#include <vector>
#include "connection.h"

namespace ug{

template<typename TValue>
class UnsortedSparseVector
{
public:
	typedef TValue value_type;
	typedef AlgebraicConnection<TValue> connection;

	typedef typename std::vector<connection>::iterator iterator;
	typedef typename std::vector<connection>::const_iterator const_iterator;

private:
	std::vector<int> posInConnections;
	std::vector<connection > con;
	size_t m_size;
public:

	UnsortedSparseVector(size_t s) : posInConnections(s, -1), m_size(s)
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

	size_t num_connections() const
	{
		return con.size();
	}

	size_t size() const
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

	const value_type &operator()(size_t c) const
	{
		assert(c < m_size);
		int p = posInConnections[c];
		if(p != -1)
			return con[p].value();
		else
			assert(0 && "const_and_not_available");
	}
	value_type &operator()(size_t c)
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
	bool has_connection(size_t c) const
	{
		assert(c < m_size);
		return posInConnections[c] != -1;
	}
};


}
#endif /* UNSORTED_VECTOR_H_ */
