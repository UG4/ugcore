/**
 * \file ugbase/lib_algebra/common/graph/old_graph.h
 *
 * \author Martin Rupp
 *
 * \date 25.11.09
 *
 * \brief a simple graph class
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG__LIB_DISC__AMG_SOLVER__GRAPH_H__
#define __H__UG__LIB_DISC__AMG_SOLVER__GRAPH_H__

#include <fstream>
#include <algorithm> // for lower_bound
#include <vector>
#include "lib_algebra/common/stl_debug.h"

#include "common/assert.h"
#include "common/log.h"



namespace ug{
//!
//! cgraph graph class
class cgraph
{
public:
	typedef stdvector<size_t>::const_iterator const_row_iterator;
	typedef stdvector<size_t>::iterator row_iterator;
public:
	/** constructor
	 */
	cgraph()
	{
	}
	
	cgraph(size_t n) : m_data(n)
	{
		for(size_t i=0; i<m_data.size(); ++i)
			m_data[i].resize(0);

#ifdef FORCE_CREATION
		FORCE_CREATION { print(); pr(0);} 
#endif
	}
	
	void resize(size_t n)
	{
		m_data.resize(n);
		for(size_t i=0; i<m_data.size(); ++i)
			m_data[i].resize(0);
	}

	//!
	//! destructor
	~cgraph()
	{
		// destructors of std::vector are getting called
	}

	//! returns nr of nodes the node "node" is connected to.
	size_t num_connections(size_t node) const
	{
		size_check(node);
		return m_data[node].size() ;
	}

	bool is_isolated(size_t i) const
	{
		size_check(i);
		return num_connections(i)==0 ||
				(num_connections(i)==1 && m_data[i][0] == i);
	}

	//! returns true if graph has connection from "from" to "to", otherwise false
	bool has_connection(size_t from, size_t to) const
	{
		size_check(from, to);
		return binary_search(begin_row(from), end_row(from), to);
	}

	//! set a connection from "from" to "to" if not already there
	void set_connection(size_t from, size_t to)
	{
		size_check(from, to);
		stdvector<size_t>::iterator it = std::lower_bound(begin_row(from), end_row(from), to);
		if(it == m_data[from].end())
			m_data[from].push_back(to);
		else if((*it) != to)
			m_data[from].insert(it, to);
	}

	row_iterator begin_row(size_t row)
	{
		size_check(row);
		return m_data[row].begin();
	}

	row_iterator end_row(size_t row)
	{
		size_check(row);
		return m_data[row].end();
	}

	const_row_iterator begin_row(size_t row) const
	{
		size_check(row);
		return m_data[row].begin();
	}

	const_row_iterator end_row(size_t row) const
	{
		size_check(row);
		return m_data[row].end();
	}

	//! tranpose this graph (by using create_as_tranpose of)
	void transpose()
	{
		cgraph G;
		G.set_as_transpose_of(*this);
		swap(m_data, G.m_data);
	}
	
	//! creates this graph as the transpose of other
	void set_as_transpose_of(const cgraph &other)
	{
		stdvector<size_t> rowSize(other.size());
		for(size_t i=0; i<other.size(); i++) rowSize[i] = 0;

		for(size_t i=0; i<other.size(); i++)
		{
			for(const_row_iterator it = other.begin_row(i); it != other.end_row(i); ++it)
				rowSize[(*it)]++;
		}
		
		m_data.resize(other.size());
		for(size_t i=0; i<other.size(); i++)
		{
			m_data[i].clear();
			m_data[i].reserve(rowSize[i]);
		}

		for(size_t i=0; i < other.size(); i++)
			for(const_row_iterator it = other.begin_row(i); it != other.end_row(i); ++it)
			{
				size_t from = (*it);
				size_t to = i;
				m_data[from].push_back(to);
			}
	}
		
	
	size_t size() const { return m_data.size(); }
	
public:
	// print functions

	//! print row i
	void pr(size_t i) const
	{
		std::cout << "graph row " << i << ", length " << num_connections(i) << ":" << std::endl;
		for(const_row_iterator it = begin_row(i); it != end_row(i); ++it)
			std::cout << (*it) << " ";
		std::cout << std::endl;
		std::cout.flush();
	}
	//! print whole graph to cout
	void print() const
	{
		std::cout << *this << std::endl;
	}

	friend std::ostream &operator << (std::ostream &out, const cgraph &g)
	{
		std::cout << "============= graph ================ " << std::endl;
		for(size_t i=0; i<g.size(); ++i)
		{
			out << "[" << i << "]:  ";
			for(const_row_iterator it = g.begin_row(i); it != g.end_row(i); ++it)
				out << (*it) << " ";
			out << std::endl;
		}
		out.flush();
		return out;
	}

	inline void size_check(size_t i) const
	{
		UG_ASSERT(i < m_data.size(), "graph contains " << m_data.size() << " nodes, but trying to access node " << i);
	}
	inline void size_check(size_t a, size_t b) const
	{
		UG_ASSERT(a < m_data.size() || b < m_data.size(), "graph contains " << m_data.size() << " nodes, but trying to access nodes " << a << " and " << b);
	}

protected:
	stdvector<stdvector<size_t> > m_data;
};


} // namespace ug

#endif // __H__UG__LIB_DISC__AMG_SOLVER__GRAPH_H__
