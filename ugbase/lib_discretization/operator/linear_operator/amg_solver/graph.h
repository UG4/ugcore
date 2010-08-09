/**
 * \file amg_rs_prolongation.h
 *
 * \author Martin Rupp
 *
 * \date 25.11.09
 *
 * \brief a simple graph class
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__GRAPH_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__GRAPH_H__

#include <fstream>
#include <algorithm> // for lower_bound
#include <vector>

namespace ug{
//!
//! cgraph graph class
class cgraph
{
public:
	typedef std::vector<size_t>::const_iterator cRowIterator;
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
		return m_data[node].size() ;
	}

	bool is_isolated(size_t i)
	{
		return num_connections(i)==0 ||
				(num_connections(i)==1 && m_data[i][0] == i);
	}

	//! returns true if graph has connection from "from" to "to", otherwise false
	bool has_connection(size_t from, size_t to)
	{
		return binary_search(m_data[from].begin(), m_data[from].end(), to);
	}

	//! set a connection from "from" to "to" if not already there
	void set_connection(size_t from, size_t to)
	{
		std::vector<size_t>::iterator it = std::lower_bound(m_data[from].begin(), m_data[from].end(), to);
		if(it == m_data[from].end())
			m_data[from].push_back(to);
		else if((*it) != to)
			m_data[from].insert(it, to);
	}

	cRowIterator begin_row(size_t row) const
	{
		return m_data[row].begin();
	}

	cRowIterator end_row(size_t row) const
	{
		return m_data[row].end();
	}

	//! tranpose this graph (by using create_as_tranpose of)
	void transpose()
	{
		cgraph G;
		G.create_as_transpose_of(*this);
		swap(m_data, G.m_data);
	}
	
	//! creates this graph as the transpose of other
	void create_as_transpose_of(const cgraph &other)
	{
		std::vector<size_t> rowSize(other.size());
		for(size_t i=0; i<other.size(); i++) rowSize[i] = 0;

		for(size_t i=0; i<other.size(); i++)
		{
			for(cRowIterator it = other.begin_row(i); it != other.end_row(i); ++it)
				rowSize[(*it)]++;
		}
		
		m_data.resize(other.size());
		for(size_t i=0; i<other.size(); i++)
		{
			m_data[i].clear();
			m_data[i].reserve(rowSize[i]);
		}

		for(size_t i=0; i < other.size(); i++)
			for(cRowIterator it = other.begin_row(i); it != other.end_row(i); ++it)
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
	void pr(size_t i)
	{
		cout << "graph row " << i << ", length " << num_connections(i) << ":" << endl;
		for(cRowIterator it = begin_row(i); it != end_row(i); ++it)
			cout << (*it) << " ";
		cout << endl;
		cout.flush();
	}
	//! print whole graph to cout
	void print()
	{
		cout << *this << endl;
	}

	friend std::ostream &operator << (std::ostream &out, const cgraph &g)
	{
		cout << "============= graph ================ " << std::endl;
		for(size_t i=0; i<g.size(); ++i)
		{
			out << "[" << i << "]:  ";
			for(cRowIterator it = g.begin_row(i); it != g.end_row(i); ++it)
				out << (*it) << " ";
			out << std::endl;
		}
		out.flush();
		return out;
	}

protected:
	std::vector<std::vector<size_t> > m_data;
};


// could be in cpp
inline void WriteAMGGraphToFile(cgraph &G, const char *filename, const cAMG_helper &h, int level)
{
	fstream file(filename, ios::out);
	file << /*CONNECTION_VIEWER_VERSION*/ 1 << endl;

	h.writePosToStream(file);
	file << 1 << endl;
	for(size_t i=0; i < G.size(); i++)
	{
		for(cgraph::cRowIterator it = G.begin_row(i); it != G.end_row(i); ++it)
			file << h.GetOriginalIndex(level, i) << " " << h.GetOriginalIndex(level, (*it)) << "  " << endl;
	}
}

} // namespace ug

#endif __H__LIB_DISCRETIZATION__AMG_SOLVER__GRAPH_H__
