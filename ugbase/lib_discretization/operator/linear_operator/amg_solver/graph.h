/*
 *  graph.h
 *  flexamg
 *
 *  Created by Martin Rupp on 25.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#pragma once

#include <fstream>

namespace ug{
//!
//! cgraph graph class
//! with fixed max nr of m_connections.
class cgraph
{
public:	
	//!
	//! constructor
	//! @param n nr of nodes
	//! @param nnz max. nr of m_connections
	cgraph() : m_rowStart(0), m_connections(0), m_iTotalNrOfConnections(0)
	{
		m_size = 0;
		m_iTotalNrOfConnections = 0;
	}
	
	cgraph(size_t n) : m_rowStart(0), m_connections(0), m_iTotalNrOfConnections(0)
	{
		m_size = n;
		
		m_iTotalNrOfConnections = n;

		m_connections.reserve(m_iTotalNrOfConnections);
		m_connections.resize(0);

		m_rowStart.resize(m_size+1);
		for(size_t i=0; i<m_size+1; i++) m_rowStart[i] = -1;
		m_iNrOfConnections.resize(m_size+1);

		FORCE_CREATION { print(); pr(0);} 
	}
	
	//!
	//! destructor
	~cgraph()
	{

	}
	
	//!
	//! @param i init thi
	void init(size_t i)
	{
		UG_ASSERT(i>=0 && i<m_size, "");
		UG_ASSERT(m_rowStart[i] == -1, "already inited");
		if(i == 0)
			m_rowStart[i] = 0;
		else 
		{
			UG_ASSERT(m_rowStart[i-1] != -1, "prev not inited!");
			m_rowStart[i] = m_rowStart[i-1]+m_iNrOfConnections[i-1];
		}
	}
	
	size_t getNrOfConnections(size_t row) const
	{
		return m_iNrOfConnections[row];
	}

	size_t getConnection(size_t row, size_t conn_nr) const
	{
		UG_ASSERT(conn_nr >= 0 && conn_nr < getNrOfConnections(row), "");
		return m_connections[m_rowStart[row] + conn_nr];
	}

	void setConnection(size_t from, size_t to)
	{
		UG_ASSERT(m_rowStart[from] != -1, "forgot an init?");
		UG_ASSERT(m_rowStart[from+1] == -1, "only from back to front!");

		m_connections.resize(m_iTotalNrOfConnections++);
		m_connections[m_rowStart[from] + m_iNrOfConnections[from]] = to;
		m_iNrOfConnections[from]++;
	}

	class cRowIterator
	{
	public:
		const cgraph &C;
		size_t m_position;
		size_t row;
	public:
		inline cRowIterator(const cgraph &C_, size_t row_) : C(C_)
		{
			row = row_;
			rewind();
		}
		inline cRowIterator(const cRowIterator &other) : C(other.C) { row = other.row; m_position = other.m_position; }

		inline size_t operator () ()const {return C.getConnection(row,m_position);}

		inline void operator ++() {	m_position++; }

		inline void rewind() { m_position = 0;}
		inline size_t getPos() const{	return m_position;}

		inline bool isEnd() const { return m_position >=C.getNrOfConnections(row); }
	};

	cRowIterator beginRow(size_t row) const
	{
		return cRowIterator(*this, row);
	}

	void pr(size_t i)
	{
		cout << "graph row " << i << ":" << endl;
		for(cRowIterator it = beginRow(i); !it.isEnd(); ++it)
			cout << it() << " ";
		cout << endl;
		cout.flush();
	}
	
	void print()
	{
		cout << "============= graph ================ " << endl;
		for(size_t i=0; i < m_size; i++)
		{
			cout << i << ": ";
			for(cRowIterator it = beginRow(i); !it.isEnd(); ++it)
				cout << it() << " ";
			cout << endl;
		}
	}


	void writeToFile(const char *filename, const cAMG_helper &h, int level) const
	{
		fstream file(filename, ios::out);
		file << /*CONNECTION_VIEWER_VERSION*/ 1 << endl;
		file << 2 << endl;
		
		h.writePosToStream(file);	
		file << 1 << endl;
		for(size_t i=0; i < m_size; i++)
		{
			for(cRowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
				file << h.GetOriginalIndex(level, i) << " " << h.GetOriginalIndex(level, conn()) << " " << 1 << endl;
		}
	}	
	
	void transpose()
	{
		cgraph G;
		G.createAsTransposeOf(*this);
		swap(m_iNrOfConnections, G.m_iNrOfConnections);
		swap(m_iTotalNrOfConnections, G.m_iTotalNrOfConnections);
		swap(m_connections, G.m_connections);
		swap(m_rowStart, G.m_rowStart);
		swap(m_size, G.m_size);
	}
	
	void createAsTransposeOf(const cgraph &other)
	{
		//ASSERT(conn == NULL && consmem == NULL);
		m_size = other.m_size;
		m_iNrOfConnections.resize(m_size+1);
		m_rowStart.resize(m_size+1);

		for(size_t i=0; i<m_size+1; i++) m_iNrOfConnections[i] = 0;

		m_iTotalNrOfConnections = 0;
		for(size_t i=0; i<m_size; i++)
		{
			for(cRowIterator it = other.beginRow(i); !it.isEnd(); ++it)
				m_iNrOfConnections[it()]++;
			m_iTotalNrOfConnections += other.getNrOfConnections(i);
		}
		
		m_connections.resize(m_iTotalNrOfConnections);
		m_rowStart[0] = 0;
		for(size_t i=1; i < m_size+1; i++)
		{
			m_rowStart[i] = m_rowStart[i-1] + m_iNrOfConnections[i-1];
			m_iNrOfConnections[i-1] = 0;
		}
		UG_ASSERT((size_t)m_rowStart[m_size] <= m_iTotalNrOfConnections, "");
		
		for(size_t i=0; i < m_size; i++)
			for(cRowIterator it = other.beginRow(i); !it.isEnd(); ++it)
			{
				size_t from = it();
				size_t to = i;
				m_connections[m_rowStart[from]+m_iNrOfConnections[from]] = to;
				m_iNrOfConnections[from]++;
			}
	}
	
		
	
	size_t size() { return m_size; }
	
private:
	size_t m_size;

	std::vector<int> m_rowStart;
	std::vector<size_t> m_connections;

	std::vector<size_t> m_iNrOfConnections;
	size_t m_iTotalNrOfConnections;
};

} // namespace ug
