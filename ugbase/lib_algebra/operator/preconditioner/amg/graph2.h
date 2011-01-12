/**
 * \file amg/graph.h
 *
 * \author Martin Rupp
 *
 * \date 25.11.09
 *
 * \brief a simple graph class
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__graph_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__graph_H__

#include <fstream>
#include <algorithm> // for lower_bound
#include <vector>

#include "common/assert.h"
#include "common/log.h"

#include "lib_algebra/common/stl_debug.h"


namespace ug{
//!
//! cgraph graph class
class cgraph
{
public:
	typedef const size_t * cRowIterator;
	typedef size_t * rowIterator;
public:
	/** constructor
	 */
	cgraph() : cons(0), consmem(0), iTotalNrOfConnections(0), iMaxTotalNrOfConnections(0)
	{
#ifdef FORCE_CREATION
		FORCE_CREATION { print();}
#endif
	}
	
	cgraph(size_t n) : cons(0), consmem(0), iTotalNrOfConnections(0), iMaxTotalNrOfConnections(0)
	{
		resize(n);
	}
	
	void resize(size_t n)
	{
		consmem.resize(10*n, 0);
		iMaxTotalNrOfConnections = consmem.size();
		cons.resize(n+1, NULL);
		iTotalNrOfConnections = 0;
	}

	//!
	//! destructor
	~cgraph()
	{
	}

	//! returns nr of nodes the node "node" is connected to.
	size_t num_connections(size_t node) const
	{
		size_check(node);
		if(cons[node+1] == NULL) return 0;
		return cons[node+1]-cons[node];
	}

	bool is_isolated(size_t i)
	{
		size_check(i);
		if(cons[i+1] == NULL) return false;
		return num_connections(i)==0 ||
				(num_connections(i)==1 && cons[i][0] == i);
	}

	//! returns true if graph has connection from "from" to "to", otherwise false
	/*bool has_connection(size_t from, size_t to)
	{
		size_check(from, to);

	}*/

	void init(size_t i)
	{
		int j;
		for(j=i; j>0; j--)
		{
			if(cons[j] != NULL)
				break;
		}

		if(j == 0)
		{
			cons[0] = &consmem[0];
			cons[1] = &consmem[0];
			j = 1;
		}

		while(j<=i && cons[j+1] == NULL)
		{
			cons[j+1] = cons[j];
			j++;
		}
	}
	//! set a connection from "from" to "to" if not already there
	void set_connection(size_t from, size_t to)
	{
		size_check(from, to);

		UG_ASSERT(from == size()-1 || cons[from+2] == NULL, "only from back to front! ( from is " << from
				<< ", cons[from+2] = " << cons[from+2] << "\n");

		if(cons[from+1]==NULL)
			init(from);
		UG_ASSERT(cons[from+1]!=NULL, "??? (from = " << from << ", size = " << size());

		iTotalNrOfConnections++;
		UG_ASSERT(iTotalNrOfConnections < iMaxTotalNrOfConnections, "too many connections, increase nnz!");
		cons[from+1][0] = to;
		cons[from+1]++;
	}


	rowIterator begin_row(size_t row)
	{
		size_check(row);
		if(cons[row+1] == NULL) return NULL;
		return cons[row];
	}

	rowIterator end_row(size_t row)
	{
		size_check(row);
		return cons[row+1];
	}

	cRowIterator begin_row(size_t row) const
	{
		size_check(row);
		if(cons[row+1] == NULL) return NULL;
		return cons[row];
	}

	cRowIterator end_row(size_t row) const
	{
		size_check(row);
		return cons[row+1];
	}

	//! tranpose this graph (by using create_as_tranpose of)
	void transpose()
	{
		cgraph G;
		G.create_as_transpose_of(*this);
		swap(cons, G.cons);
		swap(consmem, G.consmem);
		size_t t;
		t = iTotalNrOfConnections, iTotalNrOfConnections = G.iTotalNrOfConnections; G.iTotalNrOfConnections = t;
		t = iMaxTotalNrOfConnections, iMaxTotalNrOfConnections = G.iMaxTotalNrOfConnections; G.iMaxTotalNrOfConnections = t;
	}
	
	void symmetricize()
	{
		cgraph G;
		G.create_as_symmetricized(*this);
		swap(cons, G.cons);
		swap(consmem, G.consmem);
		size_t t;
		t = iTotalNrOfConnections, iTotalNrOfConnections = G.iTotalNrOfConnections; G.iTotalNrOfConnections = t;
		t = iMaxTotalNrOfConnections, iMaxTotalNrOfConnections = G.iMaxTotalNrOfConnections; G.iMaxTotalNrOfConnections = t;
	}

	void create_as_symmetricized(const cgraph &other)
	{
		stdvector<size_t> rowSize(other.size());
		for(size_t i=0; i<other.size(); i++) rowSize[i] = 0;

		for(size_t i=0; i<other.size(); i++)
		{
			for(cRowIterator it = other.begin_row(i); it != other.end_row(i); ++it)
				rowSize[(*it)]++;
			iTotalNrOfConnections += other.num_connections(i);
			rowSize[i] += other.num_connections(i);
		}
		consmem.resize(iTotalNrOfConnections);
		size_t *p = &consmem[0];
		cons.resize(other.size()+1);
		for(size_t i=0; i<other.size(); i++)
		{
			cons[i] = p;
			p += rowSize[i];
			rowSize[i] = 0;
		}
		cons[other.size()] = p;

		for(size_t i=0; i < other.size(); i++)
			for(cRowIterator it = other.begin_row(i); it != other.end_row(i); ++it)
			{
				size_t from = (*it);
				size_t to = i;
				cons[from][rowSize[from]++] = to;
				cons[to][rowSize[to]++] = from;
				UG_ASSERT(cons[from]+rowSize[from] <= cons[from+1], "");
			}
		iMaxTotalNrOfConnections = iTotalNrOfConnections;
	}

	//! creates this graph as the transpose of other
	void create_as_transpose_of(const cgraph &other)
	{
		stdvector<size_t> rowSize(other.size());
		for(size_t i=0; i<other.size(); i++) rowSize[i] = 0;

		for(size_t i=0; i<other.size(); i++)
		{
			for(cRowIterator it = other.begin_row(i); it != other.end_row(i); ++it)
				rowSize[(*it)]++;
			iTotalNrOfConnections += other.num_connections(i);
		}
		consmem.resize(iTotalNrOfConnections);
		size_t *p = &consmem[0];
		cons.resize(other.size()+1);
		for(size_t i=0; i<other.size(); i++)
		{
			cons[i] = p;
			p += rowSize[i];
			rowSize[i] = 0;
		}
		cons[other.size()] = p;

		for(size_t i=0; i < other.size(); i++)
			for(cRowIterator it = other.begin_row(i); it != other.end_row(i); ++it)
			{
				size_t from = (*it);
				size_t to = i;
				cons[from][rowSize[from]++] = to;
				UG_ASSERT(cons[from]+rowSize[from] <= cons[from+1], "");
			}
		iMaxTotalNrOfConnections = iTotalNrOfConnections;
	}
		
	
	size_t size() const { return cons.size()-1; }
	
public:
	inline void size_check(size_t i) const
	{
		UG_ASSERT(i < size(), "graph contains " << size() << " nodes, but trying to access node " << i);
	}
	inline void size_check(size_t i, size_t j) const
	{
		UG_ASSERT(i < size() && j<size(),
				"graph contains " << size() << " nodes, but trying to access nodes " << i << " and " << j);
	}

	void print()
	{
		cout << "============= graph ================ " << endl;
		for(size_t i=0; i < size(); i++)
		{
			cout << i << ": ";
			for(rowIterator it=begin_row(i); it != end_row(i); ++it)
				cout << (*it) << " ";
			cout << endl;
		}
	}


protected:
	stdvector<size_t *> cons;
	stdvector<size_t> consmem;
	size_t iTotalNrOfConnections;
	size_t iMaxTotalNrOfConnections;
};


} // namespace ug

#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__GRAPH_H__
