/*
 * Copyright (c) 2012:  G-CSC, Goethe University Frankfurt
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

/**
 * \file ugbase/lib_algebra/common/graph/new_graph.h
 *
 * \author Martin Rupp
 *
 * \date 25.11.09
 *
 * \brief a simple graph class
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG__LIB_DISC__AMG_SOLVER__graph_H__
#define __H__UG__LIB_DISC__AMG_SOLVER__graph_H__

#include <fstream>
#include <algorithm> // for lower_bound
#include <vector>

#include "common/assert.h"
#include "common/log.h"

#include "lib_algebra/common/stl_debug.h"

// consecutive memory version
// first tests with amg (1025x1025 9pt): 1138ms graph creation with old code, 333ms with new => 3x speed improvement
// old_graph is build with std::vector<std::vector< > >, new_graph uses consecutive memory,
// drawbacks at the moment:
// - construction only from 0 to n, so when you set row 3, you may not change item 1
//		amg can work with that, but perhaps we will use some pRowStart/pRowEnd mechanism like in ug::SparseMatrix
// - fixed amout of connections
//		i will change that, memory will be copied like in std::vector/reserve
namespace ug{
//!
//! cgraph graph class
class cgraph
{
public:
	using const_row_iterator = const size_t *;
	using row_iterator = size_t *;
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
		cons.resize(n+1, nullptr);
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
		if(cons[node+1] == nullptr) return 0;
		return cons[node+1]-cons[node];
	}

	bool is_isolated(size_t i) const
	{
		size_check(i);
		if(cons[i+1] == nullptr) return false;
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
			if(cons[j] != nullptr)
				break;
		}

		if(j == 0)
		{
			cons[0] = &consmem[0];
			cons[1] = &consmem[0];
			j = 1;
		}

		while(j<=i && cons[j+1] == nullptr)
		{
			cons[j+1] = cons[j];
			j++;
		}
	}
	//! set a connection from "from" to "to" if not already there
	void set_connection(size_t from, size_t to)
	{
		size_check(from, to);

		UG_ASSERT(from == size()-1 || cons[from+2] == nullptr, "only from back to front! ( from is " << from
				<< ", cons[from+2] = " << cons[from+2] << "\n");

		if(cons[from+1]==nullptr)
			init(from);
		UG_ASSERT(cons[from+1]!=nullptr, "??? (from = " << from << ", size = " << size());

		if(iTotalNrOfConnections+1 >= iMaxTotalNrOfConnections)
			increase_maxtotalnrofconnections();
		iTotalNrOfConnections++;
		cons[from+1][0] = to;
		cons[from+1]++;
	}



	row_iterator begin_row(size_t row)
	{
		size_check(row);
		if(cons[row+1] == nullptr) return nullptr;
		return cons[row];
	}

	row_iterator end_row(size_t row)
	{
		size_check(row);
		return cons[row+1];
	}

	const_row_iterator begin_row(size_t row) const
	{
		size_check(row);
		if(cons[row+1] == nullptr) return nullptr;
		return cons[row];
	}

	const_row_iterator end_row(size_t row) const
	{
		size_check(row);
		return cons[row+1];
	}

	//! tranpose this graph (by using create_as_tranpose of)
	void transpose()
	{
		cgraph G;
		G.set_as_transpose_of(*this);
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
			for(const_row_iterator it = other.begin_row(i); it != other.end_row(i); ++it)
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
			for(const_row_iterator it = other.begin_row(i); it != other.end_row(i); ++it)
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
	void set_as_transpose_of(const cgraph &other)
	{
		stdvector<size_t> rowSize(other.size());
		for(size_t i=0; i<other.size(); i++) rowSize[i] = 0;

		for(size_t i=0; i<other.size(); i++)
		{
			for(const_row_iterator it = other.begin_row(i); it != other.end_row(i); ++it)
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
			for(const_row_iterator it = other.begin_row(i); it != other.end_row(i); ++it)
			{
				size_t from = (*it);
				size_t to = i;
				cons[from][rowSize[from]++] = to;
				UG_ASSERT(cons[from]+rowSize[from] <= cons[from+1], "");
			}
		iMaxTotalNrOfConnections = iTotalNrOfConnections;
	}
		
	
	size_t size() const { return cons.size()-1; }
	
private:
	void increase_maxtotalnrofconnections()
	{
		size_t *pOld = &consmem[0];
		consmem.resize((consmem.size()+1)*2);
		size_t *pNew = &consmem[0];
		if(pOld != pNew)
		{
			for(size_t i=0; i<consmem; i++)
			{
				if(cons[i] == nullptr) continue;
				cons[i] -= pOld;
				cons[i] += pNew;
			}
		}
	}

	inline void size_check(size_t i) const
	{
		UG_ASSERT(i < size(), "graph contains " << size() << " nodes, but trying to access node " << i);
	}
	inline void size_check(size_t i, size_t j) const
	{
		UG_ASSERT(i < size() && j<size(),
				"graph contains " << size() << " nodes, but trying to access nodes " << i << " and " << j);
	}
public:
	void print() const
	{
		cout << "============= graph ================ " << endl;
		for(size_t i=0; i < size(); i++)
		{
			cout << i << ": ";
			for(row_iterator it=begin_row(i); it != end_row(i); ++it)
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

#endif