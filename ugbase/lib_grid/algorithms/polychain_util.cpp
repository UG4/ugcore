/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#include "polychain_util.h"
#include "lib_grid/callbacks/callbacks.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
std::pair<Vertex*, Edge*>
GetNextSectionOfPolyChain(Grid& grid, std::pair<Vertex*, Edge*> lastSection,
						  Grid::edge_traits::callback cbEdgeIsInPolyChain)
{
	if(!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we have to enable this option, since nothing works without it in reasonable time.
		LOG("WARNING in GetFirstVertexOfPolyChain(...): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES.\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_EDGES);
	}

//	get the vertex which is connected to the vertex in lastSection->frist
//	the edge in lastSection->second
	Vertex* nVrt = GetConnectedVertex(lastSection.second, lastSection.first);
	
//	find the next edge
	Grid::AssociatedEdgeIterator edgesEnd = grid.associated_edges_end(nVrt);
	for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(nVrt);
		iter != edgesEnd; ++iter)
	{
		if(cbEdgeIsInPolyChain(*iter) && (*iter != lastSection.second)){
		//	we got the next edge
			return make_pair(nVrt, *iter);
		}
	}
//	we couldn't find another section. Return an empty one
	return std::make_pair<Vertex*, Edge*>(nullptr, nullptr);
}

////////////////////////////////////////////////////////////////////////
bool SplitIrregularPolyChain(SubsetHandler& sh, int srcIndex, int targetIndex)
{
	if(!sh.grid())
		throw(UGError("No grid assigned to subset handler"));
	
	Grid& grid = *sh.grid();
	
	if(!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we have to enable this option, since nothing works without it in reasonable time.
		LOG("WARNING in SplitPolyChain(...): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES.\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_EDGES);
	}
	
//	we'll start at the first section of the polychain
	pair<Vertex*, Edge*> curSec = GetFirstSectionOfPolyChain(grid,
												sh.begin<Edge>(srcIndex),
												sh.end<Edge>(srcIndex),
												IsInSubset(sh, srcIndex));

//	Follow the polychain until either the end is reached or an irregular
//	vertex is encountered. Those conditions are enough.


	grid.begin_marking();
	
	size_t numEdgesEncountered = 0;
	
	while(curSec.first)
	{
		Vertex* curVrt = curSec.first;
		Edge* curEdge = curSec.second;
		grid.mark(curVrt);
		grid.mark(curEdge);
		
		numEdgesEncountered++;
		
	//	check whether the connected vertex is marked or irregular
		Vertex* cv = GetConnectedVertex(curEdge, curVrt);
	
		if(grid.is_marked(cv))
			break;

		size_t counter = 0;
		
		Grid::AssociatedEdgeIterator edgesEnd = grid.associated_edges_end(cv);
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(cv);
			iter != edgesEnd; ++iter)
		{
			if(sh.get_subset_index(*iter) == srcIndex)
				++counter;
		}
		
		if(counter > 2){
		//	the vertex is irregular
			break;
		}
		
	//	ok - everything's still regular
		curSec = GetNextSectionOfPolyChain(grid, curSec, IsInSubset(sh, srcIndex));
	}
	
//	if all edges have been encountered, we're done.
	if(numEdgesEncountered == sh.num<Edge>(srcIndex)){
		grid.end_marking();
		return false;
	}
	
//	some edges are left. assign them to the target subset
	EdgeIterator iter = sh.begin<Edge>(srcIndex);
	while(iter != sh.end<Edge>(srcIndex))
	{
		Edge* e = *iter;
		++iter;
		if(!grid.is_marked(e))
			sh.assign_subset(e, targetIndex);
	}
	
	grid.end_marking();
	
	return true;
}

}//	end of namespace
