// created by Sebastian Reiter
// y10 m12 d13
// s.b.reiter@googlemail.com

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
	return std::make_pair<Vertex*, Edge*>(NULL, NULL);	
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
