// created by Sebastian Reiter
// y10 m12 d13
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__POLYCHAIN_UTIL_IMPL__
#define __H__LIB_GRID__POLYCHAIN_UTIL_IMPL__

#include "geom_obj_util/vertex_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TEdgeIterator>
size_t
GetPolyChainType(Grid& grid, TEdgeIterator edgesBegin,
				  TEdgeIterator edgesEnd,
				  Grid::edge_traits::callback cbEdgeIsInPolyChain)
{
//	if the chain is empty, there's nothing to do
	if(edgesBegin == edgesEnd)
		return PCT_EMPTY;

//	check for each vertex to how many chain-edges it is connected
	size_t numBnd = 0;
	bool irregular = false;
	
	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
	{
		for(size_t i = 0; i < 2; ++i){
			Vertex* v = (*iter)->vertex(i);
			
			size_t counter = 0;
			for(Grid::AssociatedEdgeIterator aiter = grid.associated_edges_begin(v);
				aiter != grid.associated_edges_end(v); ++aiter)
			{
				if(cbEdgeIsInPolyChain(*aiter))
					++counter;
			}
			
			if(counter == 0){
				throw(UGError("cbEdgeIsInPolyChain does not evaluate to true for "
							  "edge between edgesBegin and edgesEnd."));
			}
							  
			if(counter == 1)
				++numBnd;
			else if(counter > 2)
				irregular = true;
		}
	}
	
//	prepare the return value
	size_t type = PCT_UNKNOWN;
	
	if(numBnd == 0)
		type = PCT_CLOSED;
	else if(numBnd < 3)
		type = PCT_OPEN;
	else
		type = PCT_SEPARATED;

	if(irregular)
		type |= PCT_IRREGULAR;
		
	return type;
}
				  
////////////////////////////////////////////////////////////////////////
template <class TEdgeIterator>
std::pair<Vertex*, Edge*>
GetFirstSectionOfPolyChain(Grid& grid, TEdgeIterator edgesBegin,
							TEdgeIterator edgesEnd,
							Grid::edge_traits::callback cbEdgeIsInPolyChain)
{
	if(edgesBegin == edgesEnd)
		return std::make_pair<Vertex*, Edge*>(NULL, NULL);
	
//	since we want to prefer vertices with local edge index 0, we'll first iterate
//	over the local indices of the edges
	for(size_t locInd = 0; locInd < 2; ++locInd){
	//	iterate through all edges.
		for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
		{
			Edge* curEdge = *iter;
			Vertex* curVrt = curEdge->vertex(locInd);
			
		//	if curVrt is a boundary vertex of the given chain, then we're done.
			if(IsBoundaryVertex1D(grid, curVrt, cbEdgeIsInPolyChain)){
				return std::make_pair(curVrt, curEdge);
			}
		}
	}
	
	return std::make_pair((*edgesBegin)->vertex(0), *edgesBegin);
}

template <class TEdgeIter>
bool CreatePolyChain(std::vector<Vertex*>& polyChainOut, Grid& grid,
					TEdgeIter edgesBegin, TEdgeIter edgesEnd)
{
	polyChainOut.clear();

	grid.begin_marking();

//	mark and count edges
	int numEdges = 0;
	for(TEdgeIter iter = edgesBegin; iter != edgesEnd; ++iter){
		grid.mark(*iter);
		++numEdges;
	}

//TODO: handle open chains.
	Edge* actEdge = *edgesBegin;
	Vertex* actVrt = actEdge->vertex(1);
	polyChainOut.push_back(actEdge->vertex(0));
	grid.mark(actEdge->vertex(0));
	polyChainOut.push_back(actEdge->vertex(1));
	grid.mark(actEdge->vertex(1));

	bool bRunning = true;
	while(bRunning){
		bRunning = false;
	//	find a connected  unmarked vertex
		Grid::AssociatedEdgeIterator assEdgesEnd = grid.associated_edges_end(actVrt);
		for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(actVrt);
			eIter != assEdgesEnd; ++eIter)
		{
			Edge* e = *eIter;
			if(grid.is_marked(e)){
			//	check whether the connected vertex is unmarked
				Vertex* cv = GetConnectedVertex(e, actVrt);
				if(!grid.is_marked(cv)){
				//	push it to the chain and mark it.
					grid.mark(cv);
					polyChainOut.push_back(cv);
				//	we're done with actVrt. cv will be the new actVrt
					actVrt = cv;
					bRunning = true;
					break;
				}
			}
		}
	}

	grid.end_marking();
/*
	if(polyChainOut.size() != numEdges)
		return false;
*/
	return true;
}

}//	end of namespace

#endif
