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
				  CB_ConsiderEdge cbEdgeIsInPolyChain)
{
//	check for each vertex to how many chain-edges it is connected
	size_t numBnd = 0;
	bool irregular = false;
	
	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
	{
		for(size_t i = 0; i < 2; ++i){
			VertexBase* v = (*iter)->vertex(i);
			
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
std::pair<VertexBase*, EdgeBase*>
GetFirstSectionOfPolyChain(Grid& grid, TEdgeIterator edgesBegin,
							TEdgeIterator edgesEnd,
							CB_ConsiderEdge cbEdgeIsInPolyChain)
{
	if(edgesBegin == edgesEnd)
		return std::make_pair<VertexBase*, EdgeBase*>(NULL, NULL);
	
//	since we want to prefer vertices with local edge index 0, we'll first iterate
//	over the local indices of the edges
	for(size_t locInd = 0; locInd < 2; ++locInd){
	//	iterate through all edges.
		for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
		{
			EdgeBase* curEdge = *iter;
			VertexBase* curVrt = curEdge->vertex(locInd);
			
		//	if curVrt is a boundary vertex of the given chain, then we're done.
			if(IsBoundaryVertex1D(grid, curVrt, cbEdgeIsInPolyChain)){
				return std::make_pair(curVrt, curEdge);
			}
		}
	}
	
	return std::make_pair((*edgesBegin)->vertex(0), *edgesBegin);
}

}//	end of namespace

#endif
