//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d14

#ifndef __H__LIB_GRID__DUAL_GRAPH__
#define __H__LIB_GRID__DUAL_GRAPH__

#include <vector>
#include "lib_grid/lg_base.h"
#include "neighbourhood.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
/**
 * The created dual graph ist described in the following form.
 * the adjacency-map-structure holds for each node an index to where
 * associated nodes are written in the adjacency-map.
 * The adjacency-map holds indices to the adjacent
 * elements for a node. The indices in the adjaceny-map are ordered so
 * that first all adjacent nodes of the first node appear, then the
 * adjacent nodes of the second node and so forth.
 *
 * Let n be the number of nodes in the graph, m be the number of edges.
 * Then adjacencyMapStructureOut will contain n+1 entries (the last one
 * is a conveniance entry), and adjacencyMapOut will contain 2m entries
 * (edges are always bidirectional).
 *
 * please note that you'll find all entries belonging to node i in the entries
 * adjacencyMapStructureOut[i] to (but not including)
 * adjacencyMapStructureOut[i+1] of the adjacencyMapOut. 
 *
 * TIndexType has to be compatible with int.
 *
 * If you are interested in a list of the TGeomBaseObj in the order
 * corresponding to adjacencyMapStructureOut, you may pass a pointer to
 * an array of TGeomBaseObj-pointers through pGeomObjsOut. Make sure, that
 * this array has the size of grid.num<TGeomBaseObj>().
 *
 * If you are interested in the indices assigned to each face, then
 * pass a pointer to an TIndexType attachment to paIndex.
 *
 * TGeomBaseObj can be either VertexBase, EdgeBase, Face or Volume
 */
template <class TGeomBaseObj, class TIndexType>
void ConstructDualGraph(std::vector<TIndexType>& adjacencyMapStructureOut,
						std::vector<TIndexType>& adjacencyMapOut,
						Grid& grid, Attachment<TIndexType>* paIndex = NULL,
						TGeomBaseObj** pGeomObjsOut = NULL,
						NeighbourhoodType nbhType = NHT_DEFAULT)
{
	using namespace std;
	typedef TGeomBaseObj Elem;
	typedef typename geometry_traits<Elem>::iterator ElemIterator;
	
//	set up index attachment and attachment accessor
	typedef Attachment<TIndexType> AIndex;
	AIndex aIndex;
	if(paIndex)
		aIndex = *paIndex;
	
	if(!grid.has_attachment<Elem>(aIndex))
		grid.attach_to<Elem>(aIndex);
		
	Grid::AttachmentAccessor<Elem, AIndex> aaInd(grid, aIndex);
	
//	init the indices
	{
		TIndexType ind = 0;
		for(ElemIterator iter = grid.begin<Elem>(); iter != grid.end<Elem>(); ++iter, ++ind)
			aaInd[*iter] = ind;
	}
	
//	init the adjacencyMapStructure
	adjacencyMapStructureOut.resize(grid.num<Elem>() + 1);
	adjacencyMapOut.clear();

//	construct the graph
	{
		vector<Elem*> vNeighbours;
		int ind = 0;

	//	iterate through all elements
		for(ElemIterator iter = grid.begin<Elem>(); iter != grid.end<Elem>(); ++iter, ++ind)
		{
		//	get all neighbours
			CollectNeighbours(vNeighbours, *iter, grid, nbhType);
			
		//	store first entry at which the connections will be written to the map
			adjacencyMapStructureOut[ind] = adjacencyMapOut.size();
			
		//	iterate over the neighbours and push adjacent indices to the adjacencyMap
			for(size_t i = 0; i < vNeighbours.size(); ++i)
				adjacencyMapOut.push_back(aaInd[vNeighbours[i]]);
		}
	}
	
//	add the final element
	adjacencyMapStructureOut[adjacencyMapStructureOut.size() - 1] = adjacencyMapOut.size();
	
//	fill pGeomObjsOut
	if(pGeomObjsOut)
	{
		int ind = 0;
		for(ElemIterator iter = grid.begin<Elem>(); iter != grid.end<Elem>(); ++iter, ++ind)
			pGeomObjsOut[ind] = *iter;
	}
	
//	clean up
	if(!paIndex)
		grid.detach_from<Elem>(aIndex);
}						

}//	end of namespace

#endif
