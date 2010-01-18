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
 * the adjacency-map-structure holds the number of adjacent elements
 * for each node. The adjacency-map holds indices to the adjacent
 * elements for a node. The indices in the adjaceny-map are ordered so
 * that first all adjacent nodes of the first node appear, then the
 * adjacent nodes of the second node and so forth.
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
						TGeomBaseObj* pGeomObjsOut = NULL,
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
	adjacencyMapStructureOut.resize(grid.num<Elem>());
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
	
//	clean up
	if(!paIndex)
		grid.detach_from<Elem>(aIndex);
}						

}//	end of namespace

#endif
