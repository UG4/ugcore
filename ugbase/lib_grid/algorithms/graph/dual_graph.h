#ifndef __H__LIB_GRID__DUAL_GRAPH__
#define __H__LIB_GRID__DUAL_GRAPH__

#include <vector>
#include "lib_grid/lg_base.h"

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
 * is a convenience entry), and adjacencyMapOut will contain 2m entries
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
 * If you are interested in the indices assigned to each element, then
 * pass a pointer to an TIndexType attachment to paIndex.
 *
 * TGeomBaseObj can be either Vertex, Edge, Face or Volume
 *
 * Through pgoc you can specify a subset of the grid which shall be put into
 * the graph. If pgoc == NULL (by default), the whole grid is used.
 *
 * If you pass a multi-grid to this method, elements will be indexed for level 0
 * first, then for level 1 and so on.
 */
template <class TGeomBaseObj, class TIndexType>
void ConstructDualGraph(std::vector<TIndexType>& adjacencyMapStructureOut,
						std::vector<TIndexType>& adjacencyMapOut,
						Grid& grid, Attachment<TIndexType>* paIndex = NULL,
						TGeomBaseObj** pGeomObjsOut = NULL,
						NeighborhoodType nbhType = NHT_DEFAULT,
						const GridObjectCollection* pgoc = NULL)
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
		grid.attach_to_dv<Elem>(aIndex, -1);
		
	Grid::AttachmentAccessor<Elem, AIndex> aaInd(grid, aIndex);
	
	GridObjectCollection goc;
	if(pgoc)
		goc = *pgoc;
	else
		goc = grid.get_grid_objects();

//	init the indices
	TIndexType ind = 0;
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
	{
		for(ElemIterator iter = goc.begin<Elem>(lvl); iter != goc.end<Elem>(lvl);
			++iter, ++ind)
		{
			aaInd[*iter] = ind;
		}
	}
	
//	init the adjacencyMapStructure
	adjacencyMapStructureOut.resize(goc.num<Elem>() + 1);
	adjacencyMapOut.clear();

//	construct the graph
	{
		vector<Elem*> vNeighbours;
		int ind = 0;

	//	iterate through all elements
		for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
		{
			for(ElemIterator iter = goc.begin<Elem>(lvl);
				iter != goc.end<Elem>(lvl); ++iter, ++ind)
			{
			//	get all neighbours
				CollectNeighbors(vNeighbours, *iter, grid, nbhType);

			//	store first entry at which the connections will be written to the map
				adjacencyMapStructureOut[ind] = adjacencyMapOut.size();

			//	iterate over the neighbours and push adjacent indices to the adjacencyMap
				for(size_t i = 0; i < vNeighbours.size(); ++i)
					adjacencyMapOut.push_back(aaInd[vNeighbours[i]]);
			}
		}
	}
	
//	add the final element
	adjacencyMapStructureOut[adjacencyMapStructureOut.size() - 1] = adjacencyMapOut.size();
	
//	fill pGeomObjsOut
	if(pGeomObjsOut)
	{
		int ind = 0;
		for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl){
			for(ElemIterator iter = goc.begin<Elem>(lvl);
				iter != goc.end<Elem>(lvl); ++iter, ++ind)
			{
				pGeomObjsOut[ind] = *iter;
			}
		}
	}
	
//	clean up
	if(!paIndex)
		grid.detach_from<Elem>(aIndex);
}						


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
 * is a convenience entry), and adjacencyMapOut will contain 2m entries
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
 * If you are interested in the indices assigned to each element, then
 * pass a pointer to an TIndexType attachment to paIndex.
 *
 * pEdgeWgtsOut may be specified with NULL. If it is given, then
 * the weight for each entry in adacencyMapOut will be written to
 * pEdgeWgtsOut. vWeight for each vertical connection (parents / children)
 * and hWeigth for each horizontal connection (neighbors).
 *
 * TGeomBaseObj can be either Vertex, Edge, Face or Volume
 *
 * Elements will be indexed for level 0 first, then for level 1 and so on.
 */
template <class TGeomBaseObj, class TIndexType>
void ConstructDualGraphMG(std::vector<TIndexType>& adjacencyMapStructureOut,
						std::vector<TIndexType>& adjacencyMapOut,
						std::vector<TIndexType>* pEdgeWeightsOut,
						MultiGrid& mg, size_t baseLevel = 0,
						int hWeight = 1, int vWeight = 1,
						Attachment<TIndexType>* paIndex = NULL,
						TGeomBaseObj** pGeomObjsOut = NULL,
						NeighborhoodType nbhType = NHT_DEFAULT)
{
	using namespace std;
	typedef TGeomBaseObj Elem;
	typedef typename geometry_traits<Elem>::iterator ElemIterator;

//	set up index attachment and attachment accessor
	typedef Attachment<TIndexType> AIndex;
	AIndex aIndex;
	if(paIndex)
		aIndex = *paIndex;

	if(!mg.has_attachment<Elem>(aIndex))
		mg.attach_to<Elem>(aIndex);

	Grid::AttachmentAccessor<Elem, AIndex> aaInd(mg, aIndex);

//	init the indices
	TIndexType ind = 0;
	size_t numElems = 0;
	for(size_t lvl = baseLevel; lvl < mg.num_levels(); ++lvl)
	{
		numElems += mg.num<Elem>(lvl);
		for(ElemIterator iter = mg.begin<Elem>(lvl); iter != mg.end<Elem>(lvl);
			++iter, ++ind)
		{
			aaInd[*iter] = ind;
		}
	}

//	init the adjacencyMapStructure
	adjacencyMapStructureOut.resize(numElems + 1);
	adjacencyMapOut.clear();
	if(pEdgeWeightsOut)
		pEdgeWeightsOut->clear();

//	construct the graph
	{
		vector<Elem*> vNeighbours;
		int ind = 0;

	//	iterate through all elements
		for(size_t lvl = baseLevel; lvl < mg.num_levels(); ++lvl)
		{
			for(ElemIterator iter = mg.begin<Elem>(lvl);
				iter != mg.end<Elem>(lvl); ++iter, ++ind)
			{
				Elem* elem = *iter;

			//	get all neighbours
				CollectNeighbors(vNeighbours, elem, mg, nbhType);

			//	store first entry at which the connections will be written to the map
				adjacencyMapStructureOut[ind] = adjacencyMapOut.size();

			//	iterate over the neighbours and push adjacent indices to the adjacencyMap
				for(size_t i = 0; i < vNeighbours.size(); ++i)
					adjacencyMapOut.push_back(aaInd[vNeighbours[i]]);

				if(pEdgeWeightsOut){
					for(size_t i = 0; i < vNeighbours.size(); ++i)
						pEdgeWeightsOut->push_back(hWeight);
				}

			//	add a connection to the parents to the list
				GridObject* parent = mg.get_parent(elem);
				if(parent && Elem::type_match(parent)){
					adjacencyMapOut.push_back(aaInd[static_cast<Elem*>(parent)]);
					if(pEdgeWeightsOut)
						pEdgeWeightsOut->push_back(vWeight);
				}

			//	add children to the list
			//todo: add edge-weights
				size_t numChildren = mg.num_children<Elem>(elem);
				for(size_t i = 0; i < numChildren; ++i){
					adjacencyMapOut.push_back(aaInd[mg.get_child<Elem>(elem, i)]);
				}

				if(pEdgeWeightsOut){
					for(size_t i = 0; i < numChildren; ++i)
						pEdgeWeightsOut->push_back(1);
				}
			}
		}
	}

//	add the final element
	adjacencyMapStructureOut[adjacencyMapStructureOut.size() - 1] = adjacencyMapOut.size();

//	fill pGeomObjsOut
	if(pGeomObjsOut)
	{
		int ind = 0;
		for(size_t lvl = baseLevel; lvl < mg.num_levels(); ++lvl){
			for(ElemIterator iter = mg.begin<Elem>(lvl);
				iter != mg.end<Elem>(lvl); ++iter, ++ind)
			{
				pGeomObjsOut[ind] = *iter;
			}
		}
	}

//	clean up
	if(!paIndex)
		mg.detach_from<Elem>(aIndex);
}


////////////////////////////////////////////////////////////////////////
/** This method creates a dual graph for the given level only.
 *
 * The i-th element of the given level will be the i-th node in the generated
 * adjacency structure.
 *
 * The created dual graph ist described in the following form.
 * the adjacency-map-structure holds for each node an index to where
 * associated nodes are written in the adjacency-map.
 * The adjacency-map holds indices to the adjacent
 * elements for a node. The indices in the adjaceny-map are ordered so
 * that first all adjacent nodes of the first node appear, then the
 * adjacent nodes of the second node and so forth.
 *
 * Note that the elements are processed one level after another.
 *
 * Let n be the number of nodes in the graph, m be the number of edges.
 * Then adjacencyMapStructureOut will contain n+1 entries (the last one
 * is a convenience entry), and adjacencyMapOut will contain 2m entries
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
 * this array has the size of grid.num<TGeomBaseObj>(level).
 *
 * If you are interested in the indices assigned to each element, then
 * pass a pointer to an TIndexType attachment to paIndex.
 *
 * TGeomBaseObj can be either Vertex, Edge, Face or Volume
 */
template <class TGeomBaseObj, class TIndexType>
void ConstructDualGraphMGLevel(
		std::vector<TIndexType>& adjacencyMapStructureOut,
		std::vector<TIndexType>& adjacencyMapOut,
		MultiGrid& mg, size_t level,
		Attachment<TIndexType>* paIndex = NULL,
		TGeomBaseObj** pGeomObjsOut = NULL,
		NeighborhoodType nbhType = NHT_DEFAULT)
{
	using namespace std;
	typedef TGeomBaseObj Elem;
	typedef typename geometry_traits<Elem>::iterator ElemIterator;

//	set up index attachment and attachment accessor
	typedef Attachment<TIndexType> AIndex;
	AIndex aIndex;
	if(paIndex)
		aIndex = *paIndex;

	if(!mg.has_attachment<Elem>(aIndex))
		mg.attach_to<Elem>(aIndex);

	Grid::AttachmentAccessor<Elem, AIndex> aaInd(mg, aIndex);

//	init the indices
	size_t numElems = mg.num<Elem>(level);
	{
		TIndexType ind = 0;
		for(ElemIterator iter = mg.begin<Elem>(level); iter != mg.end<Elem>(level);
			++iter, ++ind)
		{
			aaInd[*iter] = ind;
		}
	}

//	init the adjacencyMapStructure
	adjacencyMapStructureOut.resize(numElems + 1);
	adjacencyMapOut.clear();

//	construct the graph
	vector<Elem*> vNeighbours;
	int ind = 0;

//	generate adjacency structure first
	for(ElemIterator iter = mg.begin<Elem>(level);
		iter != mg.end<Elem>(level); ++iter, ++ind)
	{
		Elem* elem = *iter;

	//	get all neighbours
		CollectNeighbors(vNeighbours, elem, mg, nbhType);

	//	store first entry at which the connections will be written to the map
		adjacencyMapStructureOut[ind] = adjacencyMapOut.size();

	//	iterate over the neighbours and push adjacent indices to the adjacencyMap
		for(size_t i = 0; i < vNeighbours.size(); ++i)
			adjacencyMapOut.push_back(aaInd[vNeighbours[i]]);
	}

//	add the final element
	adjacencyMapStructureOut[adjacencyMapStructureOut.size() - 1] = adjacencyMapOut.size();

//	fill pGeomObjsOut
	if(pGeomObjsOut)
	{
		int ind = 0;
		for(ElemIterator iter = mg.begin<Elem>(level);
			iter != mg.end<Elem>(level); ++iter, ++ind)
		{
			pGeomObjsOut[ind] = *iter;
		}
	}

//	clean up
	if(!paIndex)
		mg.detach_from<Elem>(aIndex);
}

}//	end of namespace

#endif
