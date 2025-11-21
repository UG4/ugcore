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

#ifndef __H__LIB_GRID__DUAL_GRAPH__
#define __H__LIB_GRID__DUAL_GRAPH__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
/**
 * The created dual graph is described in the following form.
 * the adjacency-map-structure holds for each node an index to where
 * associated nodes are written in the adjacency-map.
 * The adjacency-map holds indices to the adjacent
 * elements for a node. The indices in the adjacency-map are ordered so
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
 * the graph. If pgoc == nullptr (by default), the whole grid is used.
 *
 * If you pass a multi-grid to this method, elements will be indexed for level 0
 * first, then for level 1 and so on.
 */
template <typename TGeomBaseObj, typename TIndexType>
void ConstructDualGraph(std::vector<TIndexType>& adjacencyMapStructureOut,
						std::vector<TIndexType>& adjacencyMapOut,
						Grid& grid, Attachment<TIndexType>* paIndex = nullptr,
						TGeomBaseObj** pGeomObjsOut = nullptr,
						NeighborhoodType nbhType = NHT_DEFAULT,
						const GridObjectCollection* pgoc = nullptr)
{
	using namespace std;
	using Elem = TGeomBaseObj;
	using ElemIterator = typename geometry_traits<Elem>::iterator;
	
//	set up index attachment and attachment accessor
	using AIndex = Attachment<TIndexType>;
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
 * pEdgeWgtsOut may be specified with nullptr. If it is given, then
 * the weight for each entry in adacencyMapOut will be written to
 * pEdgeWgtsOut. vWeight for each vertical connection (parents / children)
 * and hWeigth for each horizontal connection (neighbors).
 *
 * TGeomBaseObj can be either Vertex, Edge, Face or Volume
 *
 * Elements will be indexed for level 0 first, then for level 1 and so on.
 */
template <typename TGeomBaseObj, typename TIndexType>
void ConstructDualGraphMG(std::vector<TIndexType>& adjacencyMapStructureOut,
						std::vector<TIndexType>& adjacencyMapOut,
						std::vector<TIndexType>* pEdgeWeightsOut,
						MultiGrid& mg, size_t baseLevel = 0,
						int hWeight = 1, int vWeight = 1,
						Attachment<TIndexType>* paIndex = nullptr,
						TGeomBaseObj** pGeomObjsOut = nullptr,
						NeighborhoodType nbhType = NHT_DEFAULT)
{
	using namespace std;
	using Elem = TGeomBaseObj;
	using ElemIterator = typename geometry_traits<Elem>::iterator;

//	set up index attachment and attachment accessor
	using AIndex = Attachment<TIndexType>;
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


template <typename TBaseElem>
class DualGraphNeighborCollector
{
	public:
		virtual ~DualGraphNeighborCollector() = default;
		virtual void collect_neighbors(std::vector<TBaseElem*>& neighborsOut, TBaseElem* elem) = 0;
};

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
template <typename TGeomBaseObj, typename TIndexType>
void ConstructDualGraphMGLevel(
		std::vector<TIndexType>& adjacencyMapStructureOut,
		std::vector<TIndexType>& adjacencyMapOut,
		MultiGrid& mg, size_t level,
		Attachment<TIndexType>* paIndex = nullptr,
		TGeomBaseObj** pGeomObjsOut = nullptr,
		NeighborhoodType nbhType = NHT_DEFAULT,
		DualGraphNeighborCollector<TGeomBaseObj>* neighborCollector = nullptr)
{
	using namespace std;
	using Elem = TGeomBaseObj;
	using ElemIterator = typename geometry_traits<Elem>::iterator;

//	set up index attachment and attachment accessor
	using AIndex = Attachment<TIndexType>;
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
		if (neighborCollector)
			neighborCollector->collect_neighbors(vNeighbours, elem);
		else
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
