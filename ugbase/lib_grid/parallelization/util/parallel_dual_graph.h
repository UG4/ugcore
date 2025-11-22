/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__PARALLEL_DUAL_GRAPH__
#define __H__LIB_GRID__PARALLEL_DUAL_GRAPH__

#include <vector>
#include "lib_grid/lg_base.h"
#include "pcl/pcl_process_communicator.h"
#include "../distributed_grid.h"
#include "../parallelization_util.h"

namespace ug
{

///	Generates the parallel dual graph of a MultiGrid as, e.g., required by Parmetis
/** Indices on elements of TGeomBaseObj are distributed in the order in which
 * one iterates over those elements in the associated mg starting at level 0.
 * Note that only elements which were considered are associated with an index.
 * \todo	The current implementation has support for MultiGrids only.
 * 			Support for flat grids should be added.*/
template <typename TGeomBaseObj, typename TIndexType,
		  typename TConnectingObj = typename TGeomBaseObj::side>
class ParallelDualGraph{
	public:
		using element_iterator_t = typename std::vector<TGeomBaseObj*>::iterator;

		ParallelDualGraph(MultiGrid* pmg = nullptr);
		~ParallelDualGraph();

		void set_grid(MultiGrid* pmg);

	///	Access to the graph which was generated during the last call to generate_graph.
	/**	\sa	num_graph_vertices, num_graph_edges, adjacency_map_structure,
	 * 		adjacency_map, parallel_offset_map
	 * \{ */
		TIndexType num_graph_vertices();
		TIndexType num_graph_edges();
		TIndexType*	adjacency_map_structure();
		TIndexType* adjacency_map();
		TIndexType* parallel_offset_map();
	/**	\} */

	///	Some vertices are not considered for the dual graph (e.g. ghosts).
	/**	Note that this method only works while the underlying grid has not been
	 * changed since the last call to generate_graph.*/
		bool was_considered(TGeomBaseObj* o);

	///	returns the graph-vertex for the given index
		TGeomBaseObj* get_element(size_t gvrtIndex)		{return m_elems[gvrtIndex];}

	///	returns the graph-edge for the given index
		/** The i-th connection is the side_t element of the grid that represents
		 * the i-th connection stored in m_adjacencyMap.
		 */
		TConnectingObj* get_connection(size_t connIndex)		{return m_connections.at(connIndex);}

	///	returns the begin-iterator to the elements corresponding to graph vertices*/
		element_iterator_t elements_begin()				{return m_elems.begin();}

	///	returns the begin-iterator to the elements corresponding to graph vertices*/
		element_iterator_t elements_end()				{return m_elems.end();}

	///	returns the graph-vertex index of the given element
		TIndexType get_index(TGeomBaseObj* elem)		{return m_aaElemIndex[elem];}

	///	generates the graph for the given level.
	/**	Use adjacency_map_structure, adjacency_map, parallel_offset_map and
	 * num_graph_vertices to access the generated graph.
	 *
	 * \todo	add more generate_graph methods, e.g. for flat grids or complete hierarchies.*/
		void generate_graph(int level, pcl::ProcessCommunicator procCom =
											pcl::ProcessCommunicator(pcl::PCD_WORLD));

	///	returns a process communicator which only contains processes which contain an element.
	/**	\note	the parallel-offset-map is built with respect to this process-communicator.
	 * \note	this process communicator does not necessarily resemble the communicator
	 * 			specified in generate_graph.*/
		pcl::ProcessCommunicator process_communicator() const	{return m_procCom;}

	private:
		void attach_data();
		void detach_data();

	////////////////////
	//	Members
		using AElemIndex = Attachment<int>;
		using AElemIndices = Attachment<std::vector<int> >;

		pcl::ProcessCommunicator	m_procCom;
		MultiGrid*				m_pMG;
		std::vector<TGeomBaseObj*>	m_elems;
		std::vector<TConnectingObj*>	m_connections;
		std::vector<TIndexType> m_adjacencyMapStructure;
		std::vector<TIndexType> m_adjacencyMap;
		std::vector<TIndexType> m_nodeOffsetMap;
		AElemIndex				m_aElemIndex;
		AElemIndices			m_aElemIndices;
		Grid::AttachmentAccessor<TGeomBaseObj, AElemIndex>		m_aaElemIndex;
		Grid::AttachmentAccessor<TConnectingObj, AElemIndices>	m_aaElemIndices;
};

////////////////////////////////////////////////////////////////////////
/** DEPRECIATED! Use the ParallelDualGraph class instead!
 *
 * This method creates a parallel dual graph for the given level only.
 * Please refer to the parmetis documentation for a detailed overview over
 * the resulting parallel graph structure.
 *
 * The created dual graph ist described in the following form.
 * The adjacency-map-structure holds for each node an index to where
 * associated nodes are written in the adjacency-map.
 *
 * The adjacency-map holds indices to the adjacent
 * elements for a node. The indices in the adjaceny-map are ordered so
 * that first all adjacent nodes of the first node appear, then the
 * adjacent nodes of the second node and so forth. Note that those indices relate
 * to the global index of the connected node.
 *
 * The node-offset-map holds for each process the smallest global node index
 * on that process. Note that all node indices on a process are consecutive.
 * The node-offset-map holds numProcs+1 entries. The last one holds #numGlobalInds+1.
 * The number of nodes on a process is thus given by
 * (nodeOffsetMap[pInd] - nodeOffsetMap[pInd+1]).
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
 * this array has the size of grid.num<TGeomBaseObj>(level). It may however not be
 * completely filled, since ghosts are ignored in the dual graph.
 *
 * If you are interested in the indices assigned to each element, then
 * pass a pointer to an TIndexType attachment to paIndex.
 *
 * TGeomBaseObj can be either Vertex, Edge, Face or Volume and should
 * represent the elements of highest dimension in a given grid.
 */
template <typename TGeomBaseObj, typename TIndexType>
void ConstructParallelDualGraphMGLevel(
		std::vector<TIndexType>& adjacencyMapStructureOut,
		std::vector<TIndexType>& adjacencyMapOut,
		std::vector<TIndexType>& nodeOffsetMapOut,
		MultiGrid& mg, size_t level,
		pcl::ProcessCommunicator procCom,
		Attachment<TIndexType>* paIndex = nullptr,
		TGeomBaseObj** pGeomObjsOut = nullptr,
		NeighborhoodType nbhType = NHT_DEFAULT)
{
	GDIST_PROFILE_FUNC();
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

//	init the indices and count ghosts and normal elements on the fly
	size_t numElems = 0;
	{
		DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();
		for(ElemIterator iter = mg.begin<Elem>(level);
			iter != mg.end<Elem>(level); ++iter)
		{
			if(distGridMgr.is_ghost(*iter))
				aaInd[*iter] = -1;
			else
				aaInd[*iter] = numElems++;
		}
	}

//	generate the nodeOffsetMap
	nodeOffsetMapOut.clear();
	int localNodeOffset = 0;

	int numElemsTmp = (int)numElems;
	vector<int> elemCounts(procCom.size());
	procCom.allgather(&numElemsTmp, 1, PCL_DT_INT,
					  &elemCounts.front(), 1, PCL_DT_INT);

	nodeOffsetMapOut.resize(procCom.size() + 1);
	int numElemsTotal = 0;
	for(size_t i = 0; i < elemCounts.size(); ++i){
		nodeOffsetMapOut[i] = numElemsTotal;
		numElemsTotal += elemCounts[i];
	}
	nodeOffsetMapOut[elemCounts.size()] = numElemsTotal;
	localNodeOffset = nodeOffsetMapOut[procCom.get_local_proc_id()];

//	init the adjacencyMapStructure
	adjacencyMapStructureOut.resize(numElems + 1);
	adjacencyMapOut.clear();

//	construct the graph
	vector<Elem*> vNeighbours;
	int ind = 0;

//	generate adjacency structure first
	for(ElemIterator iter = mg.begin<Elem>(level);
		iter != mg.end<Elem>(level); ++iter)
	{
		Elem* elem = *iter;

		if(aaInd[elem] == -1)
			continue;

	//	get all neighbours
		CollectNeighbors(vNeighbours, elem, mg, nbhType);

	//	store first entry at which the connections will be written to the map
		adjacencyMapStructureOut[ind] = adjacencyMapOut.size();

	//	iterate over the neighbours and push adjacent indices to the adjacencyMap
		for(size_t i = 0; i < vNeighbours.size(); ++i){
			if(aaInd[vNeighbours[i]] != -1)
				adjacencyMapOut.push_back(localNodeOffset + aaInd[vNeighbours[i]]);
		}

		++ind;
	}

//	add the final element
	adjacencyMapStructureOut[adjacencyMapStructureOut.size() - 1] = adjacencyMapOut.size();

//	fill pGeomObjsOut
	if(pGeomObjsOut)
	{
		int ind = 0;
		for(ElemIterator iter = mg.begin<Elem>(level);
			iter != mg.end<Elem>(level); ++iter)
		{
			if(aaInd[*iter] != -1){
				pGeomObjsOut[ind] = *iter;
				++ind;
			}
		}
	}

//	clean up
	if(!paIndex)
		mg.detach_from<Elem>(aIndex);
}

}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "parallel_dual_graph_impl.hpp"

#endif
