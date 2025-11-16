/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__parallel_dual_graph_impl__
#define __H__UG__parallel_dual_graph_impl__

#include "parallel_dual_graph.h"
#include "compol_gather_vec_attachment.h"
#include "compol_copy_attachment.h"
#include "../distributed_grid.h"

namespace ug{

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
ParallelDualGraph(MultiGrid* pmg) :
	m_procCom(pcl::PCD_EMPTY),
	m_pMG(nullptr)
{
	if(pmg)
		set_grid(pmg);
}

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
~ParallelDualGraph()
{
	if(m_pMG)
		detach_data();
}

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
void ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
set_grid(MultiGrid* pmg)
{
	if(m_pMG == pmg)
		return;

	if(m_pMG)
		detach_data();

	m_pMG = pmg;
	if(m_pMG)
		attach_data();
}

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
TIndexType ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
num_graph_vertices()
{
	if(m_adjacencyMapStructure.empty())
		return 0;

	return (TIndexType)m_adjacencyMapStructure.size() - 1;
}

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
TIndexType ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
num_graph_edges()
{
	return (TIndexType)m_adjacencyMap.size();
}

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
TIndexType* ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
adjacency_map_structure()
{
	UG_ASSERT(!m_adjacencyMapStructure.empty(),
			  "Call generate graph before calling this method!");
	return &m_adjacencyMapStructure.front();
}

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
TIndexType* ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
adjacency_map()
{
	UG_ASSERT(!m_adjacencyMap.empty(),
			  "Call generate graph before calling this method!");
	return &m_adjacencyMap.front();
}

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
TIndexType* ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
parallel_offset_map()
{
	UG_ASSERT(!m_nodeOffsetMap.empty(),
			  "Call generate graph before calling this method!");
	return &m_nodeOffsetMap.front();
}


template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
bool ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
was_considered(TGeomBaseObj* o)
{
	UG_ASSERT(m_pMG, "A MultiGrid has to be set!");
	return m_aaElemIndex[o] != -1;
}

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
void ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
attach_data()
{
	assert(m_pMG);
	m_pMG->attach_to<TGeomBaseObj>(m_aElemIndex);
	m_pMG->attach_to<TConnectingObj>(m_aElemIndices);
	m_aaElemIndex.access(*m_pMG, m_aElemIndex);
	m_aaElemIndices.access(*m_pMG, m_aElemIndices);
}

template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
void ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
detach_data()
{
	assert(m_pMG);
	m_pMG->detach_from<TGeomBaseObj>(m_aElemIndex);
	m_pMG->detach_from<TConnectingObj>(m_aElemIndices);
	m_aaElemIndex.invalidate();
	m_aaElemIndices.invalidate();
}


template <class TGeomBaseObj, class TIndexType, class TConnectingObj>
void ParallelDualGraph<TGeomBaseObj, TIndexType, TConnectingObj>::
generate_graph(int level, pcl::ProcessCommunicator procCom)
{
	GDIST_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "ParallelDualGraph-start generate_graph\n");
	UG_ASSERT(m_pMG, "A MultiGrid has to be set!");

	using namespace std;
	using Elem = TGeomBaseObj;
	using ConElem = TConnectingObj;
	using ElemIterator = typename geometry_traits<Elem>::iterator;
	using ConElemIterator = typename geometry_traits<ConElem>::iterator;

	MultiGrid& mg = *m_pMG;
	Grid::AttachmentAccessor<TGeomBaseObj, AElemIndex>& aaInd = m_aaElemIndex;
	Grid::AttachmentAccessor<TConnectingObj, AElemIndices>& aaConInds = m_aaElemIndices;

	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();
	GridLayoutMap& glm = distGridMgr.grid_layout_map();

//	init the indices and count ghosts and normal elements on the fly
	size_t numElems = 0;
	m_elems.clear();
	{
		for(ElemIterator iter = mg.begin<Elem>(level);
			iter != mg.end<Elem>(level); ++iter)
		{
			if(distGridMgr.is_ghost(*iter))
				aaInd[*iter] = -1;
			else{
				aaInd[*iter] = numElems++;
				m_elems.push_back(*iter);
			}
		}
	}

//	generate a local procCom, which only contains processes which actually contain elements.
	m_procCom = procCom.create_sub_communicator(numElems > 0);

//	generate the nodeOffsetMap and determine the localNodeOffset
	UG_DLOG(LIB_GRID, 2, "ParallelDualGraph-generate_graph: gathering element numbers\n");
	m_nodeOffsetMap.clear();
	int localNodeOffset = 0;

	if(!m_procCom.empty()){
		int numElemsTmp = (int)numElems;
		vector<int> elemCounts(m_procCom.size());
			m_procCom.allgather(&numElemsTmp, 1, PCL_DT_INT,
							  &elemCounts.front(), 1, PCL_DT_INT);

		m_nodeOffsetMap.resize(m_procCom.size() + 1);
		int numElemsTotal = 0;
		for(size_t i = 0; i < elemCounts.size(); ++i){
			m_nodeOffsetMap[i] = numElemsTotal;
			numElemsTotal += elemCounts[i];
		}
		m_nodeOffsetMap[elemCounts.size()] = numElemsTotal;
		localNodeOffset = m_nodeOffsetMap[m_procCom.get_local_proc_id()];

		UG_DLOG(LIB_GRID, 2, "ParallelDualGraph-generate_graph: gathering indices of connected elements\n");
	//	we have to gather indices of connected elements in connecting elements.
	//	This is required since we have to communicate those indices between distributed
	//	connecting elements.
		typename Grid::traits<Elem>::secure_container elems;
		for(ConElemIterator iter = mg.begin<ConElem>(level);
			iter != mg.end<ConElem>(level); ++iter)
		{
			ConElem* ce = *iter;
			aaConInds[ce].clear();
			mg.associated_elements(elems, ce);
			for(size_t i = 0; i < elems.size(); ++i){
				if(aaInd[elems[i]] != -1)
					aaConInds[ce].push_back(localNodeOffset + aaInd[elems[i]]);
			}
		}
	}

//	communicate connected elements between horizontal interfaces on this level
	using Layout = typename GridLayoutMap::Types<ConElem>::Layout::LevelLayout;
	pcl::InterfaceCommunicator<Layout> com;

	ComPol_GatherVecAttachment<Layout, AElemIndices> compolGather(mg, m_aElemIndices);
	if(glm.has_layout<ConElem>(INT_H_SLAVE)){
		com.send_data(glm.get_layout<ConElem>(INT_H_SLAVE).layout_on_level(level),
					  compolGather);
	}
	if(glm.has_layout<ConElem>(INT_H_MASTER)){
		com.receive_data(glm.get_layout<ConElem>(INT_H_MASTER).layout_on_level(level),
						 compolGather);
	}
	com.communicate();

	ComPol_CopyAttachment<Layout, AElemIndices> compolCopy(mg, m_aElemIndices);
	if(glm.has_layout<ConElem>(INT_H_MASTER)){
		com.send_data(glm.get_layout<ConElem>(INT_H_MASTER).layout_on_level(level),
					  compolCopy);
	}
	if(glm.has_layout<ConElem>(INT_H_SLAVE)){
		com.receive_data(glm.get_layout<ConElem>(INT_H_SLAVE).layout_on_level(level),
						 compolCopy);
	}
	com.communicate();

	UG_DLOG(LIB_GRID, 2, "ParallelDualGraph-generate_graph: building adjacency structure\n");
//	init the adjacencyMapStructure
	m_adjacencyMapStructure.resize(numElems + 1);
	m_adjacencyMapStructure[0] = 0;
	m_adjacencyMap.clear();
	m_connections.clear();

//	generate adjacency structure.
	typename Grid::traits<ConElem>::secure_container conElems;
	if(ConElem::dim == Elem::dim - 1){
		int ind = 0;
		for(ElemIterator iter = mg.begin<Elem>(level); iter != mg.end<Elem>(level); ++iter)
		{
			Elem* elem = *iter;
			int eInd = aaInd[elem];
			if(eInd == -1)
				continue;
			eInd += localNodeOffset;

		//	store first entry at which the connections will be written to the map
			assert(ind < (int)m_adjacencyMapStructure.size());
			m_adjacencyMapStructure[ind] = m_adjacencyMap.size();

			mg.associated_elements(conElems, elem);
		//	iterate over the neighbors and push adjacent indices to the adjacencyMap
			for(size_t i_con = 0; i_con < conElems.size(); ++i_con){
				std::vector<int>& conInds = aaConInds[conElems[i_con]];
				for(size_t i = 0; i < conInds.size(); ++i){
					UG_ASSERT(conInds[i] != -1, "ghosts should be ignored when assigning conInds.");
					if(conInds[i] != eInd){
						m_adjacencyMap.push_back(conInds[i]);
						m_connections.push_back(conElems[i_con]);
					}
				}
			}
			++ind;
		}
		assert(ind == (int)m_adjacencyMapStructure.size() - 1);
	}
	else{
	//	if ConElem::dim < Elem::dim - 1 then we have to check whether indices have
	//	already been added...
		UG_THROW("Currently a dual graph can only be created if elements are"
				" connected via their sides. Since nearly everything is prepared,"
				" implementing this step for arbitrary connecting elements shouldn't"
				" be much work.");
	}

//	add the final element
	m_adjacencyMapStructure[m_adjacencyMapStructure.size() - 1] = m_adjacencyMap.size();
	UG_DLOG(LIB_GRID, 1, "ParallelDualGraph-stop generate_graph\n");
}

}// end of namespace

#endif
