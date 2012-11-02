//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d14

#ifndef __H__LIB_GRID__PARALLEL_DUAL_GRAPH__
#define __H__LIB_GRID__PARALLEL_DUAL_GRAPH__

#include <vector>
#include "lib_grid/lg_base.h"
#include "pcl/pcl_process_communicator.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
/** This method creates a parallel dual graph for the given level only.
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
 * this array has the size of grid.num<TGeomBaseObj>(level).
 *
 * If you are interested in the indices assigned to each element, then
 * pass a pointer to an TIndexType attachment to paIndex.
 *
 * TGeomBaseObj can be either VertexBase, EdgeBase, Face or Volume
 */
template <class TGeomBaseObj, class TIndexType>
void ConstructParallelDualGraphMGLevel(
		std::vector<TIndexType>& adjacencyMapStructureOut,
		std::vector<TIndexType>& adjacencyMapOut,
		std::vector<TIndexType>& nodeOffsetMapOut,
		MultiGrid& mg, size_t level,
		pcl::ProcessCommunicator procCom,
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
	nodeOffsetMapOut[nodeOffsetMapOut.size()] = numElemsTotal;
	localNodeOffset = nodeOffsetMapOut[procCom.get_local_proc_id()];

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
			adjacencyMapOut.push_back(localNodeOffset + aaInd[vNeighbours[i]]);
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
