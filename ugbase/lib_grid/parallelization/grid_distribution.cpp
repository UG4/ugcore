//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d08

#include <algorithm>
#include "grid_distribution.h"
#include "lib_grid/lib_grid.h"
#include "pcl/pcl.h"
#include "util/distribution_util.h"
#include "common/util/stream_pack.h"
#include "common/util/binary_stream.h"

using namespace std;

namespace ug
{
/*
void PrintData(int* data, int size)
{
	for(int i = 0; i < size; ++i)
		cout << data[i] << ", ";

	cout << endl;
}
*/

////////////////////////////////////////////////////////////////////////
/**	All nodes in the given distribution layout will be added to
 *	the INT_VERTICAL_MASTER interface of layoutMapOut.*/
template <class TGeomObj>
static void AddVerticalMasterInterfaces(GridLayoutMap& layoutMapOut,
										MultiGrid& mg,
								  		const DistributionNodeLayout<TGeomObj*>& distLayout,
								  		int targetProc)
{
//	some typedefs
	typedef DistributionNodeLayout<TGeomObj*> 			DistLayout;
	typedef typename GridLayoutMap::Types<TGeomObj>::
								Layout::LevelLayout 	TLayout;
	typedef typename TLayout::Interface					TInterface;

//	access the nodes of the distribution layout
	const typename DistLayout::NodeVec& nodeVec = distLayout.node_vec();

//	we'll cache the layout and interface.
	TLayout* pLayout = NULL;
	TInterface* pInterface = NULL;
	int currentLevel = -1;

//	iterate over the nodes
	for(size_t i = 0; i < nodeVec.size(); ++i){
		TGeomObj* node = nodeVec[i];
	//	only surface-nodes are stored in vertical interfaces
		if(!mg.has_children(node)){
		//	get the level and the matching layout and interface
			int level = (int)mg.get_level(node);
			if(level != currentLevel){
				currentLevel = level;
			//	get the layout and the interface to targetProc
				pLayout = &layoutMapOut.template get_layout<TGeomObj>(INT_VERTICAL_MASTER).
															layout_on_level(level);
				pInterface = &pLayout->interface(targetProc);
			}

		//	add node to the interface
			pInterface->push_back(node);
		}
	}
}

///	adds virtual horizontal interfaces.
/**	THIS DOCU SEEMS OUTDATED!
 *	This method creates interfaces in a GridLayoutMap from interfaces
 *	in a distribution layout.
 *	Via processMap one can specify alias processes. If a alias-process is
 *	set to -1, the associated interface is ignored.
 */
template <class TGeomObj>
static void AddHorizontalInterfaces(GridLayoutMap& layoutMapOut,
								DistributionNodeLayout<TGeomObj*>& distLayout,
								std::vector<int>* pProcessMap)
{
//	some typedefs
	typedef DistributionNodeLayout<TGeomObj*> 			DistLayout;
	typedef typename GridLayoutMap::Types<TGeomObj>::
								Layout::LevelLayout 	TLayout;
	typedef typename TLayout::Interface					TInterface;

//	access the nodes of the distribution layout
	const typename DistLayout::NodeVec& nodeVec = distLayout.node_vec();

//	we'll cache the layout and interface.
	TLayout* pLayout = NULL;

//	iterate over the nodes
	for(size_t level = 0; level < distLayout.num_levels(); ++level)
	{
		typename DistLayout::InterfaceMap& imap = distLayout.interface_map(level);
		for(typename DistLayout::InterfaceMap::iterator iter = imap.begin();
			iter != imap.end(); ++iter)
		{
			TInterface* pInterface = NULL;
			int elemType = -1;

		//	get the proc-id of this interface
			int procID = iter->first;
			if(pProcessMap)
			{
				//assert((int)pProcessMap->size() > procID && "process-map to small.");
				if((int)pProcessMap->size() > procID)
					procID = (*pProcessMap)[procID];
				else
					procID = -1;
			}

			if(procID == -1)
				continue;

		//	copy the interface-elements to the layoutMap
			typename DistLayout::Interface& distInterface = iter->second;

		//	iterate over the nodes in the interface
			for(size_t i = 0; i < distInterface.size(); ++i){
				int newElemType = distInterface[i].type;
				TGeomObj* node = nodeVec[distInterface[i].localID];
			//	get the layout and the matching interface
				if(elemType != newElemType){
					elemType = newElemType;
				//	get the layout and the interface to targetProc
					pLayout = &layoutMapOut.template get_layout<TGeomObj>(elemType).
																layout_on_level(level);
					pInterface = &pLayout->interface(procID);
				}

			//	add the node to the local interface.
				pInterface->push_back(node);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
/**	This method is only called if the source-grid is kept on the distributing
 *	process. In this case all master and slave entries have to be transformed
 *	to virtual-master and virtual-slave entries.
 */
template <class TGeomObj>
static void AdjustInterfaceElementType(DistributionNodeLayout<TGeomObj*>& distLayout)
{
//	some typedefs
	typedef DistributionNodeLayout<TGeomObj*> 			DistLayout;

//	iterate over the nodes
	for(size_t level = 0; level < distLayout.num_levels(); ++level)
	{
		typename DistLayout::InterfaceMap& imap = distLayout.interface_map(level);
		for(typename DistLayout::InterfaceMap::iterator iter = imap.begin();
			iter != imap.end(); ++iter)
		{
			typename DistLayout::Interface& distInterface = iter->second;

		//	iterate over the nodes in the interface
			for(size_t i = 0; i < distInterface.size(); ++i){
			//	change the types:
			//		MASTER -> VIRTUAL_MASTER
			//		SLAVE -> VIRTUAL_SLAVE

				int elemType = distInterface[i].type;
				switch(elemType){
					case INT_MASTER: distInterface[i].type = INT_VIRTUAL_MASTER; break;
					case INT_SLAVE: distInterface[i].type = INT_VIRTUAL_SLAVE; break;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	DistributeGrid_KeepSrcGrid
bool DistributeGrid_KeepSrcGrid(MultiGrid& mg, ISubsetHandler& sh,
								GridLayoutMap& layoutMap,
								SubsetHandler& shPartition,
								int localProcID,
								std::vector<int>* pProcessMap)
{

//	we have to store the layouts for all the processes.
	vector<DistributionVertexLayout> vVertexLayouts;
	vector<DistributionEdgeLayout> vEdgeLayouts;
	vector<DistributionFaceLayout> vFaceLayouts;
	vector<DistributionVolumeLayout> vVolumeLayouts;

//	we need some attachments that will speed up the called processes.
	AInt aInt;
	mg.attach_to_vertices(aInt);
	mg.attach_to_edges(aInt);
	mg.attach_to_faces(aInt);
	mg.attach_to_volumes(aInt);

//	the selector will help to speed things up a little.
	MGSelector msel(mg);

//	if vertical interfaces shall be created then we won't distribute
//	the whole genealogy
	CreateDistributionLayouts(vVertexLayouts, vEdgeLayouts, vFaceLayouts,
							  vVolumeLayouts, mg, shPartition,
							  false, &msel);
/*
	UG_LOG("Testing Vertex Distribution Layouts:\n");
	if(!TestDistributionLayouts(vVertexLayouts))
		return false;
	UG_LOG("Testing Edge Distribution Layouts:\n");
	if(!TestDistributionLayouts(vEdgeLayouts))
		return false;
	UG_LOG("Testing Face Distribution Layouts:\n");
	if(!TestDistributionLayouts(vFaceLayouts))
		return false;
	UG_LOG("Testing Volume Distribution Layouts:\n");
	if(!TestDistributionLayouts(vVolumeLayouts))
		return false;
*/

//	we will now fill a binary stream with all the grids.
//	this stream will receive all the data that is to be sent to other processes.
	BinaryStream globalStream;
//	this vector is required so that we can use distribute-data later on.
	vector<int>	vBlockSizes;
//	here we'll store the ids of the receiving processes.
	vector<int> vReceiverIDs;

	int numProcs = (int)shPartition.num_subsets();
	if(pProcessMap)
		numProcs = std::min((int)pProcessMap->size(), numProcs);

	for(int i = 0; i < numProcs; ++i)
	{
		int proc = i;
		if(pProcessMap)
			proc = (*pProcessMap)[i];

	//	all horizontal nodes have to be transformed to virtual horizontal nodes.
		AdjustInterfaceElementType(vVertexLayouts[i]);
		AdjustInterfaceElementType(vEdgeLayouts[i]);
		AdjustInterfaceElementType(vFaceLayouts[i]);
		AdjustInterfaceElementType(vVolumeLayouts[i]);

	//	check whether the current proc is the local proc or another proc.
		if(proc == localProcID)
		{
		//	create local horizontal interfaces
			AddHorizontalInterfaces<VertexBase>(layoutMap, vVertexLayouts[i], pProcessMap);
			AddHorizontalInterfaces<EdgeBase>(layoutMap, vEdgeLayouts[i], pProcessMap);
			AddHorizontalInterfaces<Face>(layoutMap, vFaceLayouts[i], pProcessMap);
			AddHorizontalInterfaces<Volume>(layoutMap, vVolumeLayouts[i], pProcessMap);
		}
		else
		{
		//	since the original grid is kept, we have to add vertical interfaces.
		//	all surface nodes in the layouts will be vertical interface members.
		//	Here we create the local vertical interfaces.
			AddVerticalMasterInterfaces<VertexBase>(layoutMap, mg,
													vVertexLayouts[i], proc);
			AddVerticalMasterInterfaces<EdgeBase>(layoutMap, mg,
												  vEdgeLayouts[i], proc);
			AddVerticalMasterInterfaces<Face>(layoutMap, mg,
											  vFaceLayouts[i], proc);
			AddVerticalMasterInterfaces<Volume>(layoutMap, mg,
											  vVolumeLayouts[i], proc);

		//	serialize the parts of the grid that will be send to other processes.
			int oldSize = globalStream.size();
			SerializeGridAndDistributionLayouts(
									globalStream, mg, vVertexLayouts[i],
									vEdgeLayouts[i], vFaceLayouts[i], vVolumeLayouts[i],
									aInt, aInt, aInt, aInt, &msel, pProcessMap);

		//	serialize subset indices
			SerializeSubsetHandler(mg, sh,
								   msel.get_geometric_object_collection(),
								   globalStream);

		//	serialize position attachment
			for(uint iLevel = 0; iLevel < mg.num_levels(); ++iLevel)
			{
				SerializeAttachment<VertexBase>(mg, aPosition,
												msel.begin<VertexBase>(iLevel),
												msel.end<VertexBase>(iLevel),
												globalStream);
			}
//TODO:		the user should be able to add personal data to those buffers.

			vBlockSizes.push_back((int)(globalStream.size() - oldSize));
			vReceiverIDs.push_back(proc);
		}
	}

//	send the grids to their target processes
	int numReceivers = (int)vReceiverIDs.size();
	if(numReceivers > 0)
	{
	//	every process receives the size of the data-buffer first.
		vector<int> bufferSizes(numReceivers, sizeof(int));

	//	distribute the block-sizes to the different processes
		pcl::DistributeData(localProcID, &vReceiverIDs.front(), numReceivers,
							&vBlockSizes.front(), &bufferSizes.front(), 38);

	//	distribute the grids-distribution-packs
		pcl::DistributeData(localProcID, &vReceiverIDs.front(), numReceivers,
							globalStream.buffer(), &vBlockSizes.front(), 39);
	}

//	clean up
	mg.detach_from_vertices(aInt);
	mg.detach_from_edges(aInt);
	mg.detach_from_faces(aInt);
	mg.detach_from_volumes(aInt);

	return true;
}

////////////////////////////////////////////////////////////////////////
//	DistributeGrid
bool DistributeGrid(MultiGrid& mg, ISubsetHandler& sh,
					SubsetHandler& shPartition,
					int localProcID, MultiGrid* pLocalGridOut,
					ISubsetHandler* pLocalSHOut,
					GridLayoutMap* pLocalGridLayoutMapOut,
					std::vector<int>* pProcessMap)
{

//	we have to store the layouts for all the processes.
	vector<DistributionVertexLayout> vVertexLayouts;
	vector<DistributionEdgeLayout> vEdgeLayouts;
	vector<DistributionFaceLayout> vFaceLayouts;
	vector<DistributionVolumeLayout> vVolumeLayouts;

//	we need some attachments that will speed up the called processes.
	AInt aInt;
	mg.attach_to_vertices(aInt);
	mg.attach_to_edges(aInt);
	mg.attach_to_faces(aInt);
	mg.attach_to_volumes(aInt);

//	the selector will help to speed things up a little.
	MGSelector msel(mg);

	CreateDistributionLayouts(vVertexLayouts, vEdgeLayouts, vFaceLayouts,
							  vVolumeLayouts, mg, shPartition,
							  true, &msel);

//	we will now fill a binary stream with all the grids.
//	this stream will receive the data that has to be copied to the local grid.
	BinaryStream localStream;
//	this stream will receive all the data that is to be sent to other processes.
	BinaryStream globalStream;
//	this vector is required so that we can use distribute-data later on.
	vector<int>	vBlockSizes;
//	here we'll store the ids of the receiving processes.
	vector<int> vReceiverIDs;

	int numProcs = (int)shPartition.num_subsets();
	if(pProcessMap)
		numProcs = std::min((int)pProcessMap->size(), numProcs);

	for(int i = 0; i < numProcs; ++i)
	{
/*
cout << "proc " << i << ":\n";
cout << "  layouts:\n";
cout << "    vrts: " << vVertexLayouts[i].node_vec().size() << endl;
cout << "    edges: " << vEdgeLayouts[i].node_vec().size() << endl;
cout << "    faces: " << vFaceLayouts[i].node_vec().size() << endl;
cout << "    vols: " << vVolumeLayouts[i].node_vec().size() << endl;
*/
		int proc = i;
		if(pProcessMap)
			proc = (*pProcessMap)[i];

		if(proc == localProcID)
		{
			SerializeGridAndDistributionLayouts(
								localStream, mg, vVertexLayouts[i],
								vEdgeLayouts[i], vFaceLayouts[i], vVolumeLayouts[i],
								aInt, aInt, aInt, aInt, &msel, pProcessMap);

		//	serialize subset indices
			SerializeSubsetHandler(mg, sh,
								   msel.get_geometric_object_collection(),
								   localStream);

		//	serialize position attachment
			for(uint iLevel = 0; iLevel < mg.num_levels(); ++iLevel)
			{
				SerializeAttachment<VertexBase>(mg, aPosition,
												msel.begin<VertexBase>(iLevel),
												msel.end<VertexBase>(iLevel),
												localStream);
			}

//TODO:		the user should be able to add personal data to those buffers.
		}
		else
		{
			int oldSize = globalStream.size();
			SerializeGridAndDistributionLayouts(
									globalStream, mg, vVertexLayouts[i],
									vEdgeLayouts[i], vFaceLayouts[i], vVolumeLayouts[i],
									aInt, aInt, aInt, aInt, &msel, pProcessMap);

		//	serialize subset indices
			SerializeSubsetHandler(mg, sh,
								   msel.get_geometric_object_collection(),
								   globalStream);

		//	serialize position attachment
			for(uint iLevel = 0; iLevel < mg.num_levels(); ++iLevel)
			{
				SerializeAttachment<VertexBase>(mg, aPosition,
												msel.begin<VertexBase>(iLevel),
												msel.end<VertexBase>(iLevel),
												globalStream);
			}
//TODO:		the user should be able to add personal data to those buffers.

			vBlockSizes.push_back((int)(globalStream.size() - oldSize));
			vReceiverIDs.push_back(proc);
		}
	}

//	send the grids to their target processes
	int numReceivers = (int)vReceiverIDs.size();
	if(numReceivers > 0)
	{
	//	every process receives the size of the data-buffer first.
		vector<int> bufferSizes(numReceivers, sizeof(int));

	//	distribute the block-sizes to the different processes
		pcl::DistributeData(localProcID, &vReceiverIDs.front(), numReceivers,
							&vBlockSizes.front(), &bufferSizes.front(), 38);

	//	distribute the grids-distribution-packs
		pcl::DistributeData(localProcID, &vReceiverIDs.front(), numReceivers,
							globalStream.buffer(), &vBlockSizes.front(), 39);
	}

//	fill the local grid and subset-handler
	if(pLocalGridOut && pLocalSHOut &&
		pLocalGridLayoutMapOut && (localStream.size() > 0))
	{
		if(!pLocalGridOut->has_vertex_attachment(aPosition))
			pLocalGridOut->attach_to_vertices(aPosition);

		DeserializeGridAndDistributionLayouts(
									*pLocalGridOut, *pLocalGridLayoutMapOut,
									localStream);

	//	deserialize subset handler
		if(!DeserializeSubsetHandler(*pLocalGridOut, *pLocalSHOut,
									pLocalGridOut->get_geometric_object_collection(),
									localStream))
		{
			goto bailout_false;
		}

		for(size_t i = 0; i < pLocalGridOut->num_levels(); ++i){
			if(!DeserializeAttachment<VertexBase>(*pLocalGridOut, aPosition,
											pLocalGridOut->begin<VertexBase>(i),
											pLocalGridOut->end<VertexBase>(i),
											localStream))
				goto bailout_false;
		}
	}

//	clean up
//bailout_true:
	mg.detach_from_vertices(aInt);
	mg.detach_from_edges(aInt);
	mg.detach_from_faces(aInt);
	mg.detach_from_volumes(aInt);

	return true;

bailout_false:
	mg.detach_from_vertices(aInt);
	mg.detach_from_edges(aInt);
	mg.detach_from_faces(aInt);
	mg.detach_from_volumes(aInt);

	return false;

}


////////////////////////////////////////////////////////////////////////
/**	All nodes in the given distribution layout will be added to
 *	the INT_VERTICAL_MASTER interface of layoutMapOut.
 *	This method is a little unsafe...
 *	It is absolutly crucial, that all elements have been created
 *	in the same order as they were received through the in-stream
 *	in ReceiveInitialGrid. Especially auto-create options can be
 *	problematic here. Make sure to avoid auto-element-creation.
 *	This method currently requires a non-adaptive grid (Only as long
 *	as there are dummy-parents).*/
template <class TGeomObj>
static void AddVerticalSlaveInterfaces(GridLayoutMap& layoutMapOut,
									   MultiGrid& mg,
									   int srcProc)
{
//	some typedefs
	typedef typename GridLayoutMap::Types<TGeomObj>::
								Layout::LevelLayout 	TLayout;
	typedef typename TLayout::Interface					TInterface;
	typedef typename geometry_traits<TGeomObj>::iterator	ObjIter;

//todo:	iterate over all levels of the multi-grid
//	iterate over the toplevel of the multigrid only.
//	This is required since we're using dummy-parents (THROW THEM OUT!).
	if(mg.num_levels() == 0)
		return;

	//for(size_t level = 0; level < mg.num_levels(); ++level){
	for(size_t level = mg.num_levels() - 1; level < mg.num_levels(); ++level){
	//	if there are no nodes in this level, we can return immediatly
		if(mg.num<TGeomObj>(level) == 0)
			break;

	//	get the appropriate layout and interface
		TLayout* pLayout = NULL;
		TInterface* pInterface = NULL;

	//	iterate over the elements of the level
		for(ObjIter iter = mg.begin<TGeomObj>(level);
			iter != mg.end<TGeomObj>(level); ++iter)
		{
		//TODO: this will change when we no longer send parents along with children.
		//	make sure that the node is a surface node
			if(!mg.has_children(*iter)){
			//	get the appropriate layout and interface
				if(!pInterface){
					pLayout = &layoutMapOut.template
								get_layout<TGeomObj>(INT_VERTICAL_SLAVE).
									layout_on_level(level);

					pInterface = &pLayout->interface(srcProc);
				}

				pInterface->push_back(*iter);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	ReceiveGrid
bool ReceiveGrid(MultiGrid& mgOut, ISubsetHandler& shOut,
				 GridLayoutMap& gridLayoutMapOut,
				 int srcProcID, bool createVerticalLayouts)
{
//	receive the stream-size
	int streamSize;
	pcl::ReceiveData(&streamSize, srcProcID, sizeof(int), 38);

//	receive the buffer
	BinaryStream binaryStream(streamSize);
	pcl::ReceiveData(binaryStream.buffer(), srcProcID, streamSize, 39);

//	fill the grid and the layout
	DeserializeGridAndDistributionLayouts(mgOut, gridLayoutMapOut,
											binaryStream);

//	if vertical layouts shall be created, do it now.
//	note that only surface-nodes are assigned to vertical interfaces.
	if(createVerticalLayouts){
		AddVerticalSlaveInterfaces<VertexBase>(gridLayoutMapOut, mgOut, srcProcID);
		AddVerticalSlaveInterfaces<EdgeBase>(gridLayoutMapOut, mgOut, srcProcID);
		AddVerticalSlaveInterfaces<Face>(gridLayoutMapOut, mgOut, srcProcID);
		AddVerticalSlaveInterfaces<Volume>(gridLayoutMapOut, mgOut, srcProcID);
	}

//	deserialize subset handler
	if(!DeserializeSubsetHandler(mgOut, shOut,
								mgOut.get_geometric_object_collection(),
								binaryStream))
		return false;

//	read the attached data
	if(!mgOut.has_vertex_attachment(aPosition))
		mgOut.attach_to_vertices(aPosition);

	for(size_t i = 0; i < mgOut.num_levels(); ++i){
		if(!DeserializeAttachment<VertexBase>(mgOut, aPosition,
											mgOut.begin<VertexBase>(i),
											mgOut.end<VertexBase>(i),
											binaryStream))
			return false;
	}
/*
	if(pcl::IsOutputProc())
		SaveGridToFile(mgOut, "tmpOutProcHierarchy.ugx",
						mgOut.get_hierarchy_handler());
*/
//TODO:	allow the user to read his data.

	return true;
}



//	processes which have to communicate with each other are
//	notified through this method.
//	vProcRanksInOut will grow if other processes want to
//	communicate with this process.
//	todo: pass a communicator. Currently COMM_WORLD is used.
void CommunicateInvolvedProcesses(vector<int>& vReceiveFromRanksOut,
								  vector<int>& vSendToRanks,
								  const pcl::ProcessCommunicator& procComm
								  	= pcl::ProcessCommunicator())
{
	using namespace pcl;

//	if there is only one process in the communicator, there's
//	nothing to do.
	if(procComm.size() < 2)
		return;

	const int localProcRank = GetProcRank();

//	we'll use an array in which we'll store the number of
//	processes, to which each process wants to talk.
	vector<int> vNumAssProcs(procComm.size());
	UG_LOG("  procCommSize: " << procComm.size() << endl);

//	perform allgather with proc-numbers
	int procCount = (int)vSendToRanks.size();

	procComm.allgather(&procCount, 1, PCL_DT_INT,
					   &vNumAssProcs.front(),
					   1,
					   PCL_DT_INT);

//	sum the number of processes so that the absolute list size can
//	be determined.
	size_t listSize = 0;
	vector<int> vDisplacements(vNumAssProcs.size());
	for(size_t i = 0; i < vNumAssProcs.size(); ++i){
		vDisplacements[i] = (int)listSize;
		listSize += vNumAssProcs[i];
	}

	UG_LOG("  list-size: " << listSize << endl);

//	perform allgather with proc-lists
//	this list will later on contain an adjacency-list for each proc.
//	adjacency in this case means, that processes want to communicate.
	vector<int> vGlobalProcList(listSize);

	procComm.allgatherv(&vSendToRanks.front(),
						procCount,
						PCL_DT_INT,
						&vGlobalProcList.front(),
						&vNumAssProcs.front(),
						&vDisplacements.front(),
						PCL_DT_INT);

//	we can now check for each process whether it wants to
//	communicate with this process. Simply iterate over
//	the adjacency list to do so.
	for(int i = 0; i < (int)vNumAssProcs.size(); ++i)
	{
		if(i != localProcRank){
		//	check whether the i-th proc wants to communicate with the local proc
			for(int j = 0; j < vNumAssProcs[i]; ++j)
			{
				if(vGlobalProcList[vDisplacements[i] + j] == localProcRank)
				{
				//	check whether vProcRanksInOut already contains the entry
					if(find(vReceiveFromRanksOut.begin(), vReceiveFromRanksOut.end(), i)
							== vReceiveFromRanksOut.end())
					{
						vReceiveFromRanksOut.push_back(i);
					}

				//	the i-th proc is handled completly. resume with the next.
					break;
				}
			}
		}
	}

//	vProcRanksInOut should now contain all process-ranks with which
//	the local proc should communicate.
}

///	writes interfaces-changes into a binary-stream-pack
/**
 * Data will be written into streamPackSendOut. This method won't
 * clear stremPackSendOut, since one might call it several times
 * in order to use one stream for all interface types.
 * You'll have to clear it yourself, if you want to start with
 * empty streams.
 *
 * Each stream will contain a series of interfaces afterwards.
 * Layout is as follows:
 * {proc type size level {indices}}
 *
 * Where the entries have the following meaning:
 *		- int proc: The connected proc id.
 *		- int type: The interface type on the connected proc.
 *		- int size: The number of entries in the following indexlist.
 *		- int level: The level on which the interface lies.
 *		- int-array indices: Indices refer to the interface that connects the
 *						 local process with the process from which data was received.
 *						 The interface-type of the referred interface is that, which is
 *						 connected to the type specified above.
 */
template <class TLayout, class TDistLayout>
void PrepareNewInterfaces(StreamPack& streamPackSendOut,
						DistributedGridManager& distGridMgr,
					   	TLayout& layout,
					   	int interfaceType,
					   	std::vector<TDistLayout>& distLayoutVec)
{
	typedef typename TLayout::iterator	InterfaceIter;
	typedef typename TLayout::Interface Interface;

//	the local procID
	int localProcID = pcl::GetProcRank();

//	we need a stream pack, that holds a stream for each neighbour.
	StreamPack& streamPack = streamPackSendOut;

//	touch each neighbor process
	vector<int> vNeighborProcs;
	pcl::CollectAssociatedProcesses(vNeighborProcs, layout);
	for(size_t i = 0; i < vNeighborProcs.size(); ++i){
		streamPack.get_stream(vNeighborProcs[i]);
	}

//	data in each buffer is organized as follows:
//		- int srcProc (the process on which the interface is created)
//		- int interfaceType	(the type of the interface on the srcPrc)
//		- int interfaceSize
//		- {elementIndices}

//	iterate through the processes to which interfaces have to be build.
	for(size_t iDistLayout = 0; iDistLayout < distLayoutVec.size();
		++iDistLayout)
	{
	//	the local proc can be ignored.
		if((int)iDistLayout == localProcID)
			continue;

		TDistLayout& distLayout = distLayoutVec[iDistLayout];

	//	iterate through the levels
		for(int level = 0; level < (int)distLayout.num_levels(); ++level)
		{
			typename TDistLayout::InterfaceMap imap =
				distLayout.interface_map(level);
		//	iterate through the entries of the imap
			for(typename TDistLayout::InterfaceMap::iterator iter = imap.begin();
				iter != imap.end(); ++iter)
			{
				int procID = iter->first;
			//	interfaces to the local process are ignored
				if(procID == localProcID)
					continue;

			//	check whether the interface points to a neighbour
				if(layout.interface_exists(procID, level)){
				//	we then have to write the assiciated indices into a buffer and send it
				//	along the other data to procID.
					typename TDistLayout::Interface& distInterface = iter->second;
					BinaryStream& stream = *streamPack.get_stream(procID);

				//	first write the stream-header, then write the associated indices
				//	(those indices are the indices of each node in the existing pcl-interface)
					int interfaceSize = (int)distInterface.size();
					if(interfaceSize > 0){
						stream.write((char*)&iDistLayout, sizeof(int));
						stream.write((char*)&interfaceType, sizeof(int));
						stream.write((char*)&interfaceSize, sizeof(int));
						stream.write((char*)&level, sizeof(int));

						for(size_t i = 0; i < distInterface.size(); ++i){
						//todo	In the moment all interfaceIDs are from 1, ..., n.
						//		However - as soon as we allow dynamic removal of interface
						//		elements, those ids do not match the indices.
						//		This means we have to refresh those ids before redistribution
						//		starts.
							int index = distInterface[i].associatedInterfaceID - 1;
							stream.write((char*)&index, sizeof(int));
						}
					}
				}
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////
///	communicates with other processes and removes interface entries as requested.
/**
 * \param distGridMgr: The distributed grid manager
 * \param sel: A Selector that holds all elements that stay on the process.
 * \param distLayoutVec: a vector that contains the distribution layouts,
 *						 that are to be sent to other processes.
 *
 * \returns	true, if the layout-map has changed, false if not.
 */
template <class TLayout, class TSelector, class TDistLayout>
bool UpdateInterfaces(DistributedGridManager& distGridMgr, TSelector& sel,
					  std::vector<TDistLayout>& distLayoutVec)

{
	using namespace pcl;

//	typedefs
	typedef typename TLayout::Type GeomObjType;
	typedef typename TLayout::iterator	InterfaceIter;
	typedef typename TLayout::Interface Interface;

//	retrieve the multi-grid
	if(!sel.get_assigned_grid())
		return false;

	Grid& grid = *sel.get_assigned_grid();
	GridLayoutMap& glm = distGridMgr.grid_layout_map();

////////////////////////////////
//	We have to communicate with all direct neighbours (linked by an interface),
//	how an existing interface has changed.
//	To do this we need a selector that holds all elements that stay on
//	the local proc, the old interfaces and the distribution-layouts, that
//	hold the new interfaces.
//todo: Take care of direct neighbours that receive data from this process.
//todo: Take care of vertical interfaces
//todo: selRemainingEntries should only operate on GeomObjType

	Selector selRemainingEntries(grid);
//	select all entries in sel initially
	for(size_t i = 0; i < sel.num_levels(); ++i)
	{
		selRemainingEntries.select(sel.template begin<GeomObjType>(i),
								   sel.template end<GeomObjType>(i));
	}

	SelectionCommPol<TLayout, TSelector, Selector>
		selCommPol(sel, selRemainingEntries);
	ParallelCommunicator<TLayout> communicator;

	communicator.exchange_data(glm, INT_MASTER, INT_SLAVE, selCommPol);
	communicator.exchange_data(glm, INT_SLAVE, INT_MASTER, selCommPol);


//	we have to tell the neighbours, which new interfaces they have to build.
//	indices are relative to the interface between localProc and its neighbour.
//	iterate over all interfaces
	StreamPack streamPackSend;
	StreamPack streamPackReceive;

	for(typename GridLayoutMap::Types<GeomObjType>::Map::iterator iter =
		glm.template layouts_begin<GeomObjType>();
		iter != glm.template layouts_end<GeomObjType>(); ++iter)
	{
		UG_LOG("  interface type: " << iter->first << endl);
		PrepareNewInterfaces(streamPackSend, distGridMgr,
							 iter->second, iter->first,
						  	 distLayoutVec);
	}

//	now write each stream into the communicator
//	at the same time we'll shedule a receive
//	this is possible and required, since streamPackSend contains
//	all links to direct neighbours (no other links).
	for(StreamPack::iterator iter = streamPackSend.begin();
		iter != streamPackSend.end(); ++iter)
	{
		UG_LOG("  sending raw: " << iter->second->size() << endl);

		communicator.send_raw(iter->first, iter->second->buffer(),
							  iter->second->size(), false);

		communicator.receive_raw(iter->first, *streamPackReceive.get_stream(iter->first), -1);
	}

//	communicate the data
	communicator.communicate();

//	this vector will be used to create the new interfaces
	vector<typename Interface::Element> vInterfaceElements;

//	we'll create interface entries to other processes now
//	we need a temporary GridLayoutMap, since we have to avoid problems,
//	when it comes to erasing old unused interface entries later on.
	GridLayoutMap tmpGlm;

//	iterate over the streams and process data
	for(StreamPack::iterator iter = streamPackReceive.begin();
		iter != streamPackReceive.end(); ++iter)
	{
		BinaryStream& stream = *iter->second;
//		UG_LOG("  stream.size: " << stream.size() << endl);
		while(stream.can_read_more()){
		//	read the header
			int newConnectedProcID;
			int newConnectedInterfaceType;
			int newInterfaceSize;
			int level;
			stream.read((char*)&newConnectedProcID, sizeof(int));
			stream.read((char*)&newConnectedInterfaceType, sizeof(int));
			stream.read((char*)&newInterfaceSize, sizeof(int));
			stream.read((char*)&level, sizeof(int));
/*
			UG_LOG("  newConnectedProcID: " << newConnectedProcID << endl);
			UG_LOG("  newConnectedInterfaceType: " << newConnectedInterfaceType << endl);
			UG_LOG("  newInterfaceSize: " << newInterfaceSize << endl);
*/
		//	make sure that there are nodes to put into the new interface
			if(newInterfaceSize > 0){
			//	we have to get the local interface-type
				int newInterfaceType = GetAssociatedInterfaceType(newConnectedInterfaceType);

			//	access the old and the new interface
				TLayout& oldLayout = glm.template get_layout<GeomObjType>(newInterfaceType);
				TLayout& newLayout = tmpGlm.template get_layout<GeomObjType>(newInterfaceType);
				Interface& oldInterface = oldLayout.interface(iter->first, level);
				Interface& newInterface = newLayout.interface(newConnectedProcID, level);

//TODO:	ONLY ELEMENTS THAT STAY ON THE LOCAL PROC MAY BE INSERTED INTO INTERFACES.

			//	we'll store all nodes of oldInterface in a temporary vector,
			//	since we want to access them by index
				vInterfaceElements.clear();
				for(typename Interface::iterator elemIter = oldInterface.begin();
					elemIter != oldInterface.end(); ++elemIter)
					vInterfaceElements.push_back(oldInterface.get_element(elemIter));

			//	get the nodes from oldInterface and add them to newInterface
				for(int i = 0; i < newInterfaceSize; ++i){
					int tInd;
					stream.read((char*)&tInd, sizeof(int));
					UG_LOG("interface entry: " << tInd << endl);
					newInterface.push_back(vInterfaceElements[tInd]);
				}
			}
		}
	}

//	we have to remove entries which are not part of sel.
	for(typename geometry_traits<GeomObjType>::iterator iter =
		selRemainingEntries.template begin<GeomObjType>();
		iter != selRemainingEntries.template end<GeomObjType>(); ++iter)
	{
		if(!sel.is_selected(*iter))
			selRemainingEntries.deselect(*iter);
	}

//	remove interface entries that are no longer required
	bool retVal = RemoveUnselectedInterfaceEntries<GeomObjType>(glm, selRemainingEntries);

//	now add the elements that lie in tmpGlm
	for(typename GridLayoutMap::Types<GeomObjType>::Map::iterator iter =
		tmpGlm.template layouts_begin<GeomObjType>();
		iter != tmpGlm.template layouts_end<GeomObjType>(); ++iter)
	{
		TLayout& tmpLayout = iter->second;
		TLayout& layout = glm.get_layout<GeomObjType>(iter->first);

	//	add all elements from tmpLayout to layout
		for(size_t lvl = 0; lvl < tmpLayout.num_levels(); ++lvl)
		{
		//	iterate over all interfaces
			for(InterfaceIter iIter = tmpLayout.begin(lvl);
				iIter != tmpLayout.end(lvl); ++iIter)
			{
				Interface& tmpInterface = tmpLayout.interface(iIter);
				Interface& interface = layout.interface(tmpLayout.proc_id(iIter), lvl);

			//	add the elements
				for(typename Interface::iterator iter = tmpInterface.begin();
					iter != tmpInterface.end(); ++iter)
				{
					interface.push_back(tmpInterface.get_element(iter));
				}
			}
		}
	}

	return retVal;
}

bool RedistributeGrid(DistributedGridManager& distGridMgrInOut,
					  ISubsetHandler& shInOut,
					  SubsetHandler& shPartition)
{
	using namespace pcl;

	DistributedGridManager& distGridMgr = distGridMgrInOut;

	if(!distGridMgr.get_assigned_grid()){
		UG_LOG("  WARNING in RedistributeGrid: distGridMgrInOut is not associated with a grid. Aborting.\n");
		return false;
	}

	UG_LOG("redistributing grid:\n");

//	create the communicator object
//	todo: this should be passed to the method by a parameter
	ProcessCommunicator procComm;//	WORLD by default.

	const int localProcRank = GetProcRank();

	MultiGrid& mg = *distGridMgr.get_assigned_grid();
	ISubsetHandler& sh = shInOut;
	GridLayoutMap& glm = distGridMgr.grid_layout_map();

//	we have to store the layouts for all the processes.
	vector<DistributionVertexLayout> vVertexLayouts;
	vector<DistributionEdgeLayout> vEdgeLayouts;
	vector<DistributionFaceLayout> vFaceLayouts;
	vector<DistributionVolumeLayout> vVolumeLayouts;

//	we need some attachments that will speed up the called processes.
	AInt aInt;
	mg.attach_to_vertices(aInt);
	mg.attach_to_edges(aInt);
	mg.attach_to_faces(aInt);
	mg.attach_to_volumes(aInt);

//	the selector will help to speed things up a little.
	MGSelector msel(mg);


//todo	before CreateRedistributionLayouts is called, we have to make sure,
//		that the localIDs of elements in interfaces are continuos, starting
//		from 1. This is important, since they are used as indices later on.
//		In the moment this is automatically given. As soon as dynamic
//		erasure of interface elements is allowed, this will lead to problems
//		if not adressed.

//TODO: Only prepare layouts for processes to which data is actually send.
	CreateRedistributionLayouts(vVertexLayouts, vEdgeLayouts, vFaceLayouts,
								vVolumeLayouts, distGridMgr, shPartition,
							  	false, &msel);

////////////////////////////////
//	SETUP COMMUNICATION PATTERNS
//	first we have to communicate which process has to wait for data
//	from which other processes. In order to avoid a m^2 communication,
//	we'll first communicate for each process, with how many processes it
//	wants to communicate, and in the second step we'll distribute a list
//	that contains the ranks of those processes.

//	iterate over all partitions
//	if the i-th partition is not empty, we have to communicate with
//	the process of rank i.
//TODO: use a proc-map to decouple subsets-indices and proces-ranks.

	vector<int> vSendToRanks;
	int highestSendToRank = -1;
	for(int si = 0; si < shPartition.num_subsets(); ++si)
	{
		//	todo: we have to check whether interfaces change and notify the
		//		  associated processes.
		if(si != localProcRank){
		//	check whether there is anything in the partition
			if(!shPartition.empty(si)){
				vSendToRanks.push_back(si);
				if(si > highestSendToRank)
					highestSendToRank = si;
			/*
			//	log
				UG_LOG("  sending " << shPartition.num<Face>(si)
						<< " faces to process " << si << endl);
			*/
			}
		}
	}

//	send and receive communication requests from other processes.
	vector<int> vReceiveFromRanks;
	CommunicateInvolvedProcesses(vReceiveFromRanks, vSendToRanks, procComm);
/*
//	log the comm ranks
	UG_LOG("  sending to ranks: ");
	for(size_t i = 0; i < vSendToRanks.size(); ++i){
		UG_LOG(vSendToRanks[i] << ", ");
	}
	UG_LOG(endl);

	UG_LOG("  receiving from ranks: ");
	for(size_t i = 0; i < vReceiveFromRanks.size(); ++i){
		UG_LOG(vReceiveFromRanks[i] << ", ");
	}
	UG_LOG(endl);
*/
////////////////////////////////
//	ADJUST LOCAL INTERFACES AND ERASE UNUSED GRID PARTS
//	COMMUNICATE INTERFACE CHANGES TO DIRECT NEIGHBOURS
//	adjust local interfaces and erase obsolete grid-data.
	int eraseOption = -1;
	if(vSendToRanks.size() > 0)
	{
	//	erase all elements that will not remain.
		if((int)vVertexLayouts.size() > localProcRank)
		{
			msel.clear();
			SelectNodesInLayout(msel, vVertexLayouts[localProcRank]);
			SelectNodesInLayout(msel, vEdgeLayouts[localProcRank]);
			SelectNodesInLayout(msel, vFaceLayouts[localProcRank]);
			SelectNodesInLayout(msel, vVolumeLayouts[localProcRank]);

			eraseOption = 0;
		}
		else{
		//	if there is no entry in vVertexLayouts, we can remove all elements
			msel.clear();

			eraseOption = 1;
		}
	}
	else{
	//	we have to receive interface changes from other processes.
		msel.select(mg.begin<VertexBase>(), mg.end<VertexBase>());
		msel.select(mg.begin<EdgeBase>(), mg.end<EdgeBase>());
		msel.select(mg.begin<Face>(), mg.end<Face>());
		msel.select(mg.begin<Volume>(), mg.end<Volume>());

		eraseOption = 2;
	}

	bool bErasedElementsFromLayout = false;
	bErasedElementsFromLayout |= UpdateInterfaces<VertexLayout>(distGridMgr, msel,
															vVertexLayouts);
	bErasedElementsFromLayout |= UpdateInterfaces<EdgeLayout>(distGridMgr, msel,
															vEdgeLayouts);
	bErasedElementsFromLayout |= UpdateInterfaces<FaceLayout>(distGridMgr, msel,
															vFaceLayouts);
	bErasedElementsFromLayout |= UpdateInterfaces<VolumeLayout>(distGridMgr, msel,
															vVolumeLayouts);

////////////////////////////////
//	PREPARE SEND DATA FOR GRIDS
//	pack the parts of the grid that are to be sent to other processes into
//	a binary stream.
//	WARNING: Target processes have to be empty in the current version!
//todo: Allow to send grids to processes that already have a grid.

//	iterate over all processes to which we want to send data
	BinaryStream sendStream;
	vector<int> vSendBlockSizes;

	for(size_t i = 0; i < vSendToRanks.size(); ++i){
	//	mark select the part of the grid that shall be sent to the process
		int destProc = vSendToRanks[i];
		msel.clear();
		SelectNodesInLayout(msel, vVertexLayouts[destProc]);
		SelectNodesInLayout(msel, vEdgeLayouts[destProc]);
		SelectNodesInLayout(msel, vFaceLayouts[destProc]);
		SelectNodesInLayout(msel, vVolumeLayouts[destProc]);

		size_t oldSize = sendStream.size();

	//	todo: write interface changes

	//	todo: write interfaces

	//	write the grid.
	//	during serialization the local indices are automatically generated
	//	and written to the aLocalInd... attachments.
		SerializeMultiGridElements(mg,
							msel.get_geometric_object_collection(),
							aInt, aInt,
							aInt, aInt, sendStream);

	//	serialize subset indices
		SerializeSubsetHandler(mg, sh,
							   msel.get_geometric_object_collection(),
							   sendStream);

	//	serialize interfaces
		SerializeDistributionLayoutInterfaces(sendStream, vVertexLayouts[destProc]);
		SerializeDistributionLayoutInterfaces(sendStream, vEdgeLayouts[destProc]);
		SerializeDistributionLayoutInterfaces(sendStream, vFaceLayouts[destProc]);
		SerializeDistributionLayoutInterfaces(sendStream, vVolumeLayouts[destProc]);

	//	serialize position attachment
	//TODO: take the position-attachment as a template parameter
		for(uint iLevel = 0; iLevel < mg.num_levels(); ++iLevel)
		{
			SerializeAttachment<VertexBase>(mg, aPosition,
											msel.begin<VertexBase>(iLevel),
											msel.end<VertexBase>(iLevel),
											sendStream);
		}

		vSendBlockSizes.push_back((int)(sendStream.size() - oldSize));
	}

////////////////////////////////
//	erase unused grid parts
	if(vSendToRanks.size() > 0)
	{
	//	erase all elements that will not remain.
		if((int)vVertexLayouts.size() > localProcRank)
		{
			msel.clear();
			SelectNodesInLayout(msel, vVertexLayouts[localProcRank]);
			SelectNodesInLayout(msel, vEdgeLayouts[localProcRank]);
			SelectNodesInLayout(msel, vFaceLayouts[localProcRank]);
			SelectNodesInLayout(msel, vVolumeLayouts[localProcRank]);

			eraseOption = 0;
		}
		else{
		//	if there is no entry in vVertexLayouts, we can remove all elements
			msel.clear();

			eraseOption = 1;
		}
	}
	else{
	//	we have to receive interface changes from other processes.
		msel.select(mg.begin<VertexBase>(), mg.end<VertexBase>());
		msel.select(mg.begin<EdgeBase>(), mg.end<EdgeBase>());
		msel.select(mg.begin<Face>(), mg.end<Face>());
		msel.select(mg.begin<Volume>(), mg.end<Volume>());

		eraseOption = 2;
	}

	if(bErasedElementsFromLayout){
		switch(eraseOption){
			case 0:{
			//	erase data
				InvertSelection(msel);
				EraseSelectedObjects(msel);
			}break;

			case 1:{
				mg.clear_geometry();
//todo: clear the glm.
				//glm.clear();
			}break;

			default:{
			}break;
		}
	}

////////////////////////////////
//	COMMUNICATE DATA
//	we'll use this tag for communication
	int redistTag = 23;

//	send and receive the data from associated processes
//	first we have to communicate the block sizes
	vector<int> vRecvBlockSizes(vReceiveFromRanks.size());
	{
		vector<int> vRecvSizes(vRecvBlockSizes.size(), sizeof(int));
		vector<int> vSendSizes(vSendBlockSizes.size(), sizeof(int));
		procComm.distribute_data(&vRecvBlockSizes.front(),
								 &vRecvSizes.front(),
								 &vReceiveFromRanks.front(),
								 (int)vReceiveFromRanks.size(),
								 &vSendBlockSizes.front(),
								 &vSendSizes.front(),
								 &vSendToRanks.front(),
								 (int)vSendToRanks.size(),
								 redistTag);
	}
/*
//	log the sizes of the data-bocks that will be received
	for(size_t i = 0; i < vReceiveFromRanks.size(); ++i)
	{
		UG_LOG("  receiving " << vRecvBlockSizes[i]
				<< " bytes from process " << vReceiveFromRanks[i] << endl);
	}
*/
//	prepare receive buffers
	BinaryStream recvStream;
	{
		size_t totalSize = 0;
		for(size_t i = 0; i < vRecvBlockSizes.size(); ++i)
			totalSize += vRecvBlockSizes[i];
		recvStream.resize(totalSize);
	}

//	now distribute the data
	procComm.distribute_data(recvStream.buffer(), &vRecvBlockSizes.front(),
							 &vReceiveFromRanks.front(),
							 (int)vReceiveFromRanks.size(),
							 sendStream.buffer(), &vSendBlockSizes.front(),
							 &vSendToRanks.front(),
							 (int)vSendToRanks.size(),
							 redistTag);

//	deserialize the grid
//WARNING: the current implementation only works as long as
//			the receiving processes initially contains empty grids.
//TODO: implement grid-merging
	if(recvStream.size() > 0)
	{
	//	vectors in which we'll record received elements
		vector<VertexBase*> vrts;
		vector<EdgeBase*> edges;
		vector<Face*> faces;
		vector<Volume*> vols;

	//	deserialize the multi-grid
		DeserializeMultiGridElements(mg, recvStream, &vrts, &edges, &faces, &vols);

	//	deserialize subset handler
		if(!DeserializeSubsetHandler(mg, sh,
									mg.get_geometric_object_collection(),
									recvStream))
			return false;

	//	deserialize the distribution layouts
		DeserializeDistributionLayoutInterfaces(glm, vrts, recvStream);
		DeserializeDistributionLayoutInterfaces(glm, edges, recvStream);
		DeserializeDistributionLayoutInterfaces(glm, faces, recvStream);
		DeserializeDistributionLayoutInterfaces(glm, vols, recvStream);

	//	read the attached data
		if(!mg.has_vertex_attachment(aPosition))
			mg.attach_to_vertices(aPosition);

		for(size_t i = 0; i < mg.num_levels(); ++i){
			if(!DeserializeAttachment<VertexBase>(mg, aPosition,
												mg.begin<VertexBase>(i),
												mg.end<VertexBase>(i),
												recvStream))
				return false;
		}
	}

//	create local new interfaces
	if((int)vVertexLayouts.size() > localProcRank){
	//	we have to create interfaces to the processes that we did send data to.
	//	by setting alias indices of other processes to -1, we achieve that they
	//	are ignored
		vector<int> processMap(highestSendToRank + 1, -1);
		for(size_t i = 0; i < vSendToRanks.size(); ++i){
			processMap[vSendToRanks[i]] = vSendToRanks[i];
		}

		AddHorizontalInterfaces(glm, vVertexLayouts[localProcRank], &processMap);
		AddHorizontalInterfaces(glm, vEdgeLayouts[localProcRank], &processMap);
		AddHorizontalInterfaces(glm, vFaceLayouts[localProcRank], &processMap);
		AddHorizontalInterfaces(glm, vVolumeLayouts[localProcRank], &processMap);
	}

//	update the layouts - if required.
	//if(bErasedElementsFromLayout)
	{
		distGridMgr.grid_layouts_changed(false);
	}

//	only for debugging: log the layout-map-details
	UG_LOG("-- VertexBase\n");
	pcl::LogLayoutMapStructure<VertexBase>(glm);
	UG_LOG("-- EdgeBase\n");
	pcl::LogLayoutMapStructure<EdgeBase>(glm);
	UG_LOG("-- Face\n");
	pcl::LogLayoutMapStructure<Face>(glm);
	UG_LOG("-- Volume\n");
	pcl::LogLayoutMapStructure<Volume>(glm);

//	clean up
	mg.detach_from_vertices(aInt);
	mg.detach_from_edges(aInt);
	mg.detach_from_faces(aInt);
	mg.detach_from_volumes(aInt);

	return true;
}

}//	end of namespace
