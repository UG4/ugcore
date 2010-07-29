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
/**	This method is only used if the source-grid shall be kept on
 *	the distributing process. Instead of creating horizontal interfaces
 *	it thus has to create virtual horizontal interfaces.
 *	\todo: create virtual horizontal interfaces only for surface elements.
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
				assert((int)pProcessMap->size() > procID && "process-map to small.");
				if((int)pProcessMap->size() > procID)
					procID = (*pProcessMap)[procID];
			}
			
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
 *	problematic here. Make sure to avoid auto-element-creation.*/
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
	
//	iterate over the levels of the multi-grid
	for(size_t level = 0; level < mg.num_levels(); ++level){
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
	//ISubsetHandler& sh = shInOut;
	//GridLayoutMap& glm = gridLayoutMapInOut;
	
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

//TODO: Only prepare layouts for processes to which data is actually send.
	CreateDistributionLayouts(vVertexLayouts, vEdgeLayouts, vFaceLayouts,
							  vVolumeLayouts, mg, shPartition,
							  false, &msel);
							  
	
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
	for(size_t si = 0; si < shPartition.num_subsets(); ++si)
	{
		//	todo: we have to check whether interfaces change and notify the
		//		  associated processes.
		if((int)si != localProcRank){
		//	check whether there is anything in the partition
			if(!shPartition.empty(si)){
				vSendToRanks.push_back((size_t)si);
			//	log
				UG_LOG("  sending " << shPartition.num<Face>(si)
						<< " faces to process " << si << endl);
			}
		}
	}

//	send and receive communication requests from other processes.
	vector<int> vReceiveFromRanks;
	CommunicateInvolvedProcesses(vReceiveFromRanks, vSendToRanks, procComm);

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


//	copy and delete the grid-parts
//	todo: interfaces are ignored in the moment
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
		SerializeSubsetHandler(mg, shInOut,
							   msel.get_geometric_object_collection(),
							   sendStream);
							   
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
	
//	all send data is now ready. We can now erase unrequired parts of the local grid.
//	We only have to erase elements if some are actually moved to another process.
	if(vSendBlockSizes.size() > 0)
	{
	//	erase all elements that will not remain.
		if((int)vVertexLayouts.size() > localProcRank)
		{
			msel.clear();
			SelectNodesInLayout(msel, vVertexLayouts[localProcRank]);
			SelectNodesInLayout(msel, vEdgeLayouts[localProcRank]);
			SelectNodesInLayout(msel, vFaceLayouts[localProcRank]);
			SelectNodesInLayout(msel, vVolumeLayouts[localProcRank]);
			InvertSelection(msel);
			EraseSelectedObjects(msel);
		}
		else{
		//	if there is no entry in vVertexLayouts, we can remove all elements
			mg.clear_geometry();
		}
	}
UG_LOG(endl);
	
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
	
//	log the sizes of the data-bocks that will be received
	for(size_t i = 0; i < vReceiveFromRanks.size(); ++i)
	{
		UG_LOG("  receiving " << vRecvBlockSizes[i]
				<< " bytes from process " << vReceiveFromRanks[i] << endl);
	}
	
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
//			the receiving processes initially contain empty grids.
//TODO: implement grid-merging
	if(recvStream.size() > 0)
	{
	//	deserialize the multi-grid
		DeserializeMultiGridElements(mg, recvStream);

	//	deserialize subset handler
		if(!DeserializeSubsetHandler(mg, shInOut,
									mg.get_geometric_object_collection(),
									recvStream))
			return false;
		
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
	
//	clean up
	mg.detach_from_vertices(aInt);
	mg.detach_from_edges(aInt);
	mg.detach_from_faces(aInt);
	mg.detach_from_volumes(aInt);
	
	return true;
}

}//	end of namespace
