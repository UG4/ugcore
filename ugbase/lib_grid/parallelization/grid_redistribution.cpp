// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.03.2011 (m,d,y)
 
#include <vector>
#include "grid_distribution.h"
#include "util/distribution_util.h"
#include "parallelization_util.h"

using namespace std;

namespace ug{

template <class TGeomObj, class TAATargetProcs, class TAAGlobalIDs>
static
void LogTargetProcs(MultiGrid& mg, TAATargetProcs& aaTargetProcs,
					  TAAGlobalIDs& aaIDs)
{
	typedef typename geometry_traits<TGeomObj>::iterator GeomObjIter;
	UG_LOG("  TARGET-PROCS:\n");
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		UG_LOG("    level " << lvl << ":\n");
		for(GeomObjIter iter = mg.begin<VertexBase>(lvl);
			iter != mg.end<VertexBase>(lvl); ++iter)
		{
		//	log the id
			TGeomObj* o = *iter;
			UG_LOG("      " << aaIDs[o].first << "_" << aaIDs[o].second << ":\t");

		//	log the target procs
			vector<int>& v = aaTargetProcs[o];
			for(size_t i = 0; i < v.size(); ++i){
				UG_LOG(" " << v[i]);
			}

			UG_LOG(endl);
		}
	}
}

template <class TGeomObj, class TAATransferInfos, class TAAGlobalIDs>
static
void LogTransferInfos(MultiGrid& mg, TAATransferInfos& aaTransferInfos,
					  TAAGlobalIDs& aaIDs)
{
	typedef typename geometry_traits<TGeomObj>::iterator GeomObjIter;
	UG_LOG("  TRANSFER-INFOS:\n");
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		UG_LOG("    level " << lvl << ":\n");
		for(GeomObjIter iter = mg.begin<VertexBase>(lvl);
			iter != mg.end<VertexBase>(lvl); ++iter)
		{
		//	log the id
			TGeomObj* o = *iter;
			UG_LOG("      " << aaIDs[o].first << "_" << aaIDs[o].second << ":\t");

		//	log the target procs
			vector<RedistributionNodeTransferInfo>& infos = aaTransferInfos[o];
			for(size_t i = 0; i < infos.size(); ++i){
				UG_LOG("src: " << infos[i].srcProc);
				UG_LOG(", target: " << infos[i].targetProc);
				UG_LOG(", bMove: " << infos[i].bMove << "  ||  ");
			}

			UG_LOG(endl);
		}
	}
}

bool RedistributeGrid(DistributedGridManager& distGridMgrInOut,
					  SubsetHandler& shPartition,
					  const GridDataSerializer& serializer,
					  const GridDataDeserializer& deserializer,
					  int* processMap,
					  const pcl::ProcessCommunicator& procComm)
{

//	Algorithm outline:
//	- Distribute global IDs (if they were not already computed),
//	- Prepare redistribution groups
//		* This involves immediate updates to associated interfaces:
//		  Each process has to inform neighbor processes to where elements
//		  are sent, which lie in common interfaces.
//		  Note that several migration forms are possible: move, copy
//		* GlobalIDs have to be associated with each element.
//		* Make sure not to specify interfaces to the process to which the data is sent.
//	- CommunicateInvolvedProcesses
//	- Distribute the redistribution groups
//		* This involves serializing the associated grid-parts
//		* At this point one would also want to serialize custom user-data.
//	- Update local interfaces.
//		* This can be safely done here, since all redistribution groups are complete
//		  and since all neighbors are informed which interface elements go where.
//		* Note that this is not the final update. We're only erasing entries which were
//		  sent to other processes or entries in interfaces where neighbor processes
//		  informed us that the entries will be moved to other processes.
//		  We will also create new entries for those moved neighbor elements - as long
//		  as we didn't move our elements to another process, too.
//		* If connected entries on a neighbor proc A are moved to another process B,
//		  note that if an interface IAB already exists, there may also already exist some
//		  connections. Don't readd those!
//	- Erase all parts of the grid which are no longer used.
//		* Note that associated entries are already removed from local interfaces.
//	- Create new elements and add new interface entries.
//		* Here we use the global IDs to check whether an element already exists
//		  on the given process.
//		* The check can e.g. be implemented using hashes. Make sure to add all new
//		  elements to that hash.
//		* Create temporary arrays into which the referenced elements are sorted.
//		  This will make it easy to add the new interface entries.
//		* Note that if a referenced element already existed, it is also possible,
//		  that associated interface entries already exist. This has to be checked
//		  before new interface entries are added.
//

//todo: currently the local proc also serializes and deserializes all data
//		which he sends to itself. This involves clearing the whole grid
//		and recreating it from the serialized data.
//		It was implemented this way since it simply required less code.
//		In the future one may consider changing this behavior, so that
//		parts which stay on the grid will simply be maintained.
//		Affected parts of the algorithm are
//		* COMMUNICATE INVOLVED PROCESSES
//		* INTERMEDIATE CLEANUP
//		Note that one has to carefully adjust interfaces to the ones given
//		in the local redistribution map.

	if(!distGridMgrInOut.get_assigned_grid())
		return false;

	//int localProcRank = pcl::GetProcRank();

	MultiGrid& mg = *distGridMgrInOut.get_assigned_grid();
	GridLayoutMap& glm = distGridMgrInOut.grid_layout_map();

//	The selector will be of frequent use to speed up some algorithms
	MGSelector msel(mg);

//	Since we will change huge parts of the underlying grid and the grid-layout-map,
//	we'll disable auto-insertion of elements in the distributed-grid-manager.
//	This means we carefully have to take care of all interface changes.
	distGridMgrInOut.enable_ordered_element_insertion(false);

////////////////////////////////
//	GLOBAL IDS
//todo:	only create global ids if they aren't already present
	CreateAndDistributeGlobalIDs<VertexBase>(mg, glm);
	CreateAndDistributeGlobalIDs<EdgeBase>(mg, glm);
	CreateAndDistributeGlobalIDs<Face>(mg, glm);
	CreateAndDistributeGlobalIDs<Volume>(mg, glm);

	Grid::AttachmentAccessor<VertexBase, AGeomObjID> aaIDVRT(mg, aGeomObjID);
	Grid::AttachmentAccessor<EdgeBase, AGeomObjID> 	aaIDEDGE(mg, aGeomObjID);
	Grid::AttachmentAccessor<Face, AGeomObjID> 		aaIDFACE(mg, aGeomObjID);
	Grid::AttachmentAccessor<Volume, AGeomObjID> 	aaIDVOL(mg, aGeomObjID);

////////////////////////////////
//	REDISTRIBITION LAYOUTS
//	we have to store the layouts for all the processes.
	vector<RedistributionVertexLayout> vertexLayouts;
	vector<RedistributionEdgeLayout> edgeLayouts;
	vector<RedistributionFaceLayout> faceLayouts;
	vector<RedistributionVolumeLayout> volumeLayouts;

//	create the redistribution layouts
//	we want to access the created target-procs and transfer-info-vec attachments
//	later on.
	typedef Attachment<vector<int> > AIntVec;
	AIntVec aTargetProcs;
	ARedistributionNodeTransferInfoVec aTransferInfos;

	Grid::AttachmentAccessor<VertexBase, AIntVec>
		aaTargetProcsVRT(mg, aTargetProcs, true);
	Grid::AttachmentAccessor<EdgeBase, AIntVec>
		aaTargetProcsEDGE(mg, aTargetProcs, true);
	Grid::AttachmentAccessor<Face, AIntVec>
		aaTargetProcsFACE(mg, aTargetProcs, true);
	Grid::AttachmentAccessor<Volume, AIntVec>
		aaTargetProcsVOL(mg, aTargetProcs, true);

	Grid::AttachmentAccessor<VertexBase, ARedistributionNodeTransferInfoVec>
		aaTransferInfosVRT(mg, aTransferInfos, true);
	Grid::AttachmentAccessor<EdgeBase, ARedistributionNodeTransferInfoVec>
		aaTransferInfosEDGE(mg, aTransferInfos, true);
	Grid::AttachmentAccessor<Face, ARedistributionNodeTransferInfoVec>
		aaTransferInfosFACE(mg, aTransferInfos, true);
	Grid::AttachmentAccessor<Volume, ARedistributionNodeTransferInfoVec>
		aaTransferInfosVOL(mg, aTransferInfos, true);

//	Now create the redistribution layouts. Note that this involves
//	some communication in order to synchronize interfaces between
//	connected procs.
	CreateRedistributionLayouts(vertexLayouts, edgeLayouts, faceLayouts,
					volumeLayouts, distGridMgrInOut, shPartition, false,
					&msel, processMap, &aTargetProcs, &aTransferInfos);



//BEGIN - ONLY FOR DEBUG
	LogTargetProcs<VertexBase>(mg, aaTargetProcsVRT, aaIDVRT);


	LogTransferInfos<VertexBase>(mg, aaTransferInfosVRT, aaIDVRT);

	TestRedistributionLayouts(vertexLayouts);
//END - ONLY FOR DEBUG



////////////////////////////////
//	COMMUNICATE INVOLVED PROCESSES
//	each process has to know with which other processes it
//	has to communicate.
	vector<int> sendToRanks, recvFromRanks, sendPartitionInds;

//	for each subset which is not emtpy we'll have to send data to
//	the associated process.
	for(int si = 0; si < shPartition.num_subsets(); ++si){
		if(!shPartition.empty(si)){
			int toProc = si;
		//	if a process map exists, we'll use the associated process
			if(processMap)
				toProc = processMap[si];

			sendToRanks.push_back(toProc);
			sendPartitionInds.push_back(toProc);
		}
	}

	pcl::CommunicateInvolvedProcesses(recvFromRanks, sendToRanks, procComm);


//BEGIN - ONLY FOR DEBUG
	UG_LOG("  send to ranks:");
	for(size_t i = 0; i < sendToRanks.size(); ++i){
		UG_LOG(" " << sendToRanks[i]);
	}
	UG_LOG(endl);

	UG_LOG("  receive from ranks:");
	for(size_t i = 0; i < recvFromRanks.size(); ++i){
		UG_LOG(" " << recvFromRanks[i]);
	}
	UG_LOG(endl);
//END - ONLY FOR DEBUG



////////////////////////////////
//	SERIALIZATION
//	For each send-to-process (except the local process itself) we'll fill
//	a binary stream with the data that shall be sent.
//	a temporary int-attachment is required at all geometric objects
	AInt aLocalInd;
	mg.attach_to_all(aLocalInd);

//	out and sendSegSizes will be used to distribute the grid.
	BinaryStream out;
	vector<int> outSegSizes;

	for(size_t i_to = 0; i_to < sendToRanks.size(); ++i_to){
		int partInd = sendPartitionInds[i_to];

	//	First we'll serialize the global ids of all distributed elements
		Serialize(out, vertexLayouts[partInd].m_globalIDs);
		Serialize(out, edgeLayouts[partInd].m_globalIDs);
		Serialize(out, faceLayouts[partInd].m_globalIDs);
		Serialize(out, volumeLayouts[partInd].m_globalIDs);

	//	the last size is required to calculate the size of the new segment
		size_t oldSize = out.size();

	//	Now let's serialize the grid and the redistribution layout interfaces
		SerializeGridAndDistributionLayouts(out, mg, vertexLayouts[partInd],
				edgeLayouts[partInd], faceLayouts[partInd], volumeLayouts[partInd],
				aLocalInd, aLocalInd, aLocalInd, aLocalInd, &msel);

	//	next thing is to serialize data associated with the layouts nodes
		serializer.serialize(out, vertexLayouts[partInd].node_vec().begin(),
							 vertexLayouts[partInd].node_vec().end());
		serializer.serialize(out, edgeLayouts[partInd].node_vec().begin(),
							 edgeLayouts[partInd].node_vec().end());
		serializer.serialize(out, faceLayouts[partInd].node_vec().begin(),
							 faceLayouts[partInd].node_vec().end());
		serializer.serialize(out, volumeLayouts[partInd].node_vec().begin(),
							 volumeLayouts[partInd].node_vec().end());

	//	size of the segment we just wrote to out
		outSegSizes.push_back((int)(out.size() - oldSize));
	}

//	now distribute the packs between involved processes
	BinaryStream in;
	vector<int> inSegSizes(recvFromRanks.size());

	procComm.distribute_data(in, &inSegSizes.front(),
							&recvFromRanks.front(), (int)recvFromRanks.size(),
							out.buffer(), &outSegSizes.front(),
							&sendToRanks.front(), (int)sendToRanks.size());


////////////////////////////////
//	INTERMEDIATE CLEANUP
//	we'll erase everything and deserialize even the local grid from stream
	mg.clear_geometry();
	glm.clear();

/*
	{
	//	select all elements which shall be cleared.
	//	we'll do this by first selecting all elements on the local proc
	//	and afterwards inverting the selection.
		SelectNodesInLayout(msel, vertexLayouts[localPartSI]);
		SelectNodesInLayout(msel, edgeLayouts[localPartSI]);
		SelectNodesInLayout(msel, faceLayouts[localPartSI]);
		SelectNodesInLayout(msel, volumeLayouts[localPartSI]);
		InvertSelection(msel);

	//	erase the grid parts and all selected interface entries
		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
			mg.erase(msel.vertices_begin(lvl), msel.vertices_end(lvl));
			mg.erase(msel.edges_begin(lvl), msel.edges_end(lvl));
			mg.erase(msel.faces_begin(lvl), msel.faces_end(lvl));
			mg.erase(msel.volumes_begin(lvl), msel.volumes_end(lvl));
		}
	}
*/
////////////////////////////////
//	DESERILAIZATION
//	before deserializing the received data, we'll sort the data we received
//	by the ranks of the processes from which we received data. This is
//	important to make sure that elements on different processes are added
//	to the interfaces in the right order.
//	vector of pairs <procRank, readPos>
	std::vector<pair<int, size_t> >	recvProcsSorted;
	size_t readPos = 0;
	for(size_t i = 0; i < recvFromRanks.size(); ++i){
		recvProcsSorted.push_back(make_pair(recvFromRanks[i], readPos));
		readPos += inSegSizes[i];
	}

	sort(recvProcsSorted.begin(), recvProcsSorted.end());


//BEGIN - ONLY FOR DEBUG
	UG_LOG("RecProcsSorted:");
	for(size_t i = 0; i < recvProcsSorted.size(); ++i){
		UG_LOG(" (" << recvProcsSorted[i].first << ", " << recvProcsSorted[i].second << ")");
	}
	UG_LOG(endl);
//END - ONLY FOR DEBUG


//	deserialize the redist-layouts
//	read redistribution layouts from the stream

//	this multi-grid will be used for deserialization.
//	Note that it uses GRIDOPT_NONE for maximal performance.
//	We'll attach pointers to the associated elements in the target grid.
	MultiGrid tmg(GRIDOPT_NONE);
	AVertexBase aVrt;
	AEdgeBase aEdge;
	AFace aFace;
	AVolume aVol;
	Grid::AttachmentAccessor<VertexBase, AVertexBase> aaVrt(tmg, aVrt, true);
	Grid::AttachmentAccessor<EdgeBase, AEdgeBase> aaEdge(tmg, aEdge, true);
	Grid::AttachmentAccessor<Face, AFace> aaFace(tmg, aFace, true);
	Grid::AttachmentAccessor<Volume, AVolume> aaVol(tmg, aVol, true);

	for(size_t i = 0; i < recvProcsSorted.size(); ++i){
		//int recvFrom = recvProcsSorted[i].first;

	//	set the read position of the in-stream to the data section
	//	of the process we're currently processing.
		in.reset();
		in.read_jump(recvProcsSorted[i].second);

	//	those layouts will be used to deserialize the received data
		RedistributionVertexLayout vrtLayout;
		RedistributionEdgeLayout edgeLayout;
		RedistributionFaceLayout faceLayout;
		RedistributionVolumeLayout volLayout;

	//	First we'll deserialize the global ids of all distributed elements
		Deserialize(in, vrtLayout.m_globalIDs);
		Deserialize(in, edgeLayout.m_globalIDs);
		Deserialize(in, faceLayout.m_globalIDs);
		Deserialize(in, volLayout.m_globalIDs);



		DeserializeMultiGridElements(mg, in, &vrtLayout.node_vec(),
							&edgeLayout.node_vec(), &faceLayout.node_vec(),
							&volLayout.node_vec());
//todo: use the code below instead of the one above
/*
		tmg.clear_geometry();
		DeserializeMultiGridElements(tmg, in, &vrtLayout.node_vec(),
							&edgeLayout.node_vec(), &faceLayout.node_vec(),
							&volLayout.node_vec());
*/
		DeserializeDistributionLayoutInterfaces(vrtLayout, in);
		DeserializeDistributionLayoutInterfaces(edgeLayout, in);
		DeserializeDistributionLayoutInterfaces(faceLayout, in);
		DeserializeDistributionLayoutInterfaces(volLayout, in);

	//	before we'll deserialize the associated data, we'll create the new
	//	elements in the target grid.

	//todo: create elems in mg and check for existing nodes with the same global id.
	//		Also update the layout-vecs so that they point to the elements in mg.

	//	now deserialize the data associated with the elements in the layouts node-vecs
		deserializer.deserialize(in, vrtLayout.node_vec().begin(),
							 	 vrtLayout.node_vec().end());
		deserializer.deserialize(in, edgeLayout.node_vec().begin(),
							 	 edgeLayout.node_vec().end());
		deserializer.deserialize(in, faceLayout.node_vec().begin(),
							 	 faceLayout.node_vec().end());
		deserializer.deserialize(in, volLayout.node_vec().begin(),
							 	 volLayout.node_vec().end());

//BEGIN - ONLY FOR DEBUG
		UG_LOG("received global vertex ids:");
		for(size_t i = 0; i < vrtLayout.m_globalIDs.size(); ++i){
			UG_LOG(" (" << vrtLayout.m_globalIDs[i].first << ", " << vrtLayout.m_globalIDs[i].second << ")");
		}
		UG_LOG(endl);
//END - ONLY FOR DEBUG
	}


////////////////////////////////
//	UPDATE THE DISTRIBUTED GRID MANAGER
	distGridMgrInOut.enable_ordered_element_insertion(true);
	distGridMgrInOut.grid_layouts_changed(false);

////////////////////////////////
//	CLEAN UP
	mg.detach_from_all(aLocalInd);
	mg.detach_from_all(aTargetProcs);
	mg.detach_from_all(aTransferInfos);

	return true;
}

}// end of namespace
