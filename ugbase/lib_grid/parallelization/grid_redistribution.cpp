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
					  ISubsetHandler& shInOut,
					  SubsetHandler& shPartition,
					  int* processMap)
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

	if(!distGridMgrInOut.get_assigned_grid())
		return false;

	MultiGrid& mg = *distGridMgrInOut.get_assigned_grid();
	GridLayoutMap& glm = distGridMgrInOut.grid_layout_map();

//	The selector will be of frequent use to speed up some algorithms
	MGSelector msel(mg);

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

	CreateRedistributionLayouts(vertexLayouts, edgeLayouts, faceLayouts,
					volumeLayouts, distGridMgrInOut, shPartition, false,
					&msel, processMap, &aTargetProcs, &aTransferInfos);

//BEGIN - ONLY FOR DEBUG
	LogTargetProcs<VertexBase>(mg, aaTargetProcsVRT, aaIDVRT);


	LogTransferInfos<VertexBase>(mg, aaTransferInfosVRT, aaIDVRT);

	TestRedistributionLayouts(vertexLayouts);

//END - ONLY FOR DEBUG


////////////////////////////////
//	CLEAN UP
	mg.detach_from_all(aTargetProcs);
	mg.detach_from_all(aTransferInfos);

	return true;
}

}// end of namespace
