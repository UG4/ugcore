// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 27.01.2011 (m,d,y)
 
#include <vector>
#include "parallelization_util.h"
#include "communication_policies.h"

using namespace std;

namespace ug{

int BuildOneToManyLayout(IndexLayout& masterLayoutOut,
						  IndexLayout& slaveLayoutOut,
						  int rootProcID,
						  IndexLayout& masterLayout,
						  IndexLayout& slaveLayout,
						  int highestReferencedIndex,
						  pcl::ProcessCommunicator procComm,
						  std::vector<int>* pNewMasterIDsOut)
{
//	todo: some of the gather calls can be combined.

	//int localRank = procComm.get_local_proc_id();
	vector<IndexLayout::Element> oldMasterNodes;
	vector<IndexLayout::Element> oldSlaveNodes;
	pcl::ParallelCommunicator<IndexLayout> interfaceComm;
	int algVecSize = highestReferencedIndex + 1;
	if(algVecSize < 0) algVecSize = 0;

////////////////////////
//	we have to inform rootProc where old master interfaces lie
	if(!masterLayout.empty()){
	//	collect all master nodes
		CollectUniqueElements(oldMasterNodes, masterLayout);
	}

	int numLocalMasterNodes = (int)oldMasterNodes.size();
	std::vector<int> numOldMasterNodes(procComm.size(), 0);
	procComm.gather(&numLocalMasterNodes, 1, PCL_DT_INT,
					&numOldMasterNodes.front(),
					1, PCL_DT_INT, rootProcID);

////////////////////////
//	each process that has local master nodes creates a slave interface to root proc
	if(!oldMasterNodes.empty()){
		IndexLayout::Interface& interface = slaveLayoutOut.interface(rootProcID);
		for(size_t i = 0; i < oldMasterNodes.size(); ++i){
			interface.push_back(oldMasterNodes[i]);
		}
	}

////////////////////////
//	rootproc now creates a master interface to all old master interfaces
	int nodeCounter = 0;
//	a vector that holds a unique id for each new node on the source process
	vector<int> rootIDs;

	if(pcl::GetProcRank() == rootProcID){

		for(size_t i = 0; i < procComm.size(); ++i){
			int numEntries = numOldMasterNodes[i];
			if(numEntries > 0){
			//	get the interface to the associated process
				int procID = procComm.get_proc_id(i);
				IndexLayout::Interface& interface = masterLayoutOut.interface(procID);
				for(int j = 0; j < numEntries; ++j){
					rootIDs.push_back(nodeCounter);
					interface.push_back(nodeCounter++);
				}
			}
		}

	}

////////////////////////
//	send the local ids of the new nodes on the master proc to
//	associated slaves (old masters)
	if(pcl::GetProcRank() == rootProcID){
	//	use the new interfaces to copy the localID of each node from
	//	rootProc to its slaves (the old masters).
		ComPol_VecCopy<vector<int> > compolCopy(&rootIDs);
		interfaceComm.send_data(masterLayoutOut, compolCopy);
	}

//	receive at all processes containing old master nodes.
	vector<int> masterIDs(algVecSize, -1);
	ComPol_VecCopy<vector<int> > compolCopy(&masterIDs);

	if(!oldMasterNodes.empty()){
		interfaceComm.receive_data(slaveLayoutOut, compolCopy);
	}

	interfaceComm.communicate();

////////////////////////
//	forward those masterIDs to associated slaves (old slaves)
	interfaceComm.send_data(masterLayout, compolCopy);
	interfaceComm.receive_data(slaveLayout, compolCopy);
	interfaceComm.communicate();

////////////////////////
//	now we have to build interfaces between root-proc and the old slaves.
//	inform the root process where slave interfaces lie
	if(!slaveLayout.empty()){
	//	collect all master nodes
		CollectUniqueElements(oldSlaveNodes, slaveLayout);
	}

	int numLocalSlaveNodes = (int)oldSlaveNodes.size();
	std::vector<int> numOldSlaveNodes(procComm.size(), 0);
	procComm.gather(&numLocalSlaveNodes, 1, PCL_DT_INT,
					&numOldSlaveNodes.front(),
					1, PCL_DT_INT, rootProcID);

////////////////////////
//	now that we know how many old slave nodes each process contains,
//	we can collect the localIDs associated with those slaves.

//	prepare the send vecs
	vector<int> masterIDsOfSlaves(oldSlaveNodes.size());
	for(size_t i = 0; i < oldSlaveNodes.size(); ++i){
		masterIDsOfSlaves[i] = masterIDs[oldSlaveNodes[i]];
	}

//	now gather this data on the root process
	vector<int> associatedSlaveIDs;
	vector<int> arrayDisplacements;
	size_t totalSlaveDataSize = 0;

	if(pcl::GetProcRank() == rootProcID){
	//	first calculate the total size of data that we have to receive
	//	and the displacements for each data-buffer
		arrayDisplacements.resize(numOldSlaveNodes.size());
		for(size_t i = 0; i < numOldSlaveNodes.size(); ++i){
			arrayDisplacements[i] = totalSlaveDataSize;
			totalSlaveDataSize += numOldSlaveNodes[i];
		}

		associatedSlaveIDs.resize(totalSlaveDataSize);
	}

	procComm.gatherv(&masterIDsOfSlaves.front(), (int)masterIDsOfSlaves.size(),
					 PCL_DT_INT, &associatedSlaveIDs.front(),
					 &numOldSlaveNodes.front(),
					 &arrayDisplacements.front(),
					 PCL_DT_INT, rootProcID);

////////////////////////
//	finally we build new slave interfaces to root for old slave interfaces
//	and corresponding new master interfaces on root.
	if(oldSlaveNodes.size() > 0){
		IndexLayout::Interface& interface = slaveLayoutOut.interface(rootProcID);
		for(size_t i = 0; i < oldSlaveNodes.size(); ++i){
			interface.push_back(oldSlaveNodes[i]);
		}
	}

	if(pcl::GetProcRank() == rootProcID){
		size_t counter = 0;
		for(size_t i = 0; i < numOldSlaveNodes.size(); ++i){
			int numEntries = numOldSlaveNodes[i];
			if(numEntries > 0){
				int procID = procComm.get_proc_id(i);
				IndexLayout::Interface& interface = masterLayoutOut.interface(procID);
				for(int j = 0; j < numEntries; ++j){
					interface.push_back(associatedSlaveIDs[counter++]);
				}
			}
		}
	}

////////////////////////
//	CODE BELOW (except return statement) only for debugging.
/*
//	send 1 from new masters to new slaves
	vector<int> recBuff(algVecSize, 0);
	if(pcl::GetProcRank() == rootProcID){
		vector<int> sendBuff(nodeCounter, 1);
		ComPol_VecCopy<vector<int> > compolCopy(&sendBuff);
		interfaceComm.send_data(masterLayoutOut, compolCopy);
	}

	{
		ComPol_VecCopy<vector<int> > compolCopy(&recBuff);
		interfaceComm.receive_data(slaveLayoutOut, compolCopy);
		interfaceComm.communicate();

		UG_LOG("test-communication results: ");
		for(size_t i = 0; i < recBuff.size(); ++i){
			UG_LOG(recBuff[i]);
		}
		UG_LOG(endl);
	}
*/
/*
	UG_LOG("received masterIDs: ");
	for(size_t i = 0; i < associatedSlaveIDs.size(); ++i){
		UG_LOG(associatedSlaveIDs[i] << " ");
	}
	UG_LOG(endl);
*/

/*
	UG_LOG("num slave nodes: ")
	for(size_t i = 0; i < numOldSlaveNodes.size(); ++i){
		UG_LOG(" " << numOldSlaveNodes[i]);
	}
	UG_LOG(endl);
*/
/*
	UG_LOG("One to Many master layout:\n");
	LogIndexLayout(masterLayoutOut);
	UG_LOG(endl);
	UG_LOG("One to Many slave layout:\n");
	LogIndexLayout(slaveLayoutOut);
*/

//	UG_LOG("BuildOneToManyLayout done.\n");

//	nodeCounter != 0 only on rootProc
//todo - instead of copying, the algorithm could directly work on pNewMasterIDsOut.
	if(pNewMasterIDsOut)
		*pNewMasterIDsOut = masterIDs;
	return nodeCounter;
}

}// end of namespace
