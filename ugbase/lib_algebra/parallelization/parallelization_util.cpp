// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 27.01.2011 (m,d,y)
 
#include <vector>
#include "parallelization_util.h"
#include "communication_policies.h"

using namespace std;
using namespace pcl;

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

static void CopyInterfaceEntrysToDomainDecompositionLayouts(
		IndexLayout& subdomLayoutOut, IndexLayout& processLayoutOut,
		IndexLayout& deltaNbrLayoutOut, IndexLayout& crossPointLayoutOut,
		IndexLayout& standardLayout, vector<int>& flags,
		IDomainDecompositionInfo& ddinfo)
{
	typedef IndexLayout::Interface	Interface;
	typedef IndexLayout::iterator	InterfaceIter;
	typedef Interface::Element		Element;
	typedef Interface::iterator		ElemIter;

//	the local process and subdomain id
	int localProcID = pcl::GetProcRank();
	int localSubdomID = ddinfo.map_proc_id_to_subdomain_id(localProcID);

	for(InterfaceIter iiter = standardLayout.begin();
		iiter != standardLayout.end(); ++iiter)
	{
	//	connected proc and subdomain ids
		int connProc = standardLayout.proc_id(iiter);
		int connSubdomID = ddinfo.map_proc_id_to_subdomain_id(connProc);
		Interface& stdInterface = standardLayout.interface(iiter);

	//	in order to avoid creation of unrequired interfaces, we first count
	//	how many entries of each flag-type (-2, -1, >= 0) exist in the current
	//	standard interface.
		int numMult = 0;	// flag: -2
		int numOne = 0;		// flag: >= 0
		int numNone = 0;	// flag: -1
		for(ElemIter eiter = stdInterface.begin();
			eiter != stdInterface.end(); ++eiter)
		{
			int flag = flags[stdInterface.get_element(eiter)];
			if(flag == -2)
				++numMult;
			else if(flag == -1)
				++numNone;
			else
				++numOne;
		}

		if(numOne){
		//	we have to take care of entries which are connected to exactly
		//	one other subdomain. Those can either be entries lying in
		//	an interface to another subdomain, or entries which
		//	lie in an interface to a process in the same subdomain, which
		//	also have neighbours in another subdomain.
			if(localSubdomID != connSubdomID){
			//	the interface connects two subdomains. all entries which are
			//	connected to exactly one other subdomain have to be copied
			//	to subdomLayoutOut.
				Interface& subdomInterface = subdomLayoutOut.interface(connProc);
				for(ElemIter eiter = stdInterface.begin();
					eiter != stdInterface.end(); ++eiter)
				{
					Element elem = stdInterface.get_element(eiter);
					if(flags[elem] >= 0)
						subdomInterface.push_back(elem);
				}
			}
			else{
			//	delta-neighbours lie in two interfaces - the processInterface
			//	and the deltaNbrInterface
				Interface& deltaNbrInterface = deltaNbrLayoutOut.interface(connProc);
				Interface& processInterface = processLayoutOut.interface(connProc);
				for(ElemIter eiter = stdInterface.begin();
					eiter != stdInterface.end(); ++eiter)
				{
					Element elem = stdInterface.get_element(eiter);
					if(flags[elem] >= 0){
						deltaNbrInterface.push_back(elem);
						processInterface.push_back(elem);
					}
				}
			}
		}

		if(numMult){
		//	There are elements in the interface which are connected to
		//	multiple subdomains. Those elements are thus regarded as
		//	cross points.
			Interface& crossInterface = crossPointLayoutOut.interface(connProc);

		//	iterate over the elements again and assign them to their interfaces.
			for(ElemIter eiter = stdInterface.begin();
				eiter != stdInterface.end(); ++eiter)
			{
				Element elem = stdInterface.get_element(eiter);
				if(flags[elem] == -2)
					crossInterface.push_back(elem);
			}
		}

		if(numNone){
		//	now we copy all entries which are only connected to a process in the
		//	same subdomain into the processInterfaces.
			Interface& processInterface = processLayoutOut.interface(connProc);

		//	iterate over the elements again and assign them to their interfaces.
			for(ElemIter eiter = stdInterface.begin();
				eiter != stdInterface.end(); ++eiter)
			{
				Element elem = stdInterface.get_element(eiter);
				if(flags[elem] == -1)
					processInterface.push_back(elem);
			}
		}
	}
}

void BuildDomainDecompositionLayouts(
		IndexLayout& subdomMastersOut, IndexLayout& subdomSlavesOut,
		IndexLayout& processMastersOut, IndexLayout& processSlavesOut,
		IndexLayout& deltaNbrMastersOut, IndexLayout& deltaNbrSlavesOut,
		IndexLayout& crossPointMastersOut, IndexLayout& crossPointSlavesOut,
		IndexLayout& standardMasters, IndexLayout& standardSlaves,
		int highestReferencedIndex, IDomainDecompositionInfo& ddinfo)
{
//	some typedefs
	typedef IndexLayout::Interface	Interface;
	typedef IndexLayout::iterator	InterfaceIter;
	typedef Interface::iterator		ElemIter;

//	the local process and subdomain id
	int localProcID = pcl::GetProcRank();
	int localSubdomID = ddinfo.map_proc_id_to_subdomain_id(localProcID);

//	the size which the flags vector has to have
	int flagVecSize = highestReferencedIndex + 1;
	if(flagVecSize < 0)
		flagVecSize = 0;

	pcl::ParallelCommunicator<IndexLayout> interfaceComm;

////////////////////////
//	this vector will contain an integer
//	* >= 0 for each element which has a copy in exactly one other subdomain,
//	* -2 for each element which has a copy in more than one other subdomain,
//	* -1 if the element is not connected to another subdomain.

//	delta interface and 0 for all others.
	vector<int> flags(flagVecSize, -1);

//	iterate over all standard master interfaces and fill the vector as
//	described above for master elements
	for(InterfaceIter iiter = standardMasters.begin();
		iiter != standardMasters.end(); ++iiter)
	{
	//	connected proc and subdomain ids
		int connProc = standardMasters.proc_id(iiter);
		int connSubdomID = ddinfo.map_proc_id_to_subdomain_id(connProc);
		if(localSubdomID != connSubdomID){
		//	the interface connects two subdomains.
			Interface& interface = standardMasters.interface(iiter);
			for(ElemIter eiter = interface.begin();
				eiter != interface.end(); ++eiter)
			{
				Interface::Element ind = interface.get_element(eiter);
				if(flags[ind] == -1)	// not connection was yet known.
					flags[ind] = connSubdomID;
				else if(flags[ind] != connSubdomID)
					flags[ind] = -2;	//	at least one other connection is known.
			}
		}
	}

//	communicate the flags to all standard slaves
	ComPol_VecCopy<vector<int> > compolCopy(&flags);
	interfaceComm.send_data(standardMasters, compolCopy);
	interfaceComm.receive_data(standardSlaves, compolCopy);
	interfaceComm.communicate();

////////////////////////
//	we can now build the master and slave interfaces
	CopyInterfaceEntrysToDomainDecompositionLayouts(
			subdomMastersOut, processMastersOut, deltaNbrMastersOut,
			crossPointMastersOut, standardMasters, flags, ddinfo);

	CopyInterfaceEntrysToDomainDecompositionLayouts(
			subdomSlavesOut, processSlavesOut, deltaNbrSlavesOut,
			crossPointSlavesOut, standardSlaves, flags, ddinfo);

//	done.
//DEBUG ONLY
/*
//	perform a test - communicate flags between all processes over the new layouts
	interfaceComm.send_data(subdomMastersOut, compolCopy);
	interfaceComm.send_data(processMastersOut, compolCopy);
	interfaceComm.send_data(deltaNbrMastersOut, compolCopy);
	interfaceComm.receive_data(subdomSlavesOut, compolCopy);
	interfaceComm.receive_data(processSlavesOut, compolCopy);
	interfaceComm.receive_data(deltaNbrSlavesOut, compolCopy);
	interfaceComm.communicate();
*/
}

}// end of namespace
