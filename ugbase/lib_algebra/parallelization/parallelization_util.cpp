// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 27.01.2011 (m,d,y)
//

#include <algorithm>
#include <vector>
#include "parallelization_util.h"
#include "communication_policies.h"
#include "pcl/pcl.h"
#include "common/profiler/profiler.h"

using namespace std;
using namespace pcl;

namespace ug{

template<>
size_t hash_key<AlgebraID>(const AlgebraID& key)
{
	const unsigned long factor = 1000000;
	const unsigned long ind = (unsigned long)key.index_on_master();
	return  factor * (unsigned long)key.master_proc() * ind + ind;
}

std::ostream& operator<<(std::ostream &out, const AlgebraID &ID)
{
  out << "p" << ID.first << "i" << ID.second;
  return out;
}


///	vector<int> has a variable size
template <> struct block_traits<vector<int> >
{
	enum{
		is_static = 0
	};
};

void CommunicateConnections(vector<vector<int> >& connectionsToProcsOut,
							vector<vector<int> >& connectionsToSubDomsOut,
							const IndexLayout& masterLayout,
							const IndexLayout& slaveLayout,
							int highestReferencedIndex, pcl::IDomainDecompositionInfo& ddinfo)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	typedef IndexLayout::Interface 	Interface;
	typedef IndexLayout::const_iterator 	InterfaceIter;
	typedef Interface::const_iterator		ElemIter;
	typedef Interface::Element		Element;

	int connectionVecSize = highestReferencedIndex + 1;
	if(connectionVecSize < 0) connectionVecSize = 0;

	connectionsToProcsOut.clear();
	connectionsToProcsOut.resize(connectionVecSize);

	connectionsToSubDomsOut.clear();
	connectionsToSubDomsOut.resize(connectionVecSize);

	int localProc = pcl::GetProcRank();

//	iterate over all master interfaces
	for(InterfaceIter iiter = masterLayout.begin();
		iiter != masterLayout.end(); ++iiter)
	{
		const Interface& interface = masterLayout.interface(iiter);
		int targetProc = interface.get_target_proc();

	//	iterate over all elements
		for(ElemIter eiter = interface.begin();
			eiter != interface.end(); ++eiter)
		{
			Element elem = interface.get_element(eiter);
			vector<int>& connToProc   = connectionsToProcsOut[elem];
			vector<int>& connToSubDom = connectionsToSubDomsOut[elem];

		//	the first entry in each connection is the master elements process
			if(connToProc.empty())
				connToProc.push_back(localProc);

		//	now push the slave elemens process
			connToProc.push_back(targetProc);


		//	the first entry in each connection is the subdomain id of the master elements process
			if(connToSubDom.empty())
				connToSubDom.push_back(ddinfo.map_proc_id_to_subdomain_id(localProc));

		//	now push the subdomain id of the target proc (i.e. of the slave element)
			connToSubDom.push_back(ddinfo.map_proc_id_to_subdomain_id(targetProc));
		}
	}

//	now communicate the connectionsToProcsOut to the slaves
	pcl::InterfaceCommunicator<IndexLayout> interfaceComm;
	ComPol_VecCopy<vector<vector<int> > > compolCopy(&connectionsToProcsOut);
	interfaceComm.send_data(masterLayout, compolCopy);
	interfaceComm.receive_data(slaveLayout, compolCopy);
	interfaceComm.communicate();

//	... and also the connectionsToSubDomsOut
	compolCopy.set_vector(&connectionsToSubDomsOut);
	interfaceComm.send_data(masterLayout, compolCopy);
	interfaceComm.receive_data(slaveLayout, compolCopy);
	interfaceComm.communicate();
}

int GetHighestReferencedIndex(IndexLayout& layout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	we have to find the highest referenced index.
//	In order to do so, we will iterate over all entries of the
//	given layouts and check each.
	int highestReferencedIndex = -1;
	for(IndexLayout::iterator iiter = layout.begin();
		iiter != layout.end(); ++iiter)
	{
		IndexLayout::Interface& interface = layout.interface(iiter);
		for(IndexLayout::Interface::iterator eiter = interface.begin();
			eiter != interface.end(); ++eiter)
		{
			IndexLayout::Interface::Element elem = interface.get_element(eiter);
			if((int)elem > highestReferencedIndex)
				highestReferencedIndex = (int)elem;
		}
	}

	return highestReferencedIndex;
}

int BuildOneToManyLayout(IndexLayout& masterLayoutOut,
						  IndexLayout& slaveLayoutOut,
						  int rootProcID,
						  IndexLayout& masterLayout,
						  IndexLayout& slaveLayout,
						  pcl::ProcessCommunicator procComm,
						  std::vector<int>* pNewMasterIDsOut)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	todo: some of the gather calls can be combined.

	//int localRank = procComm.get_local_proc_id();
	vector<IndexLayout::Element> oldMasterNodes;
	vector<IndexLayout::Element> oldSlaveNodes;
	pcl::InterfaceCommunicator<IndexLayout> interfaceComm;

	int highestReferencedIndex = max(GetHighestReferencedIndex(masterLayout),
									 GetHighestReferencedIndex(slaveLayout));

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
	if(pNewMasterIDsOut){
		size_t oldSize = pNewMasterIDsOut->size();
		*pNewMasterIDsOut = masterIDs;
		if(oldSize > pNewMasterIDsOut->size())
			pNewMasterIDsOut->resize(oldSize, -1);
	}
	return nodeCounter;
}

/**	Please note that this method only handles elements with flags >=0, -1 and -2.*/
static void CopyInterfaceEntrysToDomainDecompositionLayouts(
		IndexLayout& subdomLayoutOut, IndexLayout& processLayoutOut,
		IndexLayout& deltaNbrLayoutOut, IndexLayout& crossPointLayoutOut,
		const IndexLayout& standardLayout, vector<int>& flags,
		IDomainDecompositionInfo& ddinfo)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	typedef IndexLayout::Interface	Interface;
	typedef IndexLayout::iterator	InterfaceIter;
	typedef IndexLayout::const_iterator	ConstInterfaceIter;
	typedef Interface::Element		Element;
	typedef Interface::const_iterator		ConstElemIter;
	typedef Interface::iterator		ElemIter;

//	the local process and subdomain id
	int localProcID = pcl::GetProcRank();
	int localSubdomID = ddinfo.map_proc_id_to_subdomain_id(localProcID);

	for(ConstInterfaceIter iiter = standardLayout.begin();
		iiter != standardLayout.end(); ++iiter)
	{
	//	connected proc and subdomain ids
		int connProc = standardLayout.proc_id(iiter);
		int connSubdomID = ddinfo.map_proc_id_to_subdomain_id(connProc);
		const Interface& stdInterface = standardLayout.interface(iiter);

	//	in order to avoid creation of unrequired interfaces, we first count
	//	how many entries of each flag-type (-2, -1, >= 0) exist in the current
	//	standard interface.
		int numMult = 0;		// flag: -2
		int numOne = 0;			// flag: >= 0
		int numNone = 0;		// flag: -1
		for(ConstElemIter eiter = stdInterface.begin();
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
				for(ConstElemIter eiter = stdInterface.begin();
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
				for(ConstElemIter eiter = stdInterface.begin();
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
			for(ConstElemIter eiter = stdInterface.begin();
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
			for(ConstElemIter eiter = stdInterface.begin();
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
		const IndexLayout& standardMasters, const IndexLayout& standardSlaves,
		int highestReferencedIndex, IDomainDecompositionInfo& ddinfo)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	some typedefs
	typedef IndexLayout::Interface	Interface;
	typedef IndexLayout::iterator	InterfaceIter;
	typedef Interface::iterator		ElemIter;
	typedef IndexLayout::const_iterator	ConstInterfaceIter;
	typedef Interface::const_iterator	ConstElemIter;

//	the local process and subdomain id
	int localProcID = pcl::GetProcRank();
	int localSubdomID = ddinfo.map_proc_id_to_subdomain_id(localProcID);


//	first we communicate the connections of each entry to all processes
	vector<vector<int> > connectionsToProcs;
	std::vector<std::vector<int> > connectionsToSubDoms;
	CommunicateConnections(connectionsToProcs, connectionsToSubDoms, standardMasters, standardSlaves,
						   highestReferencedIndex, ddinfo);
/*
	for(size_t i = 0; i < connectionsToProcs.size(); ++i){
		UG_LOG("con_" << i << ": ");
		for(size_t j = 0; j < connectionsToProcs[i].size(); ++j){
			UG_LOG(connectionsToProcs[i][j] << " ");
		}
		UG_LOG(endl);
	}
*/
////////////////////////
//	this vector will contain an integer
//	* -1 if the element is not connected to another subdomain.
//	* >= 0 for each element which has a copy in exactly one other subdomain,
//	* -2 for each element which has a copy in more than one other subdomain,
	vector<int> flags(connectionsToProcs.size(), -1);

//	iterate over all standard master interfaces and fill the vector as
//	described above for master elements
	for(ConstInterfaceIter iiter = standardMasters.begin();
		iiter != standardMasters.end(); ++iiter)
	{
	//	connected proc and subdomain ids
		int connProc = standardMasters.proc_id(iiter);
		int connSubdomID = ddinfo.map_proc_id_to_subdomain_id(connProc);
		if(localSubdomID != connSubdomID){
		//	the interface connects two subdomains.
			const Interface& interface = standardMasters.interface(iiter);
			for(ConstElemIter eiter = interface.begin();
				eiter != interface.end(); ++eiter)
			{
				Interface::Element ind = interface.get_element(eiter);
				if(flags[ind] == -1) // not connection was yet known.
					flags[ind] = connSubdomID;
				else if(flags[ind] != connSubdomID)	//	at least one other connection is known.
					flags[ind] = -2;//	and now at least two connections are known.
			}
		}
	}

//todo: instead of communicating the master flags, one could use the
//		connections array to check whether a slave lies on a boundary
//		between different subdomains.
	vector<int> masterFlags(connectionsToProcs.size(), -1);
	ComPol_VecCopy<vector<int> > compolCopy(&masterFlags, &flags);
	InterfaceCommunicator<IndexLayout> com;
	com.send_data(standardMasters, compolCopy);
	com.receive_data(standardSlaves, compolCopy);
	com.communicate();

//	now fill the flags vector for slave entries.
//	here we'll use the connectionsToProcs, since we otherwise wouldn't know all connections.
	for(ConstInterfaceIter iiter = standardSlaves.begin();
		iiter != standardSlaves.end(); ++iiter)
	{
		const Interface& interface = standardSlaves.interface(iiter);
		for(ConstElemIter eiter = interface.begin();
			eiter != interface.end(); ++eiter)
		{
			Interface::Element elem = interface.get_element(eiter);
			if(masterFlags[elem] != -1)
			{
			//	the master connects multiple subdomains.
				if(flags[elem] == -1){
				//	check whether the element is connected to one or to more other subdomains
					vector<int>& connToProc = connectionsToProcs[elem];

					flags[elem] = masterFlags[elem];

					if((flags[elem] >= 0)){
					//	check whether the associated master is in another subdomain
					//	and whether another slave in the same subdomain exists.
						if(ddinfo.map_proc_id_to_subdomain_id(connToProc[0]) != localSubdomID){
							bool encounteredSelf = false;
							for(size_t i = 1; i < connToProc.size(); ++i){
								if(connToProc[i] == localProcID)
									encounteredSelf = true;
								else if(ddinfo.map_proc_id_to_subdomain_id(connToProc[i])
									== localSubdomID)
								{
								//	at this point we will take care of interface insertion right here.
									//flags[elem] = -3;

									if(encounteredSelf){
									//	the elem will be master (note that this if can be encountered
									//	multiple times).
										Interface& interface = processMastersOut.interface(connToProc[i]);
										interface.push_back(elem);
									}
									else{
									//	the elem will be slave to connToProc[i], since this is the first
									//	occurence of a process in this subdomain.
										Interface& interface = processSlavesOut.interface(connToProc[i]);
										interface.push_back(elem);
										break;
									}
								}
							}
						}
					}

				}
			}
		}
	}

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

/* This is an old implementation which does not create comlplete process layouts
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

	pcl::InterfaceCommunicator<IndexLayout> interfaceComm;

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
}
*/



/**
 * adds connections between slave nodes to the interfaces/layouts
 * when a master node has 2 slave nodes, this function adds a connection between these nodes
 *
 * \param communicator used InterfaceCommunicator
 * \param masterLayout master layout
 * \param slaveLayout slave layout
 * \param allToAllSend layout with slave-slave connections at the end of this function
 * \param allToAllReceive layout with slave-slave connections at the end of this function
 *
 * since we have slave->slave and slave<-slave connections, indices added to layouts
 * are always added to allToAllSend AND allToAllReceive.
 *
 * \note this function ONLY adds slave-slave connections. if you need master->slave connections, try
 *  AddLayout(allToAllSend, masterLayout);
 *  AddLayout(allToAllReceive, slaveLayout);
 * for slave->master
 *  AddLayout(allToAllSend, slaveLayout);
 *  AddLayout(allToAllReceive, masterLayout);
 *
 * \note because the order in the interfaces is important, this function is more complicate that one would expect.
 */
void AddConnectionsBetweenSlaves(pcl::InterfaceCommunicator<IndexLayout> &communicator,
		IndexLayout &masterLayout, IndexLayout &slaveLayout, IndexLayout &allToAllSend,
		IndexLayout &allToAllReceive)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	// slave -> slave
#if 1
	// 1. get for every master node where it is slave
	std::map<size_t, std::vector<int> > slaveOnProc;

	// localToLocal[pid][i] is the position of i in the interface to processor pid.
	std::map<size_t, std::map<size_t, size_t> > localToInterfaceIndex;
	//slaveOnProc.resize(A_OL2.num_rows());

	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		IndexLayout::Interface &interface = masterLayout.interface(iter);
		int pid = masterLayout.proc_id(iter);

		std::map<size_t, size_t> &localToInterfaceIndexMap = localToInterfaceIndex[pid];
		size_t i=0;
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2, ++i)
		{
			size_t index = interface.get_element(iter2);
			slaveOnProc[index].push_back(pid);
			localToInterfaceIndexMap[index] = i;
		}
	}

	typedef std::map<int, BinaryBuffer> BufferMap;
	BufferMap sendpack;
	//UG_LOG("\n\n");
	for(std::map<size_t, std::vector<int> >::iterator it = slaveOnProc.begin(); it != slaveOnProc.end(); ++it)
	{
		std::vector<int> &procs = it->second;
		if(procs.size() <= 1) continue;
		size_t index = it->first;
		//UG_LOG("index " << index << " is to processors ");
		for(size_t i=0; i<procs.size(); i++)
		{
			BinaryBuffer &stream = sendpack[procs[i]];
			size_t interfaceIndex = localToInterfaceIndex[procs[i]][index];
			//UG_LOG(procs[i] << " (interfaceIndex " << interfaceIndex << ") ");
			for(size_t j=0; j<procs.size(); j++)
			{
				if(i == j) continue;
				UG_ASSERT(procs[i] != procs[j], procs[i] << " != " << procs[j]);
				Serialize(stream, interfaceIndex);
				Serialize(stream, procs[j]);
			}
		}
		//UG_LOG("\n");
	}

	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		int pid = masterLayout.proc_id(iter);
		BinaryBuffer &stream = sendpack[pid];
		communicator.send_raw(pid, stream.buffer(), stream.write_pos(), false);
		//UG_LOG("Sending " << stream.size() << " bytes of data to processor " << pid << ":\n");
	}

	// 3. communicate
	BufferMap receivepack;
	std::vector<int> pids;

	for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
	{
		int pid = slaveLayout.proc_id(iter);
		pids.push_back(pid);
		communicator.receive_raw(pid, receivepack[pid]);
	}
	communicator.communicate();

	/* sort pids. this is important!
	 * example: processor A has two slave nodes in common with processor B, say GID 0 and GID 1
	 * so if the interface A -> B looks like 0, 1, the interface B->A must be 0, 1.
	 * to assure this, indices with lower master PID are added first AND we use the order which we got from
	 * (otherwise we could get A->B: 0,1 and B->A: 1,0)
	 */
	sort(pids.begin(), pids.end());

	// 4. process data

	for(size_t i=0; i<pids.size(); i++)
	{
		int pid = pids[i];
		BinaryBuffer &stream = receivepack[pid];
		//UG_LOG("Received " << stream.size() << " bytes of data from processor " << pid << ":\n");

		std::vector<size_t> indices;
		IndexLayout::Interface &interface = slaveLayout.interface(pid);
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			indices.push_back(interface.get_element(iter2));

		while(!stream.eof())
		{
			size_t index;
			Deserialize(stream, index);
			index = indices[index];
			int pid2;
			Deserialize(stream, pid2);

			//UG_LOG(" got " << s << " other slave connections from node " << index << ": ")
			UG_ASSERT(pid2 != pcl::GetProcRank(), "");
			allToAllReceive.interface(pid2).push_back(index);
			allToAllSend.interface(pid2).push_back(index);
			//UG_LOG("Added Index " << index << " to interface with processor " << pid2 << "\n");
		}

	}

	//UG_LOG("\n\n");
#else
	// this implementation, using CommunicateConnections, unfortunately does not work at the moment (sorting issue below)
	std::vector<std::vector<int> > connectionsToProcs;
	std::vector<std::vector<int> > connectionsToSubDoms;
	CommunicateConnections(connectionsToProcs, connectionsToSubDoms, masterLayout, slaveLayout, A_OL2.num_rows(), ddinfo);

	std::vector<int> pids;

	for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
	{
		int pid = slaveLayout.proc_id(iter);
		if(pid != pcl::GetProcRank())
			pids.push_back(pid);
	}

	/* sort pids. this is important!
	 * example: processor A has two slave nodes in common with processor B, say GID 0 and GID 1
	 * so if the interface A -> B looks like 0, 1, the interface B->A must be 0, 1.
	 * to assure this, indices with lower master PID are added first
	 * (otherwise we could get A->B: 0,1 and B->A: 1,0) */
	sort(pids.begin(), pids.end());

		// add first all connectionsToProcs from the lowest master PID
	for(std::vector<int>::iterator masterPIDit = pids.begin(); masterPIDit != pids.end(); ++masterPIDit)
	{
		int masterPID = (*masterPIDit);
		//UG_LOG("Processing pid " << masterPID << ":\n");

		// this does not work, since the indices on this node can have a different ordering as on the other processing node
		// but we need to insert them into the layout like on the other processing node
		// if you can fix this, you can use this function again
		for(size_t i=0; i<connectionsToProcs.size(); i++)
		{
			std::vector<int> &connToProc = connectionsToProcs[i];
			size_t index = i;
			//	the first entry in each connection is the master elements process
			if(connToProc.size() <= 1 || connToProc[0] != masterPID) continue;
			//UG_LOG("index " << index << "added to ");
			for(size_t j=1; j<connToProc.size(); j++)
			{
				int pid = connToProc[j];
				if(pid == pcl::GetProcRank()) continue;
				//UG_LOG(pid << " ");
				OLCoarseningReceiveLayout.interface(pid).push_back(index);
				OLCoarseningSendLayout.interface(pid).push_back(index);
			}
			//UG_LOG("\n");
		}
	}
#endif

}

void CreateAllToAllFromMasterSlave(pcl::InterfaceCommunicator<IndexLayout> &communicator,
		IndexLayout &OLCoarseningSendLayout, IndexLayout &OLCoarseningReceiveLayout,
		IndexLayout &OL1MasterLayout, IndexLayout &OL1SlaveLayout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	// master -> slave
	AddLayout(OLCoarseningSendLayout, OL1MasterLayout);
	AddLayout(OLCoarseningReceiveLayout, OL1SlaveLayout);

	// slave -> master
	AddLayout(OLCoarseningSendLayout, OL1SlaveLayout);
	AddLayout(OLCoarseningReceiveLayout, OL1MasterLayout);


	AddConnectionsBetweenSlaves(communicator,
			OL1MasterLayout, OL1SlaveLayout,
			OLCoarseningSendLayout, OLCoarseningReceiveLayout);
}

SmartPtr<AlgebraLayouts> CreateLocalAlgebraLayouts()
{
	AlgebraLayouts *p = new AlgebraLayouts;
	p->clear();
	p->proc_comm() = pcl::ProcessCommunicator(pcl::PCD_LOCAL);
	return SmartPtr<AlgebraLayouts>(p);
}


}// end of namespace
