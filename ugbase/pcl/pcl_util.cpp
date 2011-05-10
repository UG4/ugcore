// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 17.03.2011 (m,d,y)
 
#include <vector>
#include "pcl_util.h"
#include "pcl_profiling.h"
#include "common/log.h"

using namespace std;

namespace pcl{

void CommunicateInvolvedProcesses(std::vector<int>& vReceiveFromRanksOut,
								  std::vector<int>& vSendToRanks,
								  const ProcessCommunicator& procComm)
{
	PCL_PROFILE(pcl_CommunicateInvolvedProcesses);

	using namespace std;

	vReceiveFromRanksOut.clear();

//	if there is only one process in the communicator, there's
//	nothing to do.
	if(procComm.size() < 2)
		return;

	const int localProcRank = GetProcRank();

//	we'll use an array in which we'll store the number of
//	processes, to which each process wants to talk.
	vector<int> vNumAssProcs(procComm.size());

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
	//	check whether the i-th proc wants to communicate with the local proc
		for(int j = 0; j < vNumAssProcs[i]; ++j)
		{
			if(vGlobalProcList[vDisplacements[i] + j] == localProcRank)
			{
				vReceiveFromRanksOut.push_back(procComm.get_proc_id(i));
			//	the j-th proc is handled completly. resume with the next.
				break;
			}
		}
	}

//	vProcRanksInOut should now contain all process-ranks with which
//	the local proc should communicate.
}



bool StreamPacksMatch(ug::StreamPack& streamPackRecv, ug::StreamPack& streamPackSend,
					  const ProcessCommunicator& involvedProcs)
{
//	create a vector containing all recv-from ranks and a vector containing all send-to ranks.
	vector<int> recvFrom, sendTo;
	
	for(ug::StreamPack::iterator iter = streamPackRecv.begin();
		iter != streamPackRecv.end(); ++iter)
	{
		recvFrom.push_back(iter->first);
	}
	
	for(ug::StreamPack::iterator iter = streamPackSend.begin();
		iter != streamPackSend.end(); ++iter)
	{
		sendTo.push_back(iter->first);
	}
	
//	make sure that all processes know, who is sending data to them
	vector<int> sendingProcs;
	CommunicateInvolvedProcesses(sendingProcs, sendTo, involvedProcs);
	
//	check whether the list of sendingProcs and the list of recvFrom procs matches.
//	we do this by setting each entry in recvFrom, which also lies in sendingProcs to -1.
	bool sendRecvMismatch = false;
	for(size_t i = 0; i < sendingProcs.size(); ++i){
		int rank = sendingProcs[i];
		vector<int>::iterator findIter = find(recvFrom.begin(), recvFrom.end(), rank);
		if(findIter != recvFrom.end())
			*findIter = -1;
		else{
			UG_LOG("ERROR: send / receive mismatch: ");
			UG_LOG("proc " << rank << " sends to proc " << GetProcRank());
			UG_LOG(" but no matching receive is scheduled.\n");
			sendRecvMismatch = true;
		}
	}
	
//	now check whether an entry != -1 still resides in recvFrom.
	for(size_t i = 0; i < recvFrom.size(); ++i){
		if(recvFrom[i] != -1){
			UG_LOG("ERROR: receive / send mismatch: ");
			UG_LOG("proc " << GetProcRank() << " awaits data from proc " << recvFrom[i]);
			UG_LOG(", but no send was scheduled.\n");
			sendRecvMismatch = true;
		}
	}
	
//	if a mismatch occurred on one, then we'll gather all procs which failed.
	if(!AllProcsTrue(!sendRecvMismatch, involvedProcs)){
		int mismatch = (int)sendRecvMismatch;
		vector<int> buffer(involvedProcs.size(), 0);
		involvedProcs.gather(&mismatch, 1, PCL_DT_INT, &buffer.front(), 1,
							PCL_DT_INT, GetOutputProcRank());

		UG_LOG("SEND / RECEIVE MISMATCH OCCURED ON PROC:");
		for(size_t i = 0; i < buffer.size(); ++i){
			if(buffer[i] != 0){
			//	a mismatch occurred
				UG_LOG(" " << involvedProcs.get_proc_id(i));
			}
		}
		UG_LOG(endl);
		return false;
	}
	
	return true;
}


bool StreamPackBuffersMatch(ug::StreamPack &streamPackRecv,
							ug::StreamPack &streamPackSend,
							const ProcessCommunicator& involvedProcs)
{
//	number of in and out-streams.
	size_t	numOutStreams = streamPackSend.num_streams();
	size_t	numInStreams = streamPackRecv.num_streams();
	
//	used for mpi-communication.
	std::vector<MPI_Request> vSendRequests(numOutStreams);
	std::vector<MPI_Request> vReceiveRequests(numInStreams);

//	wait until data has been received
	std::vector<MPI_Status> vReceiveStates(numInStreams);
	std::vector<MPI_Status> vSendStates(numOutStreams);

	int testTag = 744444;//	an arbitrary number
	int counter = 0;

	std::vector<int> vSendBufSizes(numInStreams);
	for(ug::StreamPack::iterator iter = streamPackRecv.begin();
		iter != streamPackRecv.end(); ++iter, ++counter)
	{
		MPI_Irecv(&vSendBufSizes[counter], 1, MPI_INT, iter->first, testTag,
				  MPI_COMM_WORLD, &vReceiveRequests[counter]);
	}
	
	counter = 0;
	for(ug::StreamPack::iterator iter = streamPackSend.begin();
		iter != streamPackSend.end(); ++iter, ++counter)
	{
		int s = iter->second->size();
		MPI_Isend(&s, 1, MPI_INT, iter->first, testTag, MPI_COMM_WORLD,
				  &vSendRequests[counter]);
	}

	MPI_Waitall(numInStreams, &vReceiveRequests[0], &vReceiveStates[0]);
	MPI_Waitall(numOutStreams, &vSendRequests[0], &vSendStates[0]);

	bool bSuccess = true;
	counter = 0;
	for(ug::StreamPack::iterator iter = streamPackRecv.begin();
		iter != streamPackRecv.end(); ++iter, ++counter)
	{
		if((int)iter->second->size() != vSendBufSizes[counter])
		{
			UG_LOG("SEND / RECEIVE BUFFER MISMATCH: "
					<< "receive buffer on proc " << GetProcRank()
					<< " has "<< iter->second->size()
					<< " bytes, but send buffer on proc " << iter->first
					<< " has " << vSendBufSizes[counter] << " bytes\n");
			bSuccess = false;
		}
	}

//	if a mismatch occurred on one, then we'll gather all procs which failed.
	if(!AllProcsTrue(bSuccess, involvedProcs)){
		int mismatch = (int)(!bSuccess);
		vector<int> buffer(involvedProcs.size(), 0);
		involvedProcs.gather(&mismatch, 1, PCL_DT_INT, &buffer.front(), 1,
							PCL_DT_INT, GetOutputProcRank());

		UG_LOG("SEND / RECEIVE BUFFER MISMATCH OCCURED ON PROC:");
		for(size_t i = 0; i < buffer.size(); ++i){
			if(buffer[i] != 0){
			//	a mismatch occurred
				UG_LOG(" " << involvedProcs.get_proc_id(i));
			}
		}
		UG_LOG(endl);
		return false;
	}
	
	return true;
}

}// end of namespace
