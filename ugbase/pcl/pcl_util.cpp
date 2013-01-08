// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 17.03.2011 (m,d,y)
 
#include <vector>
#include "pcl_util.h"
#include "pcl_profiling.h"
#include "common/log.h"
#include <string>
#include "common/util/file_util.h"
#include "common/util/binary_buffer.h"
#include "common/serialization.h"

using namespace std;
using namespace ug;

namespace pcl{

////////////////////////////////////////////////////////////////////////////////
void SynchronizeProcesses()
{
	ProcessCommunicator().barrier();
}

////////////////////////////////////////////////////////////////////////////////
bool AllProcsTrue(bool bFlag, ProcessCommunicator comm)
{
	PCL_DEBUG_BARRIER(comm);
	PCL_PROFILE(pclAllProcsTrue);

//	local int bool flag
	int boolFlag = (bFlag) ? 1 : 0;

// 	local return flag
	int retBoolFlag;

//	all reduce
	comm.allreduce(&boolFlag, &retBoolFlag, 1, PCL_DT_INT, PCL_RO_LAND);

//	return global flag
	return retBoolFlag != 0;
}

////////////////////////////////////////////////////////////////////////////////
bool OneProcTrue(bool bFlag, ProcessCommunicator comm)
{
	PCL_DEBUG_BARRIER(comm);
	PCL_PROFILE(pclAllProcsTrue);

//	local int bool flag
	int boolFlag = (bFlag) ? 1 : 0;

// 	local return flag
	int retBoolFlag;

//	all reduce
	comm.allreduce(&boolFlag, &retBoolFlag, 1, PCL_DT_INT, PCL_RO_LOR);

//	return global flag
	return retBoolFlag != 0;
}

////////////////////////////////////////////////////////////////////////////////
void CommunicateInvolvedProcesses(std::vector<int>& vReceiveFromRanksOut,
								  const std::vector<int>& vSendToRanks,
								  const ProcessCommunicator& procComm)
{
	PCL_DEBUG_BARRIER(procComm);
	PCL_PROFILE(pcl_CommunicateInvolvedProcesses);

	using namespace std;

	vReceiveFromRanksOut.clear();

	if(!procComm.size())
		return;

	const int localProcRank = GetProcRank();

//	we'll use an array in which we'll store the number of
//	processes, to which each process wants to talk.
	vector<int> vNumAssProcs(procComm.size());

//	perform allgather with proc-numbers
	int procCount = (int)vSendToRanks.size();

	procComm.allgather(&procCount, 1, PCL_DT_INT,
					   GetDataPtr(vNumAssProcs),
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

	if(!listSize)
		return;

//	perform allgather with proc-lists
//	this list will later on contain an adjacency-list for each proc.
//	adjacency in this case means, that processes want to communicate.
	vector<int> vGlobalProcList(listSize);

	procComm.allgatherv(GetDataPtr(vSendToRanks),
						procCount,
						PCL_DT_INT,
						GetDataPtr(vGlobalProcList),
						GetDataPtr(vNumAssProcs),
						GetDataPtr(vDisplacements),
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

bool SendRecvListsMatch(const std::vector<int>& recvFromTmp,
						const std::vector<int>& sendTo,
						const ProcessCommunicator& involvedProcs)
{
//	we overwrite some data in recvFrom - thats why we need a copy
	std::vector<int> recvFrom = recvFromTmp;
	
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
		int root = GetLogAssistant().get_output_process();
		if(root < 0) root = 0;
		involvedProcs.gather(&mismatch, 1, PCL_DT_INT, &buffer.front(), 1,
							PCL_DT_INT, root);

		UG_LOG("SEND / RECEIVE MISMATCH OCCURED ON PROC:");
		for(size_t i = 0; i < buffer.size(); ++i){
			if(buffer[i] != 0){
			//	a mismatch occurred
				UG_LOG(" " << involvedProcs.get_proc_id(i));
			}
		}
		UG_LOG(endl);
		UG_LOG(endl);
		return false;
	}
	
	return true;
}


bool SendRecvBuffersMatch(const std::vector<int>& recvFrom, const std::vector<int>& recvBufSizes,
						  const std::vector<int>& sendTo, const std::vector<int>& sendBufSizes,
						  const ProcessCommunicator& involvedProcs)
{
	assert(recvFrom.size() == recvBufSizes.size());
	assert(sendTo.size() == sendBufSizes.size());
	
//	number of in and out-streams.
	size_t	numOutStreams = sendTo.size();
	size_t	numInStreams = recvFrom.size();
	
//	used for mpi-communication.
	std::vector<MPI_Request> vSendRequests(numOutStreams);
	std::vector<MPI_Request> vReceiveRequests(numInStreams);

	int testTag = 744444;//	an arbitrary number

	std::vector<int> vSendBufSizes(numInStreams);
	for(size_t i = 0; i < recvFrom.size(); ++i)
	{
		MPI_Irecv(&vSendBufSizes[i], 1, MPI_INT, recvFrom[i], testTag,
				  MPI_COMM_WORLD, &vReceiveRequests[i]);
	}

	for(size_t i = 0; i < sendTo.size(); ++i)
	{
		int s = sendBufSizes[i];
		MPI_Isend(&s, 1, MPI_INT, sendTo[i], testTag, MPI_COMM_WORLD,
				  &vSendRequests[i]);
	}

	Waitall(vReceiveRequests, vSendRequests);

	bool bSuccess = true;
	for(size_t i = 0; i < recvFrom.size(); ++i){
		if(recvBufSizes[i] != vSendBufSizes[i])
		{
			UG_LOG("SEND / RECEIVE BUFFER MISMATCH: "
					<< "receive buffer on proc " << GetProcRank()
					<< " has "<< recvBufSizes[i]
					<< " bytes, but send buffer on proc " << recvFrom[i]
					<< " has " << vSendBufSizes[i] << " bytes\n");
			bSuccess = false;
		}
	}

//	if a mismatch occurred on one, then we'll gather all procs which failed.
	if(!AllProcsTrue(bSuccess, involvedProcs)){
		int mismatch = (int)(!bSuccess);
		vector<int> buffer(involvedProcs.size(), 0);
		int root = GetLogAssistant().get_output_process();
		if(root < 0) root = 0;
		involvedProcs.gather(&mismatch, 1, PCL_DT_INT, &buffer.front(), 1,
							PCL_DT_INT, root);

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


bool ParallelReadFile(string &filename, vector<char> &file, bool bText, bool bDistributedLoad, const ProcessCommunicator& pc)
{
	if(bDistributedLoad == false)
		return ReadFile(filename.c_str(), file, bText);

	BinaryBuffer buf;
	bool bSuccess;
	if(GetProcRank() == pc.get_proc_id(0))
	{
		bSuccess = ReadFile(filename.c_str(), file, bText);
		Serialize(buf, bSuccess);
		if(bSuccess)
		{
			Serialize(buf, filename);
			Serialize(buf, file);
		}
	}

	pc.broadcast(buf);

	if(GetProcRank() != pc.get_proc_id(0))
	{
		Deserialize(buf, bSuccess);
		if(bSuccess)
		{
			Deserialize(buf, filename);
			Deserialize(buf, file);
		}
	}
	return bSuccess;
}

}// end of namespace
