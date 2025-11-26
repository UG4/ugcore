/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <string>
#include <vector>

#include "common/log.h"
#include "common/util/file_util.h"
#include "common/util/binary_buffer.h"
#include "common/serialization.h"

#include "pcl_util.h"
#include "pcl_profiling.h"

using namespace std;
using namespace ug;

namespace pcl{

////////////////////////////////////////////////////////////////////////////////
void SynchronizeProcesses()
{
	ProcessCommunicator().barrier();
}

////////////////////////////////////////////////////////////////////////////////
bool AllProcsTrue(bool bFlag, const ProcessCommunicator& comm)
{
	if(comm.is_local() || comm.empty()) return bFlag;
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
bool OneProcTrue(bool bFlag, const ProcessCommunicator &comm)
{
	if(comm.is_local() || comm.empty()) return bFlag;
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
	PCL_PROFILE(pcl_CommunicateInvolvedProcesses);

	using namespace std;

	vReceiveFromRanksOut.clear();

	if(procComm.empty())
		return;

	const int localProcRank = ProcRank();

//	we'll use an array in which we'll store the number of
//	processes, to which each process wants to talk.
	vector<int> vNumAssProcs(procComm.size());

//	perform allgather with proc-numbers
	int procCount = static_cast<int>(vSendToRanks.size());

	procComm.allgather(&procCount, 1, PCL_DT_INT,
					   vNumAssProcs.data(),
					   1,
					   PCL_DT_INT);

//	sum the number of processes so that the absolute list size can
//	be determined.
	size_t listSize = 0;
	vector<int> vDisplacements(vNumAssProcs.size());
	for(size_t i = 0; i < vNumAssProcs.size(); ++i){
		vDisplacements[i] = static_cast<int>(listSize);
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
	const int ssize = static_cast<int>(vNumAssProcs.size());
	for(int i = 0; i <ssize; ++i)
	{
	//	check whether the i-th proc wants to communicate with the local proc
		for(int j = 0; j < vNumAssProcs[i]; ++j)
		{
			if(vGlobalProcList[vDisplacements[i] + j] == localProcRank)
			{
				vReceiveFromRanksOut.push_back(procComm.get_proc_id(i));
			//	the j-th proc is handled completely. resume with the next.
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
	PCL_PROFILE_FUNC();
//	we overwrite some data in recvFrom - that's why we need a copy
	std::vector<int> recvFrom = recvFromTmp;
	
//	make sure that all processes know, who is sending data to them
	vector<int> sendingProcs;
	CommunicateInvolvedProcesses(sendingProcs, sendTo, involvedProcs);
	
//	check whether the list of sendingProcs and the list of recvFrom procs matches.
//	we do this by setting each entry in recvFrom, which also lies in sendingProcs to -1.
	bool sendRecvMismatch = false;
	const size_t ssize = sendingProcs.size();
	for(size_t i = 0; i < ssize; ++i){
		int rank = sendingProcs[i];
		auto findIter = find(recvFrom.begin(), recvFrom.end(), rank);
		if(findIter != recvFrom.end())
			*findIter = -1;
		else{
			UG_LOG("ERROR: send / receive mismatch: ");
			UG_LOG("proc " << rank << " sends to proc " << ProcRank());
			UG_LOG(" but no matching receive is scheduled.\n");
			sendRecvMismatch = true;
		}
	}
	
//	now check whether an entry != -1 still resides in recvFrom.
	const size_t ssize_recv = recvFrom.size();
	for(size_t i = 0; i < ssize_recv; ++i){
		if(recvFrom[i] != -1){
			UG_LOG("ERROR: receive / send mismatch: ");
			UG_LOG("proc " << ProcRank() << " awaits data from proc " << recvFrom[i]);
			UG_LOG(", but no send was scheduled.\n");
			sendRecvMismatch = true;
		}
	}
	
//	if a mismatch occurred on one, then we'll gather all procs which failed.
	if(!AllProcsTrue(!sendRecvMismatch, involvedProcs)){
		int mismatch = sendRecvMismatch;
		vector buffer(involvedProcs.size(), 0);
		int root = GetLogAssistant().get_output_process();
		if(root < 0) root = 0;
		involvedProcs.gather(&mismatch, 1, PCL_DT_INT, &buffer.front(), 1,
							PCL_DT_INT, root);

		UG_LOG("SEND / RECEIVE MISMATCH OCCURRED ON PROC:");
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
	PCL_PROFILE_FUNC();
	assert(recvFrom.size() == recvBufSizes.size());
	assert(sendTo.size() == sendBufSizes.size());
	
//	number of in and out-streams.
	const size_t	numOutStreams = sendTo.size();
	const size_t	numInStreams = recvFrom.size();
	
//	used for mpi-communication.
	std::vector<MPI_Request> vSendRequests(numOutStreams);
	std::vector<MPI_Request> vReceiveRequests(numInStreams);

	constexpr int testTag = 744444;//	an arbitrary number

	std::vector<int> vSendBufSizes(numInStreams);
	for(size_t i = 0; i < recvFrom.size(); ++i)
	{
		MPI_Irecv(&vSendBufSizes[i], 1, MPI_INT, recvFrom[i], testTag,
				  PCL_COMM_WORLD, &vReceiveRequests[i]);
	}

	for(size_t i = 0; i < sendTo.size(); ++i)
	{
		int s = sendBufSizes[i];
		MPI_Isend(&s, 1, MPI_INT, sendTo[i], testTag, PCL_COMM_WORLD,
				  &vSendRequests[i]);
	}

	Waitall(vReceiveRequests, vSendRequests);

	bool bSuccess = true;
	for(size_t i = 0; i < recvFrom.size(); ++i){
		if(recvBufSizes[i] != vSendBufSizes[i])
		{
			UG_LOG("SEND / RECEIVE BUFFER MISMATCH: "
					<< "receive buffer on proc " << ProcRank()
					<< " has "<< recvBufSizes[i]
					<< " bytes, but send buffer on proc " << recvFrom[i]
					<< " has " << vSendBufSizes[i] << " bytes\n");
			bSuccess = false;
		}
	}

//	if a mismatch occurred on one, then we'll gather all procs which failed.
	if(!AllProcsTrue(bSuccess, involvedProcs)){
		int mismatch = !bSuccess;
		vector buffer(involvedProcs.size(), 0);
		int root = GetLogAssistant().get_output_process();
		if(root < 0) root = 0;
		involvedProcs.gather(&mismatch, 1, PCL_DT_INT, &buffer.front(), 1,
							PCL_DT_INT, root);

		UG_LOG("SEND / RECEIVE BUFFER MISMATCH OCCURRED ON PROC:");
		const size_t ssize_buffer = buffer.size();
		for(size_t i = 0; i < ssize_buffer; ++i){
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
	PCL_PROFILE_FUNC();
	if(bDistributedLoad == false)
		return ReadFile(filename.c_str(), file, bText);

	BinaryBuffer buf;
	bool bSuccess=false;
	if(ProcRank() == pc.get_proc_id(0))
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

	if(ProcRank() != pc.get_proc_id(0))
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
