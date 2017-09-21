/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include <vector>
#include "pcl_methods.h"
#include "common/util/smart_pointer.h"
#include "pcl_process_communicator.h"
#include "common/log.h"
#include "common/assert.h"
#include "common/util/vector_util.h"
#include "pcl_profiling.h"
#include "pcl_datatype.h"
#include "pcl_util.h"

using namespace std;
using namespace ug;

namespace pcl
{


ProcessCommunicator::
ProcessCommunicator(ProcessCommunicatorDefaults pcd)
{
	switch(pcd)
	{
		case PCD_EMPTY:
			m_comm = SPCommWrapper(new CommWrapper(MPI_COMM_NULL, false));
			break;

		case PCD_WORLD:
			m_comm = SPCommWrapper(new CommWrapper(PCL_COMM_WORLD, false));
			break;

		case PCD_LOCAL:
			// Do nothing since m_comm.valid() == false indicates that the ProcCom operates locally only.
			break;
	}
}

size_t ProcessCommunicator::
size() const
{
	if(is_local()) return 1;
	if(m_comm->m_mpiComm == MPI_COMM_NULL)
		return 0;
	
	int size;
	if(MPI_Comm_size(m_comm->m_mpiComm, &size) == MPI_SUCCESS)
		return (size_t)size;
		
	UG_LOG("  ERROR in ProcessCommunicator::size(): Unknown MPI Error. Returning 0.\n");
	return 0;
}

int ProcessCommunicator::
get_proc_id(size_t groupIndex) const
{
	if(is_local()) return pcl::ProcRank();
	if(m_comm->m_mpiComm == PCL_COMM_WORLD)
		return (int)groupIndex;
	return m_comm->m_procs[groupIndex];
}

int ProcessCommunicator::
get_local_proc_id(int globalProcID) const
{
	if(is_local()) return 0;
	if(m_comm->m_mpiComm == PCL_COMM_WORLD)
		return globalProcID;

	const vector<int>& procs = m_comm->m_procs;
	if(globalProcID == pcl::ProcRank())
	{
		int rank;
		MPI_Comm_rank(m_comm->m_mpiComm, &rank);
		UG_ASSERT(procs[rank] == pcl::ProcRank(), "?");
		return rank;
	}

	for(size_t i = 0; i < procs.size(); ++i){
		if(globalProcID == procs[i])
			return (int)i;
	}

	UG_LOG("  ERROR in ProcessCommunicator::get_local_proc_id(): given id not contained in local proc list.\n");
	return -1;
}


ProcessCommunicator
ProcessCommunicator::
create_sub_communicator(bool participate) const
{
	UG_COND_THROW(is_local(), "not available");
	PCL_PROFILE(pcl_ProcCom_create_sub_com__participate);

//	if the current communicator is empty theres nothing to do
	if(empty())
		return ProcessCommunicator(PCD_EMPTY);

//	get the number of processes in the current communicator
	int size;
	MPI_Comm_size(m_comm->m_mpiComm, &size);

	if(size < 1)
		return ProcessCommunicator(PCD_EMPTY);

//	create a buffer and initialise it with 0
	vector<int> srcArray(size, 0);

//	if the process wants to participate, set his entry to 1
	int rank;
	MPI_Comm_rank(m_comm->m_mpiComm, &rank);
	if(participate)
		srcArray[rank] = 1;

//	synchronize the newProcs array between all processes in the communicator
	vector<int> destArray(size, 0);

	{
		PCL_PROFILE(waiting_before_allreduce);
		PCL_DEBUG_BARRIER((*this));
	}

	allreduce(&srcArray.front(), &destArray.front(),
			   size, PCL_DT_INT, PCL_RO_MAX);

	{
		PCL_PROFILE(waiting_after_allreduce);
		PCL_DEBUG_BARRIER((*this));
	}

//	build a local array that holds all the procs that shall go
//	into the new communicator
	vector<int> newProcs;
	newProcs.reserve(size);

	// note: ranks are ranks in the (group!) communicator m_comm->m_mpiComm
	// since we building the new group relative to the old,
	// we add to newProcs the group ranks
	// these are NOT the global ranks like in pcl::ProcRank
	for(size_t i = 0; i < destArray.size(); ++i){
		if(destArray[i])
			newProcs.push_back(i);
	}

//	if newProcs is not empty, we'll build a new mpi-communicator.
	if(newProcs.empty())
		return ProcessCommunicator(PCD_EMPTY);
	else
		return create_sub_communicator(newProcs);
}

/// note: ranks in newProcs are ranks in the (group!) communicator m_comm->m_mpiComm
/// @sa create_communicator
ProcessCommunicator
ProcessCommunicator::
create_sub_communicator(vector<int> &newProcs) const
{
	UG_COND_THROW(is_local(), "not available");
	PCL_PROFILE(pcl_ProcCom_create_sub_com__array);
	if(newProcs.size() == 0)
		return ProcessCommunicator(PCD_EMPTY);
	if((int)newProcs.size() == NumProcs())
		return ProcessCommunicator(PCD_WORLD);

	PCL_PROFILE(create_mpi_com);
	MPI_Group grpOld;
	MPI_Group grpNew;
	MPI_Comm commNew;

	MPI_Comm_group(m_comm->m_mpiComm, &grpOld);
	MPI_Group_incl(grpOld, (int)newProcs.size(), &newProcs.front(), &grpNew);
	MPI_Comm_create(m_comm->m_mpiComm, grpNew, &commNew);
	PCL_PROFILE_END();

//	create a new ProcessCommunicator
//	if the process is not participating, MPI_Comm_create will return MPI_COMM_NULL
	if(commNew == MPI_COMM_NULL)
		return ProcessCommunicator(PCD_EMPTY);
	else if(commNew == PCL_COMM_WORLD)
		return ProcessCommunicator(PCD_WORLD);

	PCL_PROFILE(get_global_ranks);
	// calculate global ranks for our newProcs array:
	for(size_t i = 0; i < newProcs.size(); ++i)
		newProcs[i] = get_proc_id(newProcs[i]);
	PCL_PROFILE_END();

	// note: since get_proc_rank uses newProcs, don't sort newProcs
	// otherwise above code won't work. here it is now
	// newProcs[group rank] = global rank. (!)

//	the process participates - create the ProcessCommunicator
	ProcessCommunicator newProcComm;
	newProcComm.m_comm = SPCommWrapper(new CommWrapper(commNew, true));
	newProcComm.m_comm->m_procs = newProcs;

	return newProcComm;
}

ProcessCommunicator
ProcessCommunicator::
create_communicator(vector<int> &newGlobalProcs)
{
	CommWrapper comm(PCL_COMM_WORLD, false);

	MPI_Group grpWorld;
	MPI_Group grpNew;
	MPI_Comm commNew;

	MPI_Comm_group(PCL_COMM_WORLD, &grpWorld);
	MPI_Group_incl(grpWorld, (int)newGlobalProcs.size(), &newGlobalProcs.front(), &grpNew);
	MPI_Comm_create(PCL_COMM_WORLD, grpNew, &commNew);

//	create a new ProcessCommunicator
//	if the process is not participating, MPI_Comm_create will return MPI_COMM_NULL
	if(commNew == MPI_COMM_NULL)
		return ProcessCommunicator(PCD_EMPTY);
	else if(commNew == PCL_COMM_WORLD)
		return ProcessCommunicator(PCD_WORLD);

	ProcessCommunicator newProcComm;
	newProcComm.m_comm = SPCommWrapper(new CommWrapper(commNew, true));
	newProcComm.m_comm->m_procs = newGlobalProcs;

	return newProcComm;
}

ProcessCommunicator
ProcessCommunicator::
create_communicator(size_t first, size_t num)
{
	MPI_Group grpWorld;
	MPI_Group grpNew;
	MPI_Comm commNew;

	ProcessCommunicator newProcComm;
	newProcComm.m_comm = SPCommWrapper(new CommWrapper());
	vector<int>& procs = newProcComm.m_comm->m_procs;

	procs.resize(num);
	for(size_t i = 0; i < num; ++i)
		procs[i] = first + i;

	MPI_Comm_group(PCL_COMM_WORLD, &grpWorld);
	MPI_Group_incl(grpWorld, (int)procs.size(), &procs.front(), &grpNew);
	MPI_Comm_create(PCL_COMM_WORLD, grpNew, &commNew);

	//	create a new ProcessCommunicator
//	if the process is not participating, MPI_Comm_create will return MPI_COMM_NULL
	if(commNew == MPI_COMM_NULL)
		return ProcessCommunicator(PCD_EMPTY);
	else if(commNew == PCL_COMM_WORLD)
		return ProcessCommunicator(PCD_WORLD);

	newProcComm.m_comm->m_mpiComm = commNew;
	newProcComm.m_comm->m_bReleaseCommunicator = true;

	return newProcComm;
}


void
ProcessCommunicator::
reduce(const void* sendBuf, void* recBuf, int count,
	   DataType type, ReduceOperation op, int rootProc) const
{
	PCL_PROFILE(pcl_ProcCom_reduce);
	if(is_local()) {memcpy(recBuf, sendBuf, count*GetSize(type)); return;}
	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::reduce: empty communicator.");

	MPI_Reduce(const_cast<void*>(sendBuf), recBuf, count, type, op, rootProc, m_comm->m_mpiComm);
}


size_t ProcessCommunicator::
reduce(const size_t &t, pcl::ReduceOperation op, int rootProc) const
{
	if(is_local()) return t;
	uint64 ret, tt = (uint64)t;
	reduce(&tt, &ret, 1, PCL_DT_UNSIGNED_LONG_LONG, op, rootProc);
	return (size_t)ret;
}

void
ProcessCommunicator::
allreduce(const void* sendBuf, void* recBuf, int count,
		  DataType type, ReduceOperation op) const
{
	PCL_PROFILE(pcl_ProcCom_allreduce);
	if(is_local()) {memcpy(recBuf, sendBuf, count*GetSize(type)); return;}
	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::allreduce: empty communicator.");

	MPI_Allreduce(const_cast<void*>(sendBuf), recBuf, count, type, op, m_comm->m_mpiComm);
}

size_t ProcessCommunicator::
allreduce(const size_t &t, pcl::ReduceOperation op) const
{
	if(is_local()) return t;
	uint64 ret, tt = (uint64)t;
	allreduce(&tt, &ret, 1, PCL_DT_UNSIGNED_LONG_LONG, op);
	return (size_t)ret;
}

void
ProcessCommunicator::
gather(const void* sendBuf, int sendCount, DataType sendType,
	   void* recBuf, int recCount, DataType recType, int root) const
{
	PCL_PROFILE(pcl_ProcCom_gather);
	if(is_local()) {memcpy(recBuf, sendBuf, recCount*GetSize(recType)); return;}

	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::gather: empty communicator.");
	
	MPI_Gather(const_cast<void*>(sendBuf), sendCount, sendType, recBuf,
			   recCount, recType, root, m_comm->m_mpiComm);
}


void ProcessCommunicator::
gather(ug::BinaryBuffer &buf, int root) const
{
	if(is_local()) return;

	int localSize;
	localSize = (int)buf.write_pos();
	
	if(pcl::ProcRank() == root) {
		std::vector<int> recvSizes(size());

		gather(&localSize, 1, PCL_DT_INT,
		       &recvSizes.front(), 1, PCL_DT_INT, root);

		std::vector<int> displacements(size());
		size_t totalSize = 0;
		for(size_t i = 0; i < recvSizes.size(); ++i){
			displacements[i] = (int)totalSize;
			totalSize += (size_t)recvSizes[i];
		}

		buf.reserve(totalSize);
		gatherv(MPI_IN_PLACE, localSize, PCL_DT_CHAR,
		        buf.buffer(), &recvSizes.front(), &displacements.front(),
		        PCL_DT_CHAR, root);

		buf.set_write_pos(totalSize);
	}
	else{
		gather(&localSize, 1, PCL_DT_INT,
		       NULL, 1, PCL_DT_INT, root);

		gatherv(buf.buffer(), localSize, PCL_DT_CHAR,
		        NULL, NULL, NULL,
		        PCL_DT_CHAR, root);
	}
}


void ProcessCommunicator::
scatter(const void* sendBuf, int sendCount, DataType sendType,
		 void* recBuf, int recCount, DataType recType, int root) const
{
	PCL_PROFILE(pcl_ProcCom_scatter);
	if(is_local()) {memcpy(recBuf, sendBuf, recCount*GetSize(recType)); return;}

	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::scatter: empty communicator.");
	
	MPI_Scatter(const_cast<void*>(sendBuf), sendCount, sendType, recBuf,
			   recCount, recType, root, m_comm->m_mpiComm);
}

void
ProcessCommunicator::
gatherv(const void* sendBuf, int sendCount, DataType sendType,
		void* recBuf, int* recCounts, int* displs,
		DataType recType, int root) const
{
	PCL_PROFILE(pcl_ProcCom_gatherv);
	if(is_local()) {memcpy(recBuf, sendBuf, displs[0] + recCounts[0]*GetSize(recType)); return;}

	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::gather: empty communicator.");

	MPI_Gatherv(const_cast<void*>(sendBuf), sendCount, sendType, recBuf,
				recCounts, displs, recType, root, m_comm->m_mpiComm);
}

void
ProcessCommunicator::
allgather(const void* sendBuf, int sendCount, DataType sendType,
		  void* recBuf, int recCount, DataType recType) const
{
	PCL_PROFILE(pcl_ProcCom_allgather);
	if(is_local()) {memcpy(recBuf, sendBuf, recCount*GetSize(recType)); return;}

	UG_COND_THROW(empty(), "ERROR in ProcessCommunicator::allgather: empty communicator.");
	
	MPI_Allgather(const_cast<void*>(sendBuf), sendCount, sendType, recBuf,
				  recCount, recType, m_comm->m_mpiComm);
}

void
ProcessCommunicator::
allgather(ug::BinaryBuffer &buf) const
{
	gather(buf, 0);
	broadcast(buf, 0);
}

void
ProcessCommunicator::
allgatherv(const void* sendBuf, int sendCount, DataType sendType,
			void* recBuf, int* recCounts, int* displs,
			DataType recType) const
{
	PCL_PROFILE(pcl_ProcCom_allgatherv);
	if(is_local()) {memcpy(recBuf, sendBuf, displs[0] + recCounts[0]*GetSize(recType)); return;}

	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::allgatherv: empty communicator.");
	
	MPI_Allgatherv(const_cast<void*>(sendBuf), sendCount, sendType, recBuf,
				   recCounts, displs, recType, m_comm->m_mpiComm);
}

void
ProcessCommunicator::
alltoall(const void* sendBuf, int sendCount, DataType sendType,
    void* recBuf, int recCount, DataType recType)
{
	PCL_PROFILE(pcl_ProcCom_alltoall);
	if(is_local()) {memcpy(recBuf, sendBuf, recCount*GetSize(recType)); return;}

	UG_COND_THROW(empty(), "ERROR in ProcessCommunicator::alltoall: empty communicator.");

	MPI_Alltoall(const_cast<void*>(sendBuf), sendCount, sendType, recBuf, recCount, recType, m_comm->m_mpiComm);
}

void
ProcessCommunicator::
send_data(void* pBuffer, int bufferSize, int destProc, int tag) const
{
	if(is_local()) return;
	PCL_PROFILE(pcl_ProcCom_send_data);

	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::send_data: empty communicator.");
	
	MPI_Request request;
	
	MPI_Isend(pBuffer, bufferSize, MPI_UNSIGNED_CHAR, destProc, 
			  tag, m_comm->m_mpiComm, &request);
	pcl::MPI_Wait(&request);
}

void
ProcessCommunicator::
send_data(void* pBuffer, int* pBufferSegSizes,
		  int* pRecProcMap, int numRecProcs, int tag) const
{
	if(is_local()) return;
	PCL_PROFILE(pcl_ProcCom_send_data__to_many);

	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::send_data: empty communicator.");

//	send data
	std::vector<MPI_Request> vSendRequests(numRecProcs);
		
	for(int i = 0; i < numRecProcs; ++i)
	{
		MPI_Isend(pBuffer, pBufferSegSizes[i], MPI_UNSIGNED_CHAR,
				  pRecProcMap[i], tag, m_comm->m_mpiComm, &vSendRequests[i]);
		pBuffer = (byte*)pBuffer + pBufferSegSizes[i];
	}
	
//	wait until data has been received
	Waitall(vSendRequests);
}

void
ProcessCommunicator::
receive_data(void* pBuffOut, int bufferSize, int srcProc, int tag) const
{
	if(is_local()) return;
	PCL_PROFILE(pcl_ProcCom_recv_data);

	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::receive_data: empty communicator.");
	
	MPI_Request request;
	
	MPI_Irecv(pBuffOut, bufferSize, MPI_UNSIGNED_CHAR,	
					srcProc, tag, m_comm->m_mpiComm, &request);					
	pcl::MPI_Wait(&request);
}

void ProcessCommunicator::
distribute_data(void* recvBufOut, int* recvBufSegSizesOut,
				int* recvFromRanks, int numRecvs,
				void* sendBuf, int* sendBufSegSizes,
			  	int* sendToRanks, int numSends, int tag) const
{
	if(is_local()) return;
	PCL_PROFILE(pcl_ProcCom_distribute_data);

	UG_COND_THROW(empty(),	"ERROR in ProcessCommunicator::distribute_data: empty communicator.");

//	used for mpi-communication.
	std::vector<MPI_Request> vSendRequests(numSends);
	std::vector<MPI_Request> vReceiveRequests(numRecvs);
	
//	wait until data has been received
	std::vector<MPI_Status> vSendStates(numSends);
	std::vector<MPI_Status> vReceiveStates(numRecvs);

//	shedule receives first
	for(int i = 0; i < numRecvs; ++i)
	{
		MPI_Irecv(recvBufOut, recvBufSegSizesOut[i], MPI_UNSIGNED_CHAR,	
				  recvFromRanks[i], tag, m_comm->m_mpiComm,
				  &vReceiveRequests[i]);
		recvBufOut = (byte*)recvBufOut + recvBufSegSizesOut[i];
	}

//	now send the data
	for(int i = 0; i < numSends; ++i)
	{
		MPI_Isend(sendBuf, sendBufSegSizes[i], MPI_UNSIGNED_CHAR,
				  sendToRanks[i], tag, m_comm->m_mpiComm,
				  &vSendRequests[i]);
		sendBuf = (byte*)sendBuf + sendBufSegSizes[i];
	}

	Waitall(vReceiveRequests, vSendRequests);	
}

void ProcessCommunicator::
distribute_data(ug::BinaryBuffer& recvBufOut, int* segSizesOut,
				int* recvFromRanks, int numRecvFroms,
				void* sendBuf, int* sendSegSizes,
				int* sendToRanks, int numSendTos, int tag) const
{
	if(is_local()) return;
	PCL_PROFILE(pcl_ProcCom_distribute_data__flex);

//	small helper arrays
	vector<int> tmpRecvSegSizes(numRecvFroms, sizeof(int));
	vector<int> tmpSendSegSizes(numSendTos, sizeof(int));

//	every process receives the size of the data-buffer first.
	vector<int> bufferSizes(numRecvFroms);

//	exchange buffer sizes (use an arbitrary tag)
	distribute_data(GetDataPtr(bufferSizes), GetDataPtr(tmpRecvSegSizes),
					recvFromRanks, numRecvFroms,
					sendSegSizes, GetDataPtr(tmpSendSegSizes),
					sendToRanks, numSendTos, tag);

//	calculate buffer sizes and resize the binary stream
	int totalSize = 0;
	for(int i = 0; i < numRecvFroms; ++i){
		totalSize += bufferSizes[i];
		segSizesOut[i] = bufferSizes[i];
	}

	recvBufOut.reserve(totalSize);

//	now exchange the buffers
	distribute_data(recvBufOut.buffer(), GetDataPtr(bufferSizes),
					recvFromRanks, numRecvFroms,
					sendBuf, sendSegSizes,
					sendToRanks, numSendTos, tag);
	recvBufOut.set_write_pos(totalSize);
}

void ProcessCommunicator::
distribute_data(ug::BinaryBuffer* recvBufs, int* recvFromRanks, int numRecvs,
				ug::BinaryBuffer* sendBufs, int* sendToRanks, int numSends,
				int tag) const
{
	if(is_local()) return;
	PCL_PROFILE(pcl_ProcCom_distribute_data__multi_buf);

//	small helper arrays
	vector<int> tmpRecvSegSizes(numRecvs, sizeof(int));
	vector<int> tmpSendSegSizes(numSends, sizeof(int));

//	the actual data sizes which will be sent and received
	vector<int> recvSizes(numRecvs);
	vector<int> sendSizes(numSends);
	for(int i = 0; i < numSends; ++i)
		sendSizes[i] = sendBufs[i].write_pos();

//	exchange buffer sizes (use an arbitrary tag)
	distribute_data(GetDataPtr(recvSizes), GetDataPtr(tmpRecvSegSizes),
					recvFromRanks, numRecvs,
					GetDataPtr(sendSizes), GetDataPtr(tmpSendSegSizes),
					sendToRanks, numSends, tag);

//	used for mpi-communication.
	std::vector<MPI_Request> vSendRequests(numSends);
	std::vector<MPI_Request> vReceiveRequests(numRecvs);
	
//	wait until data has been received
	std::vector<MPI_Status> vSendStates(numSends);
	std::vector<MPI_Status> vReceiveStates(numRecvs);

//	shedule receives first
	for(int i = 0; i < numRecvs; ++i){
		recvBufs[i].clear();
		recvBufs[i].reserve(recvSizes[i]);
		MPI_Irecv(recvBufs[i].buffer(), recvSizes[i], MPI_UNSIGNED_CHAR,	
				  recvFromRanks[i], tag, m_comm->m_mpiComm,
				  &vReceiveRequests[i]);
	}

//	now send the data
	for(int i = 0; i < numSends; ++i){
		MPI_Isend(sendBufs[i].buffer(), sendSizes[i], MPI_UNSIGNED_CHAR,
				  sendToRanks[i], tag, m_comm->m_mpiComm,
				  &vSendRequests[i]);
	}

	Waitall(vReceiveRequests, vSendRequests);

//	adjust write-pos in receive buffers
	for(int i = 0; i < numRecvs; ++i)
		recvBufs[i].set_write_pos(recvSizes[i]);
}

void ProcessCommunicator::
barrier() const
{
	PCL_PROFILE(pcl_ProcCom_barrier);
	if(is_local()) return;
	MPI_Barrier(m_comm->m_mpiComm);
}

void ProcessCommunicator::broadcast(void *v, size_t size, DataType type, int root) const
{
	PCL_PROFILE(pcl_ProcCom_Bcast);
	if(is_local()) return;
	//UG_LOG("broadcasting " << (root==pcl::ProcRank() ? "(sender) " : "(receiver) ") << size << " root = " << root << "\n");
	MPI_Bcast(v, size, type, root, m_comm->m_mpiComm);
}

void ProcessCommunicator::broadcast(ug::BinaryBuffer &buf, int root) const
{
	if(is_local()) return;
	if(pcl::ProcRank() == root)
	{
		long size;
		size = buf.write_pos();
		broadcast(&size, 1, PCL_DT_LONG, root);
		broadcast(buf.buffer(), size, PCL_DT_CHAR, root);
	}
	else
	{
		long size;
		size = buf.write_pos();
		broadcast(&size, 1, PCL_DT_LONG, root);
		buf.reserve(size);
		broadcast(buf.buffer(), size, PCL_DT_CHAR, root);
		buf.set_write_pos(size);
	}
}

void ProcessCommunicator::broadcast(size_t &s, int root) const
{
	if(is_local()) return;
	unsigned long l = s;
	broadcast<unsigned long>(l);
	s = l;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
ProcessCommunicator::CommWrapper::
CommWrapper() :
	m_mpiComm(PCL_COMM_WORLD),
	m_bReleaseCommunicator(false)
{}

ProcessCommunicator::CommWrapper::
CommWrapper(const MPI_Comm& comm, bool bReleaseComm) :
	m_mpiComm(comm),
	m_bReleaseCommunicator(bReleaseComm)
{}

ProcessCommunicator::CommWrapper::
~CommWrapper()
{
	if(m_bReleaseCommunicator)
		MPI_Comm_free(&m_mpiComm);
}


std::ostream &operator << (std::ostream &out, const ProcessCommunicator &processCommunicator)
{
	out << "ProcessCommunicator ";
	if(processCommunicator.is_local()) out << "LOCAL";
	else if(processCommunicator.empty()) out << "EMPTY";
	else
	{
		if(processCommunicator.is_world()) out << "(WORLD) ";
		out << "size = " << processCommunicator.size();
		out << ". involved procs: [ ";
		for(size_t i=0; i<processCommunicator.size(); i++)
			out << processCommunicator.get_proc_id(i) << " ";
		out << "]";
	}
	return out;
}


}//	end of namespace pcl
