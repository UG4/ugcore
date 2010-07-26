//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d06

#include <vector>
#include <cassert>
#include "pcl_methods.h"
#include "common/smart_pointer.h"
#include "pcl_process_communicator.h"
#include "common/log.h"
using namespace std;

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
			m_comm = SPCommWrapper(new CommWrapper(MPI_COMM_WORLD, false));
			break;
	}
}

size_t ProcessCommunicator::
size()
{
	if(m_comm->m_mpiComm == MPI_COMM_NULL)
		return 0;
	
	int size;
	if(MPI_Comm_size(m_comm->m_mpiComm, &size) == MPI_SUCCESS)
		return (size_t)size;
		
	UG_LOG("  ERROR in ProcessCommunicator::size(): Unknown MPI Error. Returning 0.\n");
	return 0;
}

ProcessCommunicator
ProcessCommunicator::
create_sub_communicator(bool participate)
{
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
	if(participate){
		int rank;
		MPI_Comm_rank(m_comm->m_mpiComm, &rank);
		srcArray[rank] = 1;
	}

//	synchronize the newProcs array between all processes in the communicator
	vector<int> destArray(size, 0);

	allreduce(&srcArray.front(), &destArray.front(),
			   size, PCL_DT_INT, PCL_RO_MAX);

//	build a local array that holds all the procs that shall go
//	into the new communicator
	vector<int> newProcs;
	newProcs.reserve(size);

	for(size_t i = 0; i < destArray.size(); ++i){
		if(destArray[i])
			newProcs.push_back(i);
	}

//	if newProcs is not empty, we'll build a new mpi-communicator.
	if(newProcs.size() == 0)
		return ProcessCommunicator(PCD_EMPTY);

	MPI_Group grpOld;
	MPI_Group grpNew;
	MPI_Comm commNew;

	MPI_Comm_group(m_comm->m_mpiComm, &grpOld);
	MPI_Group_incl(grpOld, (int)newProcs.size(), &newProcs.front(), &grpNew);
	MPI_Comm_create(m_comm->m_mpiComm, grpNew, &commNew);

//	create a new ProcessCommunicator
//	if the process is not participating, MPI_Comm_create will return MPI_COMM_NULL
	if(commNew == MPI_COMM_NULL)
		return ProcessCommunicator(PCD_EMPTY);

//	the process participates - create the ProcessCommunicator
	ProcessCommunicator newProcComm;
	newProcComm.m_comm = SPCommWrapper(new CommWrapper(commNew, true));
	return newProcComm;
}

void
ProcessCommunicator::
allreduce(void* sendBuf, void* recBuf, int count,
		  DataType type, ReduceOperation op)
{
	assert(!empty() &&
			"ERROR in ProcessCommunicator::allreduce: can't perform all_reduce on empty communicator.");
	if(empty()){
		UG_LOG("ERROR in ProcessCommunicator::allreduce: empty communicator.\n");
	}
	
	MPI_Allreduce(sendBuf, recBuf, count, type, op, m_comm->m_mpiComm);
}

void
ProcessCommunicator::
allgather(void* sendBuf, int sendCount, DataType sendType,
		  void* recBuf, int recCount, DataType recType)
{
	assert(!empty() &&
			"ERROR in ProcessCommunicator::allreduce: can't perform all_reduce on empty communicator.");
	if(empty()){
		UG_LOG("ERROR in ProcessCommunicator::allreduce: empty communicator.\n");
	}
	
	MPI_Allgather(sendBuf, sendCount, sendType, recBuf,
				  recCount, recType, m_comm->m_mpiComm);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
ProcessCommunicator::CommWrapper::
CommWrapper(const MPI_Comm& comm,
			bool bReleaseComm) :
	m_mpiComm(comm),
	m_bReleaseCommunicator(bReleaseComm)
{}

ProcessCommunicator::CommWrapper::
~CommWrapper()
{
	if(m_bReleaseCommunicator)
		MPI_Comm_free(&m_mpiComm);
}

}//	end of namespace pcl
