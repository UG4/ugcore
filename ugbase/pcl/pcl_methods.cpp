//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m06 d07


#include <mpi.h>
#include "pcl_comm.h"
#include "pcl_methods.h"
#include "common/log.h"
#include "pcl_profiling.h"

namespace pcl
{

////////////////////////////////////////////////////////////////////////
double Time()
{
	return MPI_Wtime();
}

////////////////////////////////////////////////////////////////////////
//
void SendData(ProcID destProc, void* pBuffer, int bufferSize, int tag)
{
	UG_LOG("DEPRECIATED: pcl::SendData\n");

	MPI_Request request;
	
	MPI_Isend(pBuffer, bufferSize, MPI_UNSIGNED_CHAR, destProc, tag, PCL_COMM_WORLD, &request);
	pcl::MPI_Wait(&request);
}

////////////////////////////////////////////////////////////////////////
//
void ReceiveData(void* pBuffOut, ProcID srcProc, int bufferSize, int tag)
{
	UG_LOG("DEPRECIATED: pcl::ReceiveData\n");

	MPI_Request request;
	
	MPI_Irecv(pBuffOut, bufferSize, MPI_UNSIGNED_CHAR,
					srcProc, tag, PCL_COMM_WORLD, &request);

	pcl::MPI_Wait(&request);
}

////////////////////////////////////////////////////////////////////////
//
void CollectData(ProcID thisProcID, int firstSendProc, int numSendProcs,
					void* pBuffer, int bufferSizePerProc, int tag)
{
	UG_LOG("DEPRECIATED: pcl::CollectData\n");

//	receive data
	std::vector<MPI_Request> vReceiveRequests(numSendProcs);
		
	int srcProcIndex = firstSendProc;
	for(int i = 0; i < numSendProcs; ++i, ++srcProcIndex)
	{
		if(srcProcIndex == thisProcID)
			srcProcIndex++;
	
		MPI_Irecv((byte*)pBuffer + bufferSizePerProc * i, bufferSizePerProc, MPI_UNSIGNED_CHAR,	
					srcProcIndex, tag, PCL_COMM_WORLD, &vReceiveRequests[i]);
	}
	
//	wait until data has been received
	Waitall(vReceiveRequests);
}

////////////////////////////////////////////////////////////////////////
//
void DistributeData(ProcID thisProcID, int firstRecProc, int numRecProcs,
					void* pBuffer, int* pBufferSegSizes, int tag)
{
	UG_LOG("DEPRECIATED: pcl::DistributeData\n");

//	receive data
	std::vector<MPI_Request> vSendRequests(numRecProcs);
		
	int destProcIndex = firstRecProc;

	for(int i = 0; i < numRecProcs; ++i, ++destProcIndex)
	{
		if(destProcIndex == thisProcID)
			destProcIndex++;
	
		MPI_Isend(pBuffer, pBufferSegSizes[i], MPI_UNSIGNED_CHAR, destProcIndex, tag, PCL_COMM_WORLD, &vSendRequests[i]);
		pBuffer = (byte*)pBuffer + pBufferSegSizes[i];
	}
	
//	wait until data has been received
	Waitall(vSendRequests);
}

////////////////////////////////////////////////////////////////////////
//
void DistributeData(ProcID thisProcID, int* pRecProcMap, int numRecProcs,
					void* pBuffer, int* pBufferSegSizes, int tag)
{
	UG_LOG("DEPRECIATED: pcl::DistributeData\n");

//	receive data
	std::vector<MPI_Request> vSendRequests(numRecProcs);
		
	for(int i = 0; i < numRecProcs; ++i)
	{
		MPI_Isend(pBuffer, pBufferSegSizes[i], MPI_UNSIGNED_CHAR, pRecProcMap[i], tag, PCL_COMM_WORLD, &vSendRequests[i]);
		pBuffer = (byte*)pBuffer + pBufferSegSizes[i];
	}
	
//	wait until data has been received	
	Waitall(vSendRequests);
}

////////////////////////////////////////////////////////////////////////
//
void AllReduce(void* sendBuf, void* recBuf, int count, DataType type,
				ReduceOperation op)
{
	UG_LOG("DEPRECIATED: pcl::AllReduce\n");
	MPI_Allreduce(sendBuf, recBuf, count, type, op, PCL_COMM_WORLD);
}

}//	end of namespace
