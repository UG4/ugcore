//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m06 d07

#include "mpi.h"
#include "pcl.h"

static int OUTPUT_PROC_RANK = -1;

namespace pcl
{

////////////////////////////////////////////////////////////////////////
void Init(int argc, char* argv[])
{
//	init mpi
	MPI_Init(&argc, &argv);
}

////////////////////////////////////////////////////////////////////////
void Finalize()
{
	MPI_Finalize();
}

////////////////////////////////////////////////////////////////////////
ProcID GetProcRank()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
}

int GetNumProcesses()
{
	int numProcesses;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
	return numProcesses;
}

////////////////////////////////////////////////////////////////////////
ProcID GetOutputProcRank()
{
	return OUTPUT_PROC_RANK;
}

////////////////////////////////////////////////////////////////////////
void SetOutputProcRank(ProcID rank)
{
	OUTPUT_PROC_RANK = rank;
}

////////////////////////////////////////////////////////////////////////
double Time()
{
	return MPI_Wtime();
}

////////////////////////////////////////////////////////////////////////
//
void SendData(ProcID destProc, void* pBuffer, int bufferSize, int tag)
{
	MPI_Request request;
	MPI_Status	status;
	
	MPI_Isend(pBuffer, bufferSize, MPI_UNSIGNED_CHAR, destProc, tag, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
}

////////////////////////////////////////////////////////////////////////
//
void ReceiveData(void* pBuffOut, ProcID srcProc, int bufferSize, int tag)
{
	MPI_Request request;
	MPI_Status	status;
	
	MPI_Irecv(pBuffOut, bufferSize, MPI_UNSIGNED_CHAR,	
					srcProc, tag, MPI_COMM_WORLD, &request);
					
	MPI_Wait(&request, &status);
}

////////////////////////////////////////////////////////////////////////
//
void CollectData(ProcID thisProcID, int firstSendProc, int numSendProcs,
					void* pBuffer, int bufferSizePerProc, int tag)
{
//	receive data
	std::vector<MPI_Request> vReceiveRequests(numSendProcs);
		
	int srcProcIndex = firstSendProc;
	for(int i = 0; i < numSendProcs; ++i, ++srcProcIndex)
	{
		if(srcProcIndex == thisProcID)
			srcProcIndex++;
	
		MPI_Irecv((byte*)pBuffer + bufferSizePerProc * i, bufferSizePerProc, MPI_UNSIGNED_CHAR,	
					srcProcIndex, tag, MPI_COMM_WORLD, &vReceiveRequests[i]);
	}
	
//	wait until data has been received
	std::vector<MPI_Status> vReceiveStatii(numSendProcs);//TODO: fix spelling!
	MPI_Waitall(numSendProcs, &vReceiveRequests[0], &vReceiveStatii[0]);
}

////////////////////////////////////////////////////////////////////////
//
void DistributeData(ProcID thisProcID, int firstRecProc, int numRecProcs,
					void* pBuffer, int* pBufferSegSizes, int tag)
{
//	receive data
	std::vector<MPI_Request> vSendRequests(numRecProcs);
		
	int destProcIndex = firstRecProc;

	for(int i = 0; i < numRecProcs; ++i, ++destProcIndex)
	{
		if(destProcIndex == thisProcID)
			destProcIndex++;
	
		MPI_Isend(pBuffer, pBufferSegSizes[i], MPI_UNSIGNED_CHAR, destProcIndex, tag, MPI_COMM_WORLD, &vSendRequests[i]);
		pBuffer = (byte*)pBuffer + pBufferSegSizes[i];
	}
	
//	wait until data has been received
	std::vector<MPI_Status> vSendStatii(numRecProcs);//TODO: fix spelling!
	MPI_Waitall(numRecProcs, &vSendRequests.front(), &vSendStatii.front());
}

////////////////////////////////////////////////////////////////////////
//
void DistributeData(ProcID thisProcID, int* pRecProcMap, int numRecProcs,
					void* pBuffer, int* pBufferSegSizes, int tag)
{
//	receive data
	std::vector<MPI_Request> vSendRequests(numRecProcs);
		
	for(int i = 0; i < numRecProcs; ++i)
	{
		MPI_Isend(pBuffer, pBufferSegSizes[i], MPI_UNSIGNED_CHAR, pRecProcMap[i], tag, MPI_COMM_WORLD, &vSendRequests[i]);
		pBuffer = (byte*)pBuffer + pBufferSegSizes[i];
	}
	
//	wait until data has been received
	std::vector<MPI_Status> vSendStatii(numRecProcs);//TODO: fix spelling!
	MPI_Waitall(numRecProcs, &vSendRequests.front(), &vSendStatii.front());
}

////////////////////////////////////////////////////////////////////////
//
void AllReduce(void* sendBuf, void* recBuf, int count, DataType type,
				ReduceOperation op)
{
	MPI_Allreduce(sendBuf, recBuf, count, type, op, MPI_COMM_WORLD);
}

}//	end of namespace
