/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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


#include <mpi.h>
#include "pcl_comm_world.h"
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
