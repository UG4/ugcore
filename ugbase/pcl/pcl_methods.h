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

#ifndef __H__PCL_METHODS__
#define __H__PCL_METHODS__

#include <vector>
#include <list>
#include <iostream>

//	Don't rely on mpi being included.
//	It is only included to allow us to define some constants.
//	This include will most likely be removed in future versions.
#include <mpi.h>
#include "pcl_comm_world.h"
#include "pcl_base.h"
#include "common/types.h"
#include "common/profiler/profiler.h"
#include "pcl_datatype.h"
#include <cassert>


namespace pcl
{

/// \addtogroup pcl
/// \{

typedef int ProcID;

//	ReduceOperation
#define PCL_RO_MAX 		MPI_MAX
#define PCL_RO_MIN 		MPI_MIN
#define PCL_RO_SUM 		MPI_SUM
#define PCL_RO_PROD 	MPI_PROD
#define PCL_RO_LAND 	MPI_LAND
#define PCL_RO_BAND 	MPI_BAND
#define PCL_RO_LOR 		MPI_LOR
#define PCL_RO_BOR 		MPI_BOR
#define PCL_RO_LXOR 	MPI_LXOR
#define PCL_RO_BXOR 	MPI_BXOR
#define PCL_RO_MAXLOC 	MPI_MAXLOC
#define PCL_RO_MINLOC 	MPI_MINLOC

typedef MPI_Op ReduceOperation;



////////////////////////////////////////////////////////////////////////

///	returns the time in seconds
double Time();

///	sends data to another process. data may be received using \sa ReceiveData or \sa CollectData.
void SendData(ProcID destProc, void* pBuffer, int bufferSize, int tag);

///	receives the data that was send with \sa SendData or \sa DistributeData.
void ReceiveData(void* pBuffOut, ProcID srcProc, int bufferSize, int tag);

///	collect the data send with send_data from proc firstSendProc to numSendProcs excluding destProc.
void CollectData(ProcID thisProcID, int firstSendProc, int numSendProcs,
					void* pBuffer, int bufferSizePerProc, int tag);

///	sends the data in the different sections of the buffer to the specified processes.
/**
 * pBufferSegSizes has to contain numRecProcs elements.
 * Buffer-Segments are send to the processes
 * firstRecProc, ..., firstRecProc + numRecProcs.
 * If thisProcID lies in this range, then the buffer-segments are
 * sent to firstRecProc, ..., firstRecProc + numRecProcs + 1,
 * excluding thisProcID.
 */
void DistributeData(ProcID thisProcID, int firstRecProc, int numRecProcs,
					void* pBuffer, int* pBufferSegSizes, int tag);

///	sends the data in the different sections of the buffer to the specified processes.
/**
 * pRecProcMap and pBufferSegSizes have to contain numRecProcs elements.
 * Entries in pRecProcMap specify the target-process of the i-th buffer
 * segment.
 */
void DistributeData(ProcID thisProcID, int* pRecProcMap, int numRecProcs,
					void* pBuffer, int* pBufferSegSizes, int tag);

///	reduces the data to a single buffer using the specified ReduceOperation and distributes the result to all processes.
/**
 * Both, sendBuf and recBuf have to hold count elements of the specified type.
 * \param sendBuf	Sending buffer
 * \param recBuf	Recieve buffer
 * \param count
 * \param type		Data type
 * \param op has to be one of the
 */
void AllReduce(void* sendBuf, void* recBuf, int count, DataType type,
				ReduceOperation op);

//void StartWait(), StopWait();

inline void MPI_Waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses)
{
//	StartWait();
	PROFILE_FUNC_GROUP("mpi");
	::MPI_Waitall(count, array_of_requests, array_of_statuses);
//	StopWait();
}

inline void Waitall(std::vector<MPI_Request> &requests, std::vector<MPI_Status> &statuses)
{
//	StartWait();
	PROFILE_FUNC_GROUP("mpi");
	assert(requests.size() == statuses.size());
	pcl::MPI_Waitall(requests.size(), &requests[0], &statuses[0]);
//	StopWait();
}

inline void Waitall(std::vector<MPI_Request> &requests)
{	
//	StartWait();
	PROFILE_FUNC_GROUP("mpi");
	if(requests.size() > 0) 
		pcl::MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
//	StopWait();
}

inline void Waitall(std::vector<MPI_Request> &requests, std::vector<MPI_Request> &requests2)
{	
//	StartWait();
	PROFILE_FUNC_GROUP("mpi");
	if(requests.size() > 0) pcl::MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
	if(requests2.size() > 0) pcl::MPI_Waitall(requests2.size(), &requests2[0], MPI_STATUSES_IGNORE);
//	StopWait();
}

inline int MPI_Wait(MPI_Request *request, MPI_Status *status=MPI_STATUS_IGNORE)
{
//	StartWait();
	PROFILE_FUNC_GROUP("mpi");
	int i=::MPI_Wait(request, status);
//	StopWait();
	return i;
}



/*inline int Irecv(void *buf, int count, MPI_Datatype datatype, int source,
                             int tag, MPI_Comm comm, MPI_Request *request)
{
	PROFILE_FUNC_GROUP("mpi");
	return ::MPI_Irecv(buf, count, datatype, source, tag, comm, request);
}

template<typename T>
inline int IRecv(T *buf, int count, int source, int tag, MPI_Comm comm, MPI_Request *request)
{
	PROFILE_FUNC_GROUP("mpi");
	return ::MPI_IRecv(buf, count, DataTypeTraits<T>::get_data_type(), source, tag, comm, request);
}*/


// end group pcl
/// \}


}//	end of namespace

#endif
