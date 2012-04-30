//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m05 d08

#ifndef __H__PCL_METHODS__
#define __H__PCL_METHODS__

#include <vector>
#include <list>
#include <iostream>

//	Don't rely on mpi being included.
//	It is only included to allow us to define some constants.
//	This include will most likely be removed in future versions.
#include "mpi.h"
#include "pcl_base.h"
#include "common/types.h"
#include "pcl_datatype.h"

namespace pcl
{

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
 
}//	end of namespace

#endif
