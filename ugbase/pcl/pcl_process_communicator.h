//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d06

#ifndef __H__PCL__PCL_PROCESS_COMMUNICATOR__
#define __H__PCL__PCL_PROCESS_COMMUNICATOR__

#include <iostream>
#include <map>
#include "pcl_methods.h"
#include "common/smart_pointer.h"

namespace pcl
{

///	values that can be passed to a ProcessCommunicators constructor.
enum ProcessCommunicatorDefaults
{
	PCD_EMPTY = 0,
	PCD_WORLD = 1
};


/**
 * A ProcessCommunicator is a very lightweight object that can be passed
 * by value. Creation using the constructor is a lightweight operation too.
 * Creating a new communicator using create_sub_communicator however requires
 * communication and thus shouldn't be unnecessarily called.
 *
 * Please note that ProcessCommunicators created on different processes via
 * create_sub_communicator(...) should be deleted at the same time - even
 * if the creating processes are not part of the communicator.
 */
class ProcessCommunicator
{
	public:
	///	creates a communicator.
	/**	By default a communicator for all processes is generated.*/
		ProcessCommunicator(ProcessCommunicatorDefaults pcd = PCD_WORLD);
		
	///	returns true if the communicator is empty, false if not.
		inline bool empty() const		{return m_comm->m_mpiComm == MPI_COMM_NULL;}
		
	///	returns the size of the communicator
		size_t size() const;
		
	///	creates a new communicator containing a subset of the current communicator
	/**	Note that this method has to be called by all processes in the current
	 *	communicator - even if they don't want to participate in the new one.*/
		ProcessCommunicator create_sub_communicator(bool participate);
		
	///	performs MPI_Allreduce on the processes of the communicator.
		void allreduce(void* sendBuf, void* recBuf, int count,
					   DataType type, ReduceOperation op) const;
		
	///	performs MPI_Allgather on the processes of the communicator.
	/**
	 * \param sendBuf starting address of send buffer (choice)
	 * \param sendCount number of elements in send buffer (integer)
	 * \param sendType data type of send buffer elements (handle)
	 * \param recCount number of elements received from any process (integer)
	 * \param recType data type of receive buffer elements (handle)
	 */
		void allgather(void* sendBuf, int sendCount, DataType sendType,
					   void* recBuf, int recCount, DataType recType) const;

	///	performs MPI_Allgatherv on the processes of the communicator.
	/**
	 * \param sendBuf 	starting address of send buffer (choice)
	 * \param sendCount 	number of elements in send buffer (integer)
	 * \param sendType 	data type of send buffer elements (handle)
	 * \param recCounts 	integer array (of length group size) containing the number of elements that are received from each process
	 * \param displs 	integer array (of length group size). Entry i specifies the displacement (relative to recvbuf ) at which to place the incoming data from process i
	 * \param recType 	data type of receive buffer elements (handle)
	 */
		void allgatherv(void* sendBuf, int sendCount, DataType sendType,
						void* recBuf, int* recCounts, int* displs,
						DataType recType) const;
	private:
	///	holds an mpi-communicator.
	/**	A variable stores whether the communicator has to be freed when the
	 *	the wrapper is deleted.*/
		struct CommWrapper{
			CommWrapper(const MPI_Comm& comm,
						bool bReleaseComm);
			~CommWrapper();
			
			MPI_Comm	m_mpiComm;
			bool		m_bReleaseCommunicator;
		};
		
	///	Smart-pointer that encapsulates a CommWrapper.
		typedef SmartPtr<CommWrapper> SPCommWrapper;
		
	private:
	///	smart-pointer to an instance of a CommWrapper.
		SPCommWrapper	m_comm;
};

}//	end of namespace pcl


#endif
