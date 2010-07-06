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
		inline bool empty()		{return m_comm->m_mpiComm == MPI_COMM_NULL;}
		
	///	creates a new communicator containing a subset of the current communicator
	/**	Note that this method has to be called by all processes in the current
	 *	communicator - even if they don't want to participate in the new one.*/
		ProcessCommunicator create_sub_communicator(bool participate);
		
	///	performs MPI_Allreduce on the processes of the communicator.
		void allreduce(void* sendBuf, void* recBuf, int count,
					   DataType type, ReduceOperation op);
		
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
