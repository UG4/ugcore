//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d06

#ifndef __H__PCL__PCL_PROCESS_COMMUNICATOR__
#define __H__PCL__PCL_PROCESS_COMMUNICATOR__

#include <iostream>
#include <map>
#include <vector>
#include "pcl_methods.h"
#include "common/util/smart_pointer.h"

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
		
	///	returns the i-th process in the communicator
		int get_proc_id(size_t index) const;

	///	returns the proc-id relative to this communicator
	/**	This method has a worst time complexity of O(n)*/
		int get_local_proc_id() const;

	///	creates a new communicator containing a subset of the current communicator
	/**	Note that this method has to be called by all processes in the current
	 *	communicator - even if they don't want to participate in the new one.*/
		ProcessCommunicator create_sub_communicator(bool participate);
		
	///	performs MPI_Allreduce on the processes of the communicator.
	/**	This method synchronises involved processes.
	 */	
		void allreduce(void* sendBuf, void* recBuf, int count,
					   DataType type, ReduceOperation op) const;
		
	///	performs MPI_Gather on the processes of the communicator.
	/**	This method synchronises involved processes.
	 * \param sendBuf starting address of send buffer (choice)
	 * \param sendCount number of elements in send buffer (integer)
	 * \param sendType data type of send buffer elements (handle)
	 * \param recBuf only significant for root process. All gathered data is written here.
	 * \param recCount number of elements received from each process (integer)
	 * \param recType data type of receive buffer elements (handle)
	 * \param root The rank of the process that receives all the data.
	 */
		void gather(void* sendBuf, int sendCount, DataType sendType,
					void* recBuf, int recCount, DataType recType, int root) const;

	///	performs MPI_Gatherv on the processes of the communicator.
	/**	This method synchronises involved processes.
	 * \param sendBuf 	starting address of send buffer (choice)
	 * \param sendCount 	number of elements in send buffer (integer)
	 * \param sendType 	data type of send buffer elements (handle)
	 * \param recCounts 	integer array (of length group size) containing the
	 * 						number of elements that are received from each process.
	 * 						(Only significant at root)
	 * \param recBuf 		starting address of receive buffer (choice)
	 * 						(Only significant at root)
	 * \param displs 	integer array (of length group size).
	 * 					Entry i specifies the displacement (relative to recvbuf )
	 * 					at which to place the incoming data from process i.
	 * 					(Only significant at root)
	 * \param recType 	data type of receive buffer elements (handle).
	 * 					(Only significant at root)
	 * \param root The rank of the process that receives all the data.
	 */
		void gatherv(void* sendBuf, int sendCount, DataType sendType,
						void* recBuf, int* recCounts, int* displs,
						DataType recType, int root) const;

	///	performs MPI_Allgather on the processes of the communicator.
	/**	This method synchronises involved processes.
	 * \param sendBuf starting address of send buffer (choice)
	 * \param sendCount number of elements in send buffer (integer)
	 * \param sendType data type of send buffer elements (handle)
	 * \param recBuf starting address of receive buffer (choice)
	 * \param recCount number of elements received from any process (integer)
	 * \param recType data type of receive buffer elements (handle)
	 */
		void allgather(void* sendBuf, int sendCount, DataType sendType,
					   void* recBuf, int recCount, DataType recType) const;

	///	performs MPI_Allgatherv on the processes of the communicator.
	/**	This method synchronises involved processes.
	 * \param sendBuf 	starting address of send buffer (choice)
	 * \param sendCount 	number of elements in send buffer (integer)
	 * \param sendType 	data type of send buffer elements (handle)
	 * \param recCounts 	integer array (of length group size) containing the number of elements that are received from each process
	 * \param recBuf starting address of receive buffer (choice)
	 * \param displs 	integer array (of length group size). Entry i specifies the displacement (relative to recvbuf ) at which to place the incoming data from process i
	 * \param recType 	data type of receive buffer elements (handle)
	 */
		void allgatherv(void* sendBuf, int sendCount, DataType sendType,
						void* recBuf, int* recCounts, int* displs,
						DataType recType) const;

	///	sends data with the given tag to the specified process.
	/**	This method waits until the data has been sent.*/
		void send_data(void* pBuffer, int bufferSize, int destProc, int tag);
		
	///	sends data in different blocks of pBuffer to the processes in pRecProcMap.
	/**	This method synchronises involved processes.
	 *	Call receive_data on the processes in pRecProcMap to receive the sent data.
	 *
	 *	\param pBuffer: Blocks of data. The i-th block is send to the i-th process
	 *					of pRecProcMap.
	 *	\param pBufferSegSizes: The i-th entry holds the size of the i-th block in pBuffer.
	 *	\param pRecProcMap: The i-th entry holds the process-rank to which the i-th
	 *						block shall be send.
	 *	\param numRecProcs: The number of processes to which data shall be send.
	 *						Note that pBufferSegSizes and pRecProcMap have to have
	 *						numRecProcs entries.
	 *	\param tag: A tag that tags the message. Use the same tag in receive_data.
	 */
		void send_data(void* pBuffer, int* pBufferSegSizes,
					   int* pRecProcMap, int numRecProcs, int tag);
							 
	///	receives data from srcPrc with the specified tag.
	/**	This method waits until the data has been received
	 */
		void receive_data(void* pBuffOut, int bufferSize, int srcProc, int tag);

	///	sends and receives data to / from multiple processes.
	/**
	 * \param pBufferOut: Received data is written to this buffer.
	 *					  Make sure that this buffer is big enough (sum of all seg-sizes).
	 * \param pBufferOutSegSizes: i-th entry corresponds to the size of the i-th
	 *							  segment of pBufferOut.
	 * \param pSenderProcMap: The processes from which data is received.
	 * \param numSenderProcs: The number of processes from which data is received.
	 *						  Has to be the same as the size of pBufferOutSegSizes and
	 *						  pSenderProcMap.
	 *
	 * \param pBuffer: Holds the data that is to be send to other processes.
	 * \param pBufferSegSizes: i-th entry corresponds to the size of the i-th
	 *						   segment in pBuffer.
	 * \param pRecvProcMap: ranks of processes to which data will be send.
	 * \param numRecvProcs: Number of processes in pRecvProcMap. Also corresponds
	 *						to the size of pBufferSegSizes and to the number of
	 *						segments in pBuffer.
	 * \param tag: This tag will be used to identify send and received messages.
	 *			   Default: 1
	 */
		void distribute_data(void* pBufferOut, int* pBufferOutSegSizes,
							 int* pSenderProcMap, int numSenderProcs,
							 void* pBuffer, int* pBufferSegSizes,
							 int* pRecvProcMap, int numRecvProcs,
							 int tag = 1);
	private:
	///	holds an mpi-communicator.
	/**	A variable stores whether the communicator has to be freed when the
	 *	the wrapper is deleted.*/
		struct CommWrapper{
			CommWrapper(const MPI_Comm& comm,
						bool bReleaseComm);
			~CommWrapper();
			
			MPI_Comm			m_mpiComm;
			bool				m_bReleaseCommunicator;

		///	only contains data if m_mpiComm != MPI_COMM_WORLD
			std::vector<int>	m_procs;
		};
		
	///	Smart-pointer that encapsulates a CommWrapper.
		typedef SmartPtr<CommWrapper> SPCommWrapper;
		
	private:
	///	smart-pointer to an instance of a CommWrapper.
		SPCommWrapper	m_comm;
};

}//	end of namespace pcl


#endif
