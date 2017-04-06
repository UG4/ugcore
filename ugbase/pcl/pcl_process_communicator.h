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

#ifndef __H__PCL__PCL_PROCESS_COMMUNICATOR__
#define __H__PCL__PCL_PROCESS_COMMUNICATOR__

#include <iostream>
#include <map>
#include <vector>
#include "pcl_methods.h"
#include "common/util/smart_pointer.h"
#include "common/util/binary_stream.h"
#include "common/util/binary_buffer.h"
#include "common/error.h"

namespace pcl
{

/// \addtogroup pcl
/// \{

///	values that can be passed to a ProcessCommunicators constructor.
enum ProcessCommunicatorDefaults
{
	PCD_EMPTY = 0, //< an empty communicator
	PCD_WORLD = 1, //< a communicator to all available processes
	PCD_LOCAL = 2  //< a local communicator, simulating current proc is the only proc
};


/** A ProcessCommunicator is a very lightweight object that can be passed
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
		inline bool empty() const		{return !is_local() && m_comm->m_mpiComm == MPI_COMM_NULL;}

	/// return true if the communicator is PCD_WORLD
		inline bool is_world() const	{ return !is_local() && m_comm->m_mpiComm == PCL_COMM_WORLD; }
		
	/// return true if the communicator is local, simulating current proc is the only proc
		inline bool is_local() const {return m_comm.valid() == false;}

	///	returns the size of the communicator
		size_t size() const;
		
	///	returns the i-th process in the communicator
		int get_proc_id(size_t index) const;

	/// returns true if we are the i-th process in the communicator
		bool is_proc_id(size_t index) const
		{
			return get_proc_id(index) == ProcRank();
		}

	///	returns the proc-id relative to this communicator
	/**	This method has a worst time complexity of O(n)*/
		int get_local_proc_id(int globalProcID = pcl::ProcRank()) const;

	///	returns the mpi-communicator, in case someone needs it
		MPI_Comm get_mpi_communicator()	{return m_comm->m_mpiComm;}

	///	creates a new communicator containing a subset of the current communicator
	/**	Note that this method has to be called by all processes in the current
	 *	communicator - even if they don't want to participate in the new one.*/
		ProcessCommunicator create_sub_communicator(bool participate) const;
		
	/**	Make sure that all processes which participate in the current communicator
	 * call this method with the same parameters! Process indices are to be specified
	 * relative to the current communicator.*/
		ProcessCommunicator create_sub_communicator(std::vector<int> &newProcs) const;

	/**	Make sure that all processes call this method with the same parameters!
	 * \{ */
		static ProcessCommunicator create_communicator(std::vector<int> &newGlobalProcs);
		static ProcessCommunicator create_communicator(size_t first, size_t num);
	/** \} */
		

	///	performs MPI_Gather on the processes of the communicator.
	/**	This method synchronises involved processes.
	 * \param sendBuf starting address of send buffer (choice)
	 * \param sendCount number of elements in send buffer (integer)
	 * \param sendType data type of send buffer elements (handle)
	 * \param recBuf only significant for root process. All gathered data is written here.
	 * \param recCount number of elements received from each process (integer)
	 * \param recType data type of receive buffer elements (handle)
	 * \param root The rank of the process that receives all the data.*/
		void gather(const void* sendBuf, int sendCount, DataType sendType,
					void* recBuf, int recCount, DataType recType, int root) const;

	///	performs MPI_Scatter on the processes of the communicator
	/** \param sendBuf starting address of send buffer (choice)
	 * \param sendCount number of elements in send buffer (integer)
	 * \param sendType data type of send buffer elements (handle)
	 * \param recBuf only significant for root process. All gathered data is written here.
	 * \param recCount number of elements received from each process (integer)
	 * \param recType data type of receive buffer elements (handle)
	 * \param root The rank of the process that receives all the data.*/
		void scatter(const void* sendBuf, int sendCount, DataType sendType,
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
	 * \param root The rank of the process that receives all the data.*/
		void gatherv(const void* sendBuf, int sendCount, DataType sendType,
						void* recBuf, int* recCounts, int* displs,
						DataType recType, int root) const;

	///	gathers variable arrays on one process.
	/**	The arrays specified in sendBuf will be copied to recBufOut
	 * on root. recBufOut is thus only important on root. The order
	 * of the arrays is the same as the order of the sending processes
	 * in this ProcessCommnicator.
	 * You can optionally specify an array pSizesOut, which will be
	 * filled with the size of the array that was received from each
	 * process (only relevant on root) and an array in which the
	 * offsets of each array are stored - that means the first index
	 * of each array - again in the order of processes in this ProcessCommunicator
	 * (only relevant on root, too).
	 *
	 * Note that recBufOut, pSizesOut and pOffsetsOut will be resized as required.
	 *
	 * Please note that this method communicates twice.*/
		template<class TValue>
		void gatherv(std::vector<TValue>& recBufOut,
					 std::vector<TValue>& sendBuf, int root,
					 std::vector<int>* pSizesOut = NULL,
					 std::vector<int>* pOffsetsOut = NULL) const;

	///	performs MPI_Allgather on the processes of the communicator.
	/**	This method synchronises involved processes.
	 * \param sendBuf starting address of send buffer (choice)
	 * \param sendCount number of elements in send buffer (integer)
	 * \param sendType data type of send buffer elements (handle)
	 * \param recBuf starting address of receive buffer (choice)
	 * \param recCount number of elements received from any process (integer)
	 * \param recType data type of receive buffer elements (handle)*/
		void allgather(const void* sendBuf, int sendCount, DataType sendType,
					   void* recBuf, int recCount, DataType recType) const;

	///	performs MPI_Allgatherv on the processes of the communicator.
	/**	This method synchronises involved processes.
	 * \param sendBuf 	starting address of send buffer (choice)
	 * \param sendCount 	number of elements in send buffer (integer)
	 * \param sendType 	data type of send buffer elements (handle)
	 * \param recCounts 	integer array (of length group size) containing the number of elements that are received from each process
	 * \param recBuf starting address of receive buffer (choice)
	 * \param displs 	integer array (of length group size). Entry i specifies the displacement (relative to recvbuf ) at which to place the incoming data from process i
	 * \param recType 	data type of receive buffer elements (handle)*/
		void allgatherv(const void* sendBuf, int sendCount, DataType sendType,
						void* recBuf, int* recCounts, int* displs,
						DataType recType) const;

		void alltoall(const void* sendBuf, int sendCount, DataType sendType,
		    		  void* recBuf, int recCount, DataType recType);

	///	gathers variable arrays on all processes.
	/**	The arrays specified in sendBuf will be copied to all processes
	 * in the ProcessCommunicator. The order of the arrays in recBufOut
	 * is the same as the order of the sending processes in this
	 * ProcessCommnicator.
	 * You can optionally specify an array pSizesOut, which will be
	 * filled with the size of the array that was received from each
	 * process and an array in which the offsets of each array are stored -
	 * that means the first index of each array - again in the order of
	 * processes in this ProcessCommunicator.
	 *
	 * Note that recBufOut, pSizesOut and pOffsetsOut will be resized as required.
	 *
	 * Please note that this method communicates twice.*/
		template<class TValue>
		void allgatherv(std::vector<TValue>& recBufOut,
						std::vector<TValue>& sendBuf,
						std::vector<int>* pSizesOut = NULL,
						std::vector<int>* pOffsetsOut = NULL) const;

		
	///	performs MPI_Reduce on the processes of the communicator.
	/**	This method synchronises involved processes.*/
		void reduce(const void* sendBuf, void* recBuf, int count,
					DataType type, ReduceOperation op, int rootProc) const;

	/** simplified reduce for size=1. calls reduce for parameter t, and then returns the result
	 * compiler error for unsupported datatypes
	 * \param t the input parameter
	 * \param op the Reduce Operation
	 * \param rootProc	the process onto which the result will be reduced.
	 * \return the reduced result*/
		template<typename T>
		T reduce(const T &t, pcl::ReduceOperation op, int rootProc) const;

	/** simplified reduce for buffers.
	 * \param pSendBuff the input buffer
	 * \param pReceiveBuff the output buffer
	 * \param count number of elements in the input/output buffers
	 * \param op the Reduce Operation
	 * \param rootProc	the process onto which the result will be reduced.*/
		template<typename T>
		void reduce(const T *pSendBuff, T *pReceiveBuff, size_t count,
					pcl::ReduceOperation op, int rootProc) const;

	/** simplified reduce for std::vector. The receive-buffer is automatically
	 * resized as required.*/
		template<typename T>
		void reduce(const std::vector<T> &send, std::vector<T> &receive,
					   pcl::ReduceOperation op, int rootProc) const;

	///	overload for size_t
		size_t reduce(const size_t &t, pcl::ReduceOperation op, int rootProc) const;

	///	performs MPI_Allreduce on the processes of the communicator.
	/**	This method synchronises involved processes.*/
		void allreduce(const void* sendBuf, void* recBuf, int count,
					   DataType type, ReduceOperation op) const;

	/** simplified allreduce for size=1. calls allreduce for parameter t,
	 * and then returns the result.
	 * \param t the input parameter
	 * \param op the Reduce Operation
	 * \return the reduced result*/
		template<typename T>
		T allreduce(const T &t, pcl::ReduceOperation op) const;

	///	overload for size_t
		size_t allreduce(const size_t &t, pcl::ReduceOperation op) const;

	/** simplified allreduce for buffers.
	 * \param pSendBuff the input buffer
	 * \param pReceiveBuff the output buffer
	 * \param count number of elements in the input/output buffers
	 * \param op the Reduce Operation*/
		template<typename T>
		void allreduce(const T *pSendBuff, T *pReceiveBuff, size_t count,
					   pcl::ReduceOperation op) const;

	/** simplified allreduce for std::vector. The receive-buffer is automatically
	 * resized as required.*/
		template<typename T>
		void allreduce(const std::vector<T> &send, std::vector<T> &receive,
					   pcl::ReduceOperation op) const;


	/** performs a MPI_Bcast
	 * @param v		pointer to data
	 * @param size	size of data
	 * @param type	type of data
	 * @param root	root process, that distributes data*/
		void broadcast(void *v, size_t size, DataType type, int root=0) const;

	/** simplified broadcast for supported datatypes
	 * compiler error for unsupported datatypes
	 * you can cast to unsigned char to broadcast arbitrary fixed data
	 * @param p		pointer to data
	 * @param size	number of T elements the pointer p is pointing to. default 1
	 * @param root	process that distributes data (default 0)*/
		template<typename T>
		inline void broadcast(T *p, size_t size=1, int root=0) const;

	/** Bcast for objects
	 * @param v		object to be broadcasted (uses Serialize/Deserialize)
	 * @param root	process that distributes data (default 0)
	 * @sa Serialize, Deserialize*/
		template<typename T>
		inline void broadcast(T &t, int root=0) const;

	/// broadcast function for directly supported types
		template<typename T>
		inline void broadcast(T &t, int root, DataTypeDirectlySupported d) const;
	/// broadcast function for indirectly supported types (using Serialize/Deserialize)
		template<typename T>
		void broadcast(T &t, int root, DataTypeIndirectlySupported d) const;

	/// overload for size_t
		void broadcast(size_t &s, int root=0) const;

	/** broadcast of BinaryBuffers
	 * @param buf	Binary buffer in/out
	 * @param root	root processor*/
		void broadcast(ug::BinaryBuffer &buf, int root=0) const;


	///	this method will not return until all processes in the communicator have called it.
		void barrier() const;


	///	sends data with the given tag to the specified process.
	/**	This method waits until the data has been sent.*/
		void send_data(void* pBuffer, int bufferSize, int destProc, int tag) const;

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
	 *	\param tag: A tag that tags the message. Use the same tag in receive_data.*/
		void send_data(void* pBuffer, int* pBufferSegSizes,
					   int* pRecProcMap, int numRecProcs, int tag) const;

	///	receives data from srcPrc with the specified tag.
	/**	This method waits until the data has been received*/
		void receive_data(void* pBuffOut, int bufferSize, int srcProc, int tag) const;

	///	sends and receives data to / from multiple processes.
	/** \param pBufferOut: Received data is written to this buffer.
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
	 *			   Default: 1*/
		void distribute_data(void* pBufferOut, int* pBufferOutSegSizes,
							 int* pSenderProcMap, int numSenderProcs,
							 void* pBuffer, int* pBufferSegSizes,
							 int* pRecvProcMap, int numRecvProcs,
							 int tag = 1) const;

	///	sends and receives data to/from multiple processes
	/**	This method automatically determines the size of the distributed
	 * data and writes it to a binary stream.
	 * Note that it has to communicate twice, since the buffer-sizes also
	 * have to be communicated.
	 *
	 * \param recvBufOut	Received data will be written to this buffer.
	 * \param segSizesOut	Array to which the block-sizes received from
	 * 						each process will be written.
	 * 						Has to have size numRecvFroms.
	 * \param recvFromRanks	Array containing the ranks which send data to
	 * 						this process. Has to have size numRecvFroms.
	 * \param numRecvFroms	Specifies from how many processes this process
	 * 						will receive data.
	 * \param sendBuf		Contains the data which will be send to other
	 * 						processes. Make sure that it is big enough
	 * 						(sum of all sendSegSizes).
	 * \param sendSegSizes	The i-th entry corresponds to the block-size
	 * 						which will be send to the i-th process in
	 * 						sendToRanks. Has to have size numSendTos.
	 * \param sendToRanks	An array of process ids, which defines to where
	 * 						data shall be sent. Has to have size numSendTos.
	 * \param numSendTos	Specifies to how many processes data will be sent.*/
		void distribute_data(ug::BinaryBuffer& recvBufOut, int* segSizesOut,
							int* recvFromRanks, int numRecvFroms,
							void* sendBuf, int* sendSegSizes,
							int* sendToRanks, int numSendTos,
							 int tag = 1) const;

	///	sends and receives data to/from multiple processes
	/** Note that it has to communicate twice, since the buffer-sizes also
	 * have to be communicated.
	 *
	 * \param recvBufs		Received data will be written to this buffers.
	 *						This array has to be of the size numRecvs
	 * \param recvFromRanks	Array containing the ranks from which data
	 * 						shall be received. Has to have size numRecvs.
	 * \param numRecvs		Specifies from how many processes this process
	 * 						will receive data.
	 * \param sendBufs		Array of buffers whose data will be send to other
	 * 						processes. Has to be of size numSends
	 * \param sendToRanks	An array of process ids, which defines to where
	 * 						data shall be sent. Has to have size numSends.
	 * \param numSendTos	Specifies to how many processes data will be sent.*/
		void distribute_data(ug::BinaryBuffer* recvBufs, int* recvFromRanks, int numRecvs,
							 ug::BinaryBuffer* sendBufs, int* sendToRanks, int numSendTos,
							 int tag = 1) const;
	private:
	///	holds an mpi-communicator.
	/**	A variable stores whether the communicator has to be freed when the
	 *	the wrapper is deleted.*/
		struct CommWrapper{
		///	initializes the commWrapper with PCL_COMM_WORLD
			CommWrapper();
			CommWrapper(const MPI_Comm& comm,
						bool bReleaseComm);
			~CommWrapper();
			
			MPI_Comm			m_mpiComm;
			bool				m_bReleaseCommunicator;

		///	only contains data if m_mpiComm != PCL_COMM_WORLD
			std::vector<int>	m_procs;
		};
		
	///	Smart-pointer that encapsulates a CommWrapper.
		typedef SmartPtr<CommWrapper> SPCommWrapper;
		
	private:
	///	smart-pointer to an instance of a CommWrapper.
		SPCommWrapper	m_comm;

};

std::ostream &operator << (std::ostream &out, const ProcessCommunicator &processCommunicator);

// end group pcl
/// \}

}//	end of namespace pcl

////////////////////////////////
//	include implementation
#include "pcl_process_communicator_impl.hpp"

#endif
