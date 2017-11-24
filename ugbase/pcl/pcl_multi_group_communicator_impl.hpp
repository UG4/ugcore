/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__PCL_multi_group_communicator_impl
#define __H__PCL_multi_group_communicator_impl

#include "pcl_methods.h"
#include "pcl_datatype.h"
#include "pcl_reduce_traits.h"

namespace pcl {
	
template<typename T>
void MultiGroupCommunicator::
allreduce (const T *sendBuf, T *recvBuf,
		   size_t countPerGroup, ReduceOperation op) const
{
// UG_LOG("<dbg> allreduce STARTS\n");
//todo: If groups are big (e.g. >= 0.25*m_com.size()) an MPI allreduce
//		possibly on sub-communicators could be beneficial.

//todo:	Performance may be gained if multiple sends from one proc to another
//		would be handled in just one send.
	using namespace std;

	if(m_memberships.empty()){
		// UG_LOG("<dbg> allreduce ENDS EARLY\n");
		return;
	}

	const int tag = 4378;
	const int procRank = m_com.get_local_proc_id ();
	
	vector<T> receivedVals;
	size_t numSends = 0;
	size_t numRecvs = 0;

// UG_LOG("<dbg> 1\n");
//	(1) collect values of all procs on the root of each group (first proc in group)
//	get number of requried sends/receives
	{

		for(size_t imem = 0; imem < m_memberships.size(); ++imem){
			const size_t o = m_groupOffsets [imem];
			const size_t s = m_groupOffsets [imem + 1] - o;

			if (s == 1)	continue;

			if (m_groupMembers [0] == procRank)
				++ numRecvs;
			else
				++ numSends;
		}

		receivedVals.resize (numRecvs * countPerGroup);

	//	send / receive
		std::vector<MPI_Request> sendRequests(numSends);
		std::vector<MPI_Request> recvRequests(numRecvs);
		size_t sendCount = 0;
		size_t recvCount = 0;

		for(size_t imem = 0; imem < m_memberships.size(); ++imem){
			const size_t o = m_groupOffsets [imem];
			const size_t s = m_groupOffsets [imem + 1] - o;

			if (s == 1)	continue;

			if (m_groupMembers [0] == procRank){
				for(size_t i = 1; i < s; ++i){
					MPI_Irecv(&receivedVals.front() + recvCount * countPerGroup,
					          (int) countPerGroup * sizeof (T),
					          MPI_UNSIGNED_CHAR,	
							  m_com.get_proc_id ((size_t)m_groupMembers [o + i]),
							  tag,
							  m_com.get_mpi_communicator(),
							  &recvRequests[recvCount]);	
					++recvCount;
				}
			}
			else {
				//note: const_cast required for some MPI implementations...
				MPI_Isend(const_cast<T*>(sendBuf) + imem * countPerGroup,
				          (int) countPerGroup * sizeof (T),
				          MPI_UNSIGNED_CHAR,
				  		  m_com.get_proc_id ((size_t)m_groupMembers [0]),
				  		  tag,
				  		  m_com.get_mpi_communicator(),
				  		  &sendRequests[sendCount]);
				++sendCount;
			}
		}

		Waitall (sendRequests);
		Waitall (recvRequests);
	}
// UG_LOG("<dbg> 2\n");
//	(2) apply the reduce operation
	{
		Reducer<T> reducer (op);

		size_t recvHandled = 0;

		for(size_t imem = 0; imem < m_memberships.size(); ++imem){
			const size_t o = m_groupOffsets [imem];
			const size_t s = m_groupOffsets [imem + 1] - o;

			if (m_groupMembers [0] == procRank){
			//	copy data local data to recvBuf
				const size_t targetBaseInd = imem * countPerGroup;

				for(size_t i = 0; i < countPerGroup; ++i)
					recvBuf[targetBaseInd + i] = sendBuf[targetBaseInd + i];
				
			//	perform reduce operation for the data of each group member
				for(size_t iblock = 1; iblock < s; ++iblock){
					size_t srcBaseInd = recvHandled * countPerGroup;

					for(size_t i = 0; i < countPerGroup; ++i){
						reducer (recvBuf[targetBaseInd + i],
						         receivedVals[srcBaseInd + i]);
					}
					++recvHandled;
				}
			}
		}
	}
// UG_LOG("<dbg> 3\n");
//	(3) send values from group-roots to all group members
	{
	//	we have to send data to each process from which we received data and vice versa
		const size_t numSendBack = numRecvs;
		const size_t numRecvBack = numSends;

		std::vector<MPI_Request> sendRequests(numSendBack);
		std::vector<MPI_Request> recvRequests(numRecvBack);

		size_t sendCount = 0;
		size_t recvCount = 0;

		for(size_t imem = 0; imem < m_memberships.size(); ++imem){
			const size_t o = m_groupOffsets [imem];
			const size_t s = m_groupOffsets [imem + 1] - o;

			if (s == 1)	continue;

			if (m_groupMembers [0] == procRank){
				for(size_t i = 1; i < s; ++i){
					MPI_Isend(recvBuf + imem * countPerGroup,
				          (int) countPerGroup * sizeof (T),
				          MPI_UNSIGNED_CHAR,
				  		  m_com.get_proc_id ((size_t)m_groupMembers [o + i]),
				  		  tag,
				  		  m_com.get_mpi_communicator(),
				  		  &sendRequests[sendCount]);
					++sendCount;
				}
			}
			else {
				MPI_Irecv(recvBuf + imem * countPerGroup,
				          (int) countPerGroup * sizeof (T),
				          MPI_UNSIGNED_CHAR,
				  		  m_com.get_proc_id ((size_t)m_groupMembers [0]),
				  		  tag,
				  		  m_com.get_mpi_communicator(),
				  		  &recvRequests[recvCount]);
				++recvCount;
			}
		}

		Waitall (sendRequests);
		Waitall (recvRequests);
	}
// UG_LOG("<dbg> allreduce ENDS\n");
}

}//	end of namespace

#endif	//__H__PCL_multi_group_communicator_impl
