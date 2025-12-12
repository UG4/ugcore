/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__PCL__pcl_process_communicator_impl__
#define __H__PCL__pcl_process_communicator_impl__

#include "pcl_process_communicator.h"

#include "common/util/vector_util.h"

namespace pcl {

template<typename TValue>
void ProcessCommunicator::
gatherv(std::vector<TValue>& recBufOut,
		 std::vector<TValue>& sendBuf, int root,
		 std::vector<int>* pSizesOut, std::vector<int>* pOffsetsOut) const
{
	if(is_local()) return;
	using namespace ug;

//todo: One could declare a special MPI_Datatype and could thus
//	directly work on pSizesOut and pOffsetsOut.
	std::vector sizes(this->size(), 0);
	std::vector offsets(this->size(), 0);

//	gather the sizes on root. Note that we send the actual data
//	in bytes later on.
	int localSize = static_cast<int>(sendBuf.size()) * sizeof(TValue);
	gather(&localSize, 1, PCL_DT_INT, GetDataPtr(sizes),
			1, PCL_DT_INT, root);

	int totalSize = 0;
	for(size_t i = 0; i < sizes.size(); ++i){
		offsets[i] = totalSize;
		totalSize += sizes[i];
	}

	//	root now knows all sizes. We can now gather the connections on root.
	recBufOut.resize(totalSize / sizeof(TValue));
	gatherv(GetDataPtr(sendBuf),
			 localSize,
			 PCL_DT_BYTE,
			 GetDataPtr(recBufOut),
			 GetDataPtr(sizes),
			 GetDataPtr(offsets),
			 PCL_DT_BYTE,
			 root);

//	send is complete now. If pSizesOut or pOffsetsOut was specified, we
//	fill them now.
//todo: if sizes and offsets would contain the actual number of entries
//		instead of bytes, this step could be skipped (see above).
	if(pSizesOut){
		std::vector<int>& sizesOut = *pSizesOut;
		sizesOut.resize(sizes.size());
		for(size_t i = 0; i < sizes.size(); ++i)
			sizesOut[i] = sizes[i] / sizeof(TValue);
	}

	if(pOffsetsOut){
		std::vector<int>& offsetsOut = *pOffsetsOut;
		offsetsOut.resize(offsets.size());
		for(size_t i = 0; i < offsets.size(); ++i)
			offsetsOut[i] = offsets[i] / sizeof(TValue);
	}
}

template<typename TValue>
void ProcessCommunicator::
allgatherv(std::vector<TValue>& recBufOut,
			std::vector<TValue>& sendBuf,
			std::vector<int>* pSizesOut,
			std::vector<int>* pOffsetsOut) const
{
	if(is_local()) {
		recBufOut.resize (sendBuf.size());
		for(size_t i = 0; i < sendBuf.size(); ++i){
			recBufOut[i] = sendBuf[i];
		}
		if(pSizesOut){
			pSizesOut->resize (1);
			pSizesOut->at(0) = sendBuf.size();
		}
		if(pOffsetsOut){
			pOffsetsOut->resize (1);
			pOffsetsOut->at(0) = 0;
		}
		return;
	}
	using namespace ug;

//todo: One could declare a special MPI_Datatype and could thus
//	directly work on pSizesOut and pOffsetsOut.
	std::vector<int> sizes(this->size(), 0);
	std::vector<int> offsets(this->size(), 0);

//	gather the sizes on root. Note that we send the actual data
//	in bytes later on.
	int localSize = static_cast<int>(sendBuf.size()) * sizeof(TValue);
	allgather(&localSize, 1, PCL_DT_INT, &sizes.front(),
			  1, PCL_DT_INT);

	int totalSize = 0;
	for(size_t i = 0; i < sizes.size(); ++i){
		offsets[i] = totalSize;
		totalSize += sizes[i];
	}

	if(totalSize == 0){
		recBufOut.resize(0);
		if(pSizesOut)
			pSizesOut->resize (0);
		if(pOffsetsOut)
			pOffsetsOut->resize (0);
		return;
	}

	//	all procs now know all sizes. We can now gather the connections.
	recBufOut.resize(totalSize / sizeof(TValue));
	allgatherv(GetDataPtr(sendBuf), localSize, PCL_DT_BYTE,
				GetDataPtr(recBufOut), GetDataPtr(sizes),
				GetDataPtr(offsets), PCL_DT_BYTE);

//	send is complete now. If pSizesOut or pOffsetsOut was specified, we
//	fill them now.
//todo: if sizes and offsets would contain the actual number of entries
//		instead of bytes, this step could be skipped (see above).
	if(pSizesOut){
		std::vector<int>& sizesOut = *pSizesOut;
		sizesOut.resize(sizes.size());
		for(size_t i = 0; i < sizes.size(); ++i)
			sizesOut[i] = sizes[i] / sizeof(TValue);
	}

	if(pOffsetsOut){
		std::vector<int>& offsetsOut = *pOffsetsOut;
		offsetsOut.resize(offsets.size());
		for(size_t i = 0; i < offsets.size(); ++i)
			offsetsOut[i] = offsets[i] / sizeof(TValue);
	}
}


template<typename T>
T ProcessCommunicator::
reduce(const T &t, ReduceOperation op, int rootProc) const
{
	T ret;
	reduce(&t, &ret, 1, DataTypeTraits<T>::get_data_type(), op, rootProc);
	return ret;
}

template<typename T>
void ProcessCommunicator::
reduce(const T *pSendBuff, T *pReceiveBuff, size_t count,
		  ReduceOperation op, int rootProc) const
{
	reduce(pSendBuff, pReceiveBuff, count, DataTypeTraits<T>::get_data_type(),
		   op, rootProc);
}

template<typename T>
void ProcessCommunicator::
reduce(const std::vector<T> &send, std::vector<T> &receive,
	   ReduceOperation op, int rootProc) const
{
	if(!send.empty()){
		receive.resize(send.size());
		reduce(&send[0], &receive[0], send.size(), op, rootProc);
	}
}


template<typename T>
T ProcessCommunicator::
allreduce(const T &t, ReduceOperation op) const
{
	T ret;
	allreduce(&t, &ret, 1, DataTypeTraits<T>::get_data_type(), op);
	return ret;
}


template<typename T>
void ProcessCommunicator::
allreduce(const T *pSendBuff, T *pReceiveBuff, size_t count, ReduceOperation op) const
{
	allreduce(pSendBuff, pReceiveBuff, count, DataTypeTraits<T>::get_data_type(), op);
}

template<typename T>
void ProcessCommunicator::
allreduce(const std::vector<T> &send, std::vector<T> &receive,
		  ReduceOperation op) const
{
	if(!send.empty()){
		receive.resize(send.size());
		allreduce(&send[0], &receive[0], send.size(), op);
	}
}



template<typename T>
void ProcessCommunicator::
broadcast(T &t, int root) const
{
	broadcast(t, root, typename DataTypeTraits<T>::supported());
}


template<typename T>
void ProcessCommunicator::
broadcast(T &t, int root, DataTypeDirectlySupported d) const
{
	broadcast(&t, 1, DataTypeTraits<T>::get_data_type(), root);
}

template<typename T>
void ProcessCommunicator::
broadcast(T &t, int root, DataTypeIndirectlySupported d) const
{
	ug::BinaryBuffer buf;
	if(ProcRank() == root)
	{
		Serialize(buf, t);
		broadcast(buf, root);
	}
	else
	{
		broadcast(buf, root);
		Deserialize(buf, t);
	}
}


template<typename T>
void ProcessCommunicator::
broadcast(T *p, size_t size, int root) const
{
	broadcast(p, size, DataTypeTraits<T>::get_data_type(), root);
}

}//	end of namespace

#endif
