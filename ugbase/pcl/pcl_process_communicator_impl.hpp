// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 03.02.2011 (m,d,y)

#ifndef __H__PCL__pcl_process_communicator_impl__
#define __H__PCL__pcl_process_communicator_impl__

#include "common/util/vector_util.h"

namespace pcl
{

template<class TValue>
void ProcessCommunicator::
gatherv(std::vector<TValue>& recBufOut,
		 std::vector<TValue>& sendBuf, int root,
		 std::vector<int>* pSizesOut, std::vector<int>* pOffsetsOut) const
{
	if(is_local()) return;
	using namespace ug;

//todo: One could declare a special MPI_Datatype and could thus
//	directly work on pSizesOut and pOffsetsOut.
	std::vector<int> sizes(this->size(), 0);
	std::vector<int> offsets(this->size(), 0);

//	gather the sizes on root. Note that we send the actual data
//	in bytes later on.
	int localSize = (int)sendBuf.size() * sizeof(TValue);
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

template<class TValue>
void ProcessCommunicator::
allgatherv(std::vector<TValue>& recBufOut,
			std::vector<TValue>& sendBuf,
			std::vector<int>* pSizesOut,
			std::vector<int>* pOffsetsOut) const
{
	if(is_local()) return;
	using namespace ug;

//todo: One could declare a special MPI_Datatype and could thus
//	directly work on pSizesOut and pOffsetsOut.
	std::vector<int> sizes(this->size(), 0);
	std::vector<int> offsets(this->size(), 0);

//	gather the sizes on root. Note that we send the actual data
//	in bytes later on.
	int localSize = (int)sendBuf.size() * sizeof(TValue);
	allgather(&localSize, 1, PCL_DT_INT, &sizes.front(),
			  1, PCL_DT_INT);

	int totalSize = 0;
	for(size_t i = 0; i < sizes.size(); ++i){
		offsets[i] = totalSize;
		totalSize += sizes[i];
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
reduce(const T &t, pcl::ReduceOperation op, int rootProc) const
{
	T ret;
	reduce(&t, &ret, 1, DataTypeTraits<T>::get_data_type(), op, rootProc);
	return ret;
}

template<typename T>
void ProcessCommunicator::
reduce(const T *pSendBuff, T *pReceiveBuff, size_t count,
		  pcl::ReduceOperation op, int rootProc) const
{
	reduce(pSendBuff, pReceiveBuff, count, DataTypeTraits<T>::get_data_type(),
		   op, rootProc);
}

template<typename T>
void ProcessCommunicator::
reduce(const std::vector<T> &send, std::vector<T> &receive,
	   pcl::ReduceOperation op, int rootProc) const
{
	if(send.size() > 0){
		receive.resize(send.size());
		reduce(&send[0], &receive[0], send.size(), op, rootProc);
	}
}


template<typename T>
T ProcessCommunicator::
allreduce(const T &t, pcl::ReduceOperation op) const
{
	T ret;
	allreduce(&t, &ret, 1, DataTypeTraits<T>::get_data_type(), op);
	return ret;
}


template<typename T>
void ProcessCommunicator::
allreduce(const T *pSendBuff, T *pReceiveBuff, size_t count, pcl::ReduceOperation op) const
{
	allreduce(pSendBuff, pReceiveBuff, count, DataTypeTraits<T>::get_data_type(), op);
}

template<typename T>
void ProcessCommunicator::
allreduce(const std::vector<T> &send, std::vector<T> &receive,
		  pcl::ReduceOperation op) const
{
	if(send.size() > 0){
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
	if(pcl::ProcRank() == root)
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
