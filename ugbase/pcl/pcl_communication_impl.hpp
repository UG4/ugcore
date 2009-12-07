//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

#ifndef __H__PCL__PCL_COMMUNICATION_IMPL__
#define __H__PCL__PCL_COMMUNICATION_IMPL__

#include <iostream>
#include "mpi.h"
#include "pcl_base.h"
#include "pcl_communication.h"

namespace pcl
{
////////////////////////////////////////////////////////////////////////
template <class TElementGroup>
void Communicator<TElementGroup>::
collect_data(int targetProc, Interface& interface,
			  ICollector<TElementGroup>& collector)
{
	ug::BinaryStream& stream = m_streamMapOut[targetProc];
	collector.collect(stream, interface);
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup>
void Communicator<TElementGroup>::
collect_data(Layout& layout, ICollector<TElementGroup>& collector)
{
	typename Layout::iterator iter = layout.begin();
	typename Layout::iterator end = layout.end();

	for(; iter != end; ++iter)
	{
		ug::BinaryStream& stream = m_streamMapOut[layout.proc_id(iter)];
		collector.collect(stream, layout.interface(iter));
	}
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup>
void Communicator<TElementGroup>::
await_data(int srcProc, Interface& interface,
			IExtractor<TElementGroup>& extractor)
{
	m_extractorInfos.push_back(ExtractorInfo(&extractor, srcProc, &interface, NULL));
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup>
void Communicator<TElementGroup>::
await_data(Layout& layout, IExtractor<TElementGroup>& extractor)
{
	m_extractorInfos.push_back(ExtractorInfo(&extractor, -1, NULL, &layout));
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup>
bool Communicator<TElementGroup>::
communicate()
{
//	prepare receive streams
//TODO:	This should be done in a way so that the least possible amount
//		of data has to be reallocated (very often the map won't change
//		compared to the last communication step).
//	clear the map
	m_streamMapIn = StreamMap();

//	iterate through all registered extractors and create entries for
//	the source-processes in the map (by simply 'touching' the entry).
	for(typename ExtractorInfoList::iterator iter = m_extractorInfos.begin();
		iter != m_extractorInfos.end(); ++iter)
	{
		ExtractorInfo& info = *iter;
		if(info.m_srcProc > -1)
			m_streamMapIn[info.m_srcProc];
		else
		{
			for(typename Layout::iterator li = info.m_layout->begin();
				li != info.m_layout->end(); ++li)
				m_streamMapIn[info.m_layout->proc_id(li)];
		}
	}

//	number of in and out-streams.
	size_t	numOutStreams = m_streamMapOut.size();
	size_t	numInStreams = m_streamMapIn.size();

//	used for mpi-communication.
	std::vector<MPI_Request> vSendRequests(numOutStreams);
	std::vector<MPI_Request> vReceiveRequests(numInStreams);

////////////////////////////////////////////////
//	communicate buffer sizes.
//TODO:	this could be avoided if extractors would calculate the
//		size of the data that they receive.

	int sizeTag = 189345;//	an arbitrary number
	int counter;

//	first send buffer sizes
	counter = 0;
	for(StreamMap::iterator iter = m_streamMapOut.begin();
		iter != m_streamMapOut.end(); ++iter, ++counter)
	{
		int streamSize = (int)iter->second.size();
		int retVal = MPI_Isend(&streamSize, sizeof(int), MPI_UNSIGNED_CHAR,
				iter->first, sizeTag, MPI_COMM_WORLD, &vSendRequests[counter]);
	}
	
//	now receive buffer sizes
	std::vector<int> vBufferSizes(numInStreams);
	counter = 0;
	for(StreamMap::iterator iter = m_streamMapIn.begin();
		iter != m_streamMapIn.end(); ++iter, ++counter)
	{
		MPI_Irecv(&vBufferSizes[counter], sizeof(int), MPI_UNSIGNED_CHAR,	
				iter->first, sizeTag, MPI_COMM_WORLD, &vReceiveRequests[counter]);
	}
	
//	wait until data has been received
	std::vector<MPI_Status> vReceiveStatii(numInStreams);//TODO: fix spelling!
	
//	TODO: this can be improved:
//		instead of waiting for all, one could wait until one has finished and directly
//		start copying the data to the local receive buffer. Afterwards on could continue
//		by waiting for the next one etc...
	MPI_Waitall(numInStreams, &vReceiveRequests[0], &vReceiveStatii[0]);

//	resize receive streams
	counter = 0;
	for(StreamMap::iterator iter = m_streamMapIn.begin();
		iter != m_streamMapIn.end(); ++iter, ++counter)
	{
		iter->second.resize(vBufferSizes[counter]);
	}

////////////////////////////////////////////////
//	communicate data.
	int dataTag = 749345;//	an arbitrary number

//	first send data
	counter = 0;
	for(StreamMap::iterator iter = m_streamMapOut.begin();
		iter != m_streamMapOut.end(); ++iter, ++counter)
	{
		int retVal = MPI_Isend(iter->second.buffer(), vBufferSizes[counter], MPI_UNSIGNED_CHAR,
				iter->first, dataTag, MPI_COMM_WORLD, &vSendRequests[counter]);
	}
	
//	now receive data
	counter = 0;
	for(StreamMap::iterator iter = m_streamMapIn.begin();
		iter != m_streamMapIn.end(); ++iter, ++counter)
	{
		MPI_Irecv(iter->second.buffer(), vBufferSizes[counter], MPI_UNSIGNED_CHAR,	
				iter->first, dataTag, MPI_COMM_WORLD, &vReceiveRequests[counter]);
	}

//	TODO: this can be improved:
//		instead of waiting for all, one could wait until one has finished and directly
//		start copying the data to the local receive buffer. Afterwards on could continue
//		by waiting for the next one etc...
	MPI_Waitall(numInStreams, &vReceiveRequests[0], &vReceiveStatii[0]);
}

}//	end of namespace pcl

#endif
