//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

#ifndef __H__PCL__PCL_COMMUNICATION_IMPL__
#define __H__PCL__PCL_COMMUNICATION_IMPL__

#include <iostream>
#include "mpi.h"
#include "pcl_methods.h"
#include "pcl_base.h"
#include "pcl_communication.h"

namespace pcl
{
////////////////////////////////////////////////////////////////////////
template <class TLayout>
void Communicator<TLayout>::
send_data(int targetProc, Interface& interface,
			  ICommunicationPolicy<TLayout>& commPol)
{
	ug::BinaryStream& stream = *m_streamPackOut.get_stream(targetProc);
	commPol.collect(stream, interface);
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void Communicator<TLayout>::
send_data(Layout& layout, ICommunicationPolicy<TLayout>& commPol)
{
	typename Layout::iterator iter = layout.begin();
	typename Layout::iterator end = layout.end();

	for(; iter != end; ++iter)
	{
		ug::BinaryStream& stream = *m_streamPackOut.get_stream(layout.proc_id(iter));
		commPol.collect(stream, layout.interface(iter));
	}
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void Communicator<TLayout>::
receive_data(int srcProc, Interface& interface,
			ICommunicationPolicy<TLayout>& commPol)
{
	m_extractorInfos.push_back(ExtractorInfo(&commPol, srcProc, &interface, NULL));
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void Communicator<TLayout>::
receive_data(Layout& layout, ICommunicationPolicy<TLayout>& commPol)
{
	m_extractorInfos.push_back(ExtractorInfo(&commPol, -1, NULL, &layout));
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
bool Communicator<TLayout>::
communicate()
{
//	prepare receive streams
//TODO:	This should be done in a way so that the least possible amount
//		of data has to be reallocated (very often the map won't change
//		compared to the last communication step).
//	clear the map
	m_streamPackIn.clear();

//	iterate through all registered extractors and create entries for
//	the source-processes in the map (by simply 'touching' the entry).
	for(typename ExtractorInfoList::iterator iter = m_extractorInfos.begin();
		iter != m_extractorInfos.end(); ++iter)
	{
		ExtractorInfo& info = *iter;
		if(info.m_srcProc > -1)
			m_streamPackIn.get_stream(info.m_srcProc);
		else
		{
			for(typename Layout::iterator li = info.m_layout->begin();
				li != info.m_layout->end(); ++li)
			{
				m_streamPackIn.get_stream(info.m_layout->proc_id(li));
			}
		}
	}

//	number of in and out-streams.
	size_t	numOutStreams = m_streamPackOut.num_streams();
	size_t	numInStreams = m_streamPackIn.num_streams();

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
	for(ug::StreamPack::iterator iter = m_streamPackOut.begin();
		iter != m_streamPackOut.end(); ++iter, ++counter)
	{
		int streamSize = (int)iter->second->size();
		int retVal = MPI_Isend(&streamSize, sizeof(int), MPI_UNSIGNED_CHAR,
				iter->first, sizeTag, MPI_COMM_WORLD, &vSendRequests[counter]);
	}
	
//	now receive buffer sizes
	std::vector<int> vBufferSizesIn(numInStreams);
	counter = 0;
	for(ug::StreamPack::iterator iter = m_streamPackIn.begin();
		iter != m_streamPackIn.end(); ++iter, ++counter)
	{
		MPI_Irecv(&vBufferSizesIn[counter], sizeof(int), MPI_UNSIGNED_CHAR,	
				iter->first, sizeTag, MPI_COMM_WORLD, &vReceiveRequests[counter]);
	}
	
//	wait until data has been received
	std::vector<MPI_Status> vReceiveStatii(numInStreams);//TODO: fix spelling!
	
//	TODO: this can be improved:
//		instead of waiting for all, one could wait until one has finished and directly
//		start copying the data to the local receive buffer. Afterwards on could continue
//		by waiting for the next one etc...
	MPI_Waitall(numInStreams, &vReceiveRequests[0], &vReceiveStatii[0]);

////////////////////////////////////////////////
//	communicate data.
	int dataTag = 749345;//	an arbitrary number

//	first send data
	counter = 0;
	for(ug::StreamPack::iterator iter = m_streamPackOut.begin();
		iter != m_streamPackOut.end(); ++iter, ++counter)
	{
		int retVal = MPI_Isend(iter->second->buffer(), iter->second->size(), MPI_UNSIGNED_CHAR,
				iter->first, dataTag, MPI_COMM_WORLD, &vSendRequests[counter]);
	}

//	now receive data
	counter = 0;
	for(ug::StreamPack::iterator iter = m_streamPackIn.begin();
		iter != m_streamPackIn.end(); ++iter, ++counter)
	{
	//	resize the buffer
		iter->second->resize(vBufferSizesIn[counter]);
		
	//	receive the data
		MPI_Irecv(iter->second->buffer(), vBufferSizesIn[counter], MPI_UNSIGNED_CHAR,	
				iter->first, dataTag, MPI_COMM_WORLD, &vReceiveRequests[counter]);
	}

//	TODO: this can be improved:
//		instead of waiting for all, one could wait until one has finished and directly
//		start copying the data to the local receive buffer. Afterwards on could continue
//		by waiting for the next one etc...
	MPI_Waitall(numInStreams, &vReceiveRequests[0], &vReceiveStatii[0]);

//	call the extractors with the received data
	for(typename ExtractorInfoList::iterator iter = m_extractorInfos.begin();
		iter != m_extractorInfos.end(); ++iter)
	{
		ExtractorInfo& info = *iter;
		if(info.m_srcProc > -1)
		{
		//	extract the data for single proc
			info.m_extractor->extract(*m_streamPackIn.get_stream(info.m_srcProc),
										*info.m_interface);
		}
		else
		{
		//	extract the data for a layout
			Layout& layout = *info.m_layout;
		//	notify the extractor that extraction for a layout begins.
			info.m_extractor->begin_layout_extraction(info.m_layout);
		//	extract data for the layouts interfaces
			for(typename Layout::iterator li = layout.begin();
				li != layout.end(); ++li)
			{
				info.m_extractor->extract(*m_streamPackIn.get_stream(layout.proc_id(li)),
										layout.interface(li));
			}
			
		//	notify the extractor that extraction is complete.
			info.m_extractor->end_layout_extraction();
		}
	}

//	clean up
	m_streamPackOut.clear();
	m_extractorInfos.clear();
}

}//	end of namespace pcl

#endif
