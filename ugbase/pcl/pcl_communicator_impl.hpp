//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

#ifndef __H__PCL__PCL_COMMUNICATOR_IMPL__
#define __H__PCL__PCL_COMMUNICATOR_IMPL__

#include <iostream>
#include <cassert>
#include "mpi.h"
#include "pcl_methods.h"
#include "pcl_communication_structs.h"
#include "pcl_communicator.h"
#include "common/profiler/profiler.h"

namespace pcl
{

template <class TLayout>
ParallelCommunicator<TLayout>::
ParallelCommunicator() :
	m_bSendBuffersFixed(true)
{
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
send_raw(int targetProc, void* pBuff, int bufferSize,
	     bool bSizeKnownAtTarget)
{
	assert(targetProc == -1 || targetProc >= 0 && targetProc < pcl::GetNumProcesses());

	ug::BinaryStream& stream = *m_streamPackOut.get_stream(targetProc);
	if(!bSizeKnownAtTarget)
		stream.write((char*)&bufferSize, sizeof(int));
		
	stream.write((char*)pBuff, bufferSize);
	m_bSendBuffersFixed = m_bSendBuffersFixed
						&& bSizeKnownAtTarget;
}
			   
////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
send_data(int targetProc, Interface& interface,
			  ICommunicationPolicy<TLayout>& commPol)
{
	assert(targetProc == -1 || targetProc >= 0 && targetProc < pcl::GetNumProcesses());

	ug::BinaryStream& stream = *m_streamPackOut.get_stream(targetProc);
	commPol.collect(stream, interface);
	m_bSendBuffersFixed = m_bSendBuffersFixed
						&& (commPol.get_required_buffer_size(interface) >= 0);
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
send_data(Layout& layout, ICommunicationPolicy<TLayout>& commPol)
{
//	through the the category_tag we're able to find the correct send method.
	send_data(layout, commPol, typename TLayout::category_tag());
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
send_data(Layout& layout,
		  ICommunicationPolicy<TLayout>& commPol,
		  const layout_tags::single_level_layout_tag&)
{
	typename Layout::iterator iter = layout.begin();
	typename Layout::iterator end = layout.end();
	
	commPol.begin_layout_collection(&layout);
	
	for(; iter != end; ++iter)
	{
		ug::BinaryStream& stream = *m_streamPackOut.get_stream(layout.proc_id(iter));
		commPol.collect(stream, layout.interface(iter));
		m_bSendBuffersFixed = m_bSendBuffersFixed
			&& (commPol.get_required_buffer_size(layout.interface(iter)) >= 0);
	}
	
	commPol.end_layout_collection();
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
send_data(Layout& layout,
		  ICommunicationPolicy<TLayout>& commPol,
		  const layout_tags::multi_level_layout_tag&)
{
	commPol.begin_layout_collection(&layout);
	
	for(size_t i = 0; i < layout.num_levels(); ++i)
	{
		typename Layout::iterator iter = layout.begin(i);
		typename Layout::iterator end = layout.end(i);
		
		for(; iter != end; ++iter)
		{
			ug::BinaryStream& stream = *m_streamPackOut.get_stream(layout.proc_id(iter));
			commPol.collect(stream, layout.interface(iter));
			m_bSendBuffersFixed = m_bSendBuffersFixed
				&& (commPol.get_required_buffer_size(layout.interface(iter)) >= 0);
		}
	}
	
	commPol.end_layout_collection();
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
receive_raw(int srcProc, ug::BinaryStream& binStreamOut,
			int bufferSize)
{
	m_extractorInfos.push_back(ExtractorInfo(srcProc, NULL,
											 NULL, NULL,
											 NULL, &binStreamOut,
											 bufferSize));
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
receive_raw(int srcProc, void* buffOut,
			int bufferSize)
{
	m_extractorInfos.push_back(ExtractorInfo(srcProc, NULL,
											 NULL, NULL,
											 buffOut, NULL,
											 bufferSize));
}
			
////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
receive_data(int srcProc, Interface& interface,
			ICommunicationPolicy<TLayout>& commPol)
{
	m_extractorInfos.push_back(ExtractorInfo(srcProc, &commPol
											 &interface, NULL,
											 NULL, NULL, 0));
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
receive_data(Layout& layout, ICommunicationPolicy<TLayout>& commPol)
{
	m_extractorInfos.push_back(ExtractorInfo(-1, &commPol,
											 NULL, &layout,
											 NULL, NULL, 0));
}

template <class TLayout>
template <class TLayoutMap>
void ParallelCommunicator<TLayout>::
exchange_data(TLayoutMap& layoutMap,
				const typename TLayoutMap::Key& keyFrom,
				const typename TLayoutMap::Key& keyTo,
				ICommunicationPolicy<TLayout>& commPol)
{
	if(layoutMap.template has_layout<Type>(keyFrom)){
		send_data(layoutMap.template get_layout<Type>(keyFrom), commPol);
	}
		
	if(layoutMap.template has_layout<Type>(keyTo)){
		receive_data(layoutMap.template get_layout<Type>(keyTo), commPol);
	}
}
							
////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
prepare_receiver_stream_pack(ug::StreamPack& streamPack,
							TLayout& layout)
{
	prepare_receiver_stream_pack(streamPack, layout,
								typename TLayout::category_tag());
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
prepare_receiver_stream_pack(ug::StreamPack& streamPack,
							TLayout& layout,
							const layout_tags::single_level_layout_tag&)
{
//	simply 'touch' the stream to make sure it's in the pack.
	for(typename TLayout::iterator li = layout.begin();
		li != layout.end(); ++li)
	{
		streamPack.get_stream(layout.proc_id(li));
	}
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
prepare_receiver_stream_pack(ug::StreamPack& streamPack,
							TLayout& layout,
							const layout_tags::multi_level_layout_tag&)
{
//	simply 'touch' the stream to make sure it's in the pack.
	for(size_t i = 0; i < layout.num_levels(); ++i)
	{
		for(typename TLayout::iterator li = layout.begin(i);
			li != layout.end(i); ++li)
		{
			streamPack.get_stream(layout.proc_id(li));
		}
	}
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
bool ParallelCommunicator<TLayout>::
collect_layout_buffer_sizes(TLayout& layout,
							ICommunicationPolicy<TLayout>& commPol,
							std::map<int, int>* pMapBuffSizesOut,
							const layout_tags::single_level_layout_tag&)
{
	for(typename TLayout::iterator li = layout.begin();
		li != layout.end(); ++li)
	{
	//	get the buffer size
		int buffSize = commPol.get_required_buffer_size(
								layout.interface(li));
		if(buffSize < 0){
		//	buffer sizes can't be determined
			return false;
		}
		else if(pMapBuffSizesOut){
		//	find the entry in the map
			std::map<int, int>::iterator iter = pMapBuffSizesOut->find(layout.proc_id(li));
			if(iter != pMapBuffSizesOut->end())
				iter->second += buffSize;
			else
				(*pMapBuffSizesOut)[layout.proc_id(li)] = buffSize;
		}
	}
	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
bool ParallelCommunicator<TLayout>::
collect_layout_buffer_sizes(TLayout& layout,
							ICommunicationPolicy<TLayout>& commPol,
							std::map<int, int>* pMapBuffSizesOut,
							const layout_tags::multi_level_layout_tag&)
{
//	iterate through all interfaces
	for(size_t i = 0; i < layout.num_levels(); ++i){
		for(typename TLayout::iterator li = layout.begin(i);
			li != layout.end(i); ++li)
		{
		//	get the buffer size
			int buffSize = commPol.get_required_buffer_size(
									layout.interface(li));
			if(buffSize < 0){
			//	buffer sizes can't be determined
				return false;
			}
			else if(pMapBuffSizesOut){
			//	find the entry in the map
				std::map<int, int>::iterator iter = pMapBuffSizesOut->find(layout.proc_id(li));
				if(iter != pMapBuffSizesOut->end())
					iter->second += buffSize;
				else
					(*pMapBuffSizesOut)[layout.proc_id(li)] = buffSize;
			}
		}
	}
	return true;
}
										
////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
extract_data(TLayout& layout, ug::StreamPack& streamPack, CommPol& extractor)
{
	extract_data(layout, streamPack,
				extractor,
				typename TLayout::category_tag());
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
extract_data(TLayout& layout, ug::StreamPack& streamPack, CommPol& extractor,
				const layout_tags::single_level_layout_tag&)
{
	extractor.begin_level_extraction(0);
//	extract data for the layouts interfaces
	for(typename Layout::iterator li = layout.begin();
		li != layout.end(); ++li)
	{
		extractor.extract(*streamPack.get_stream(layout.proc_id(li)),
						layout.interface(li));
	}
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void ParallelCommunicator<TLayout>::
extract_data(TLayout& layout, ug::StreamPack& streamPack, CommPol& extractor,
				const layout_tags::multi_level_layout_tag&)
{
//	extract data for the layouts interfaces
	for(size_t i = 0; i < layout.num_levels(); ++i)
	{
		extractor.begin_level_extraction(i);
		for(typename Layout::iterator li = layout.begin(i);
			li != layout.end(i); ++li)
		{
			extractor.extract(*streamPack.get_stream(layout.proc_id(li)),
							layout.interface(li));
		}
	}
}
				
////////////////////////////////////////////////////////////////////////
template <class TLayout>
bool ParallelCommunicator<TLayout>::
communicate()
{
//	PROFILE_FUNC();
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
			prepare_receiver_stream_pack(m_streamPackIn, *info.m_layout);
		}
	}

//	number of in and out-streams.
	size_t	numOutStreams = m_streamPackOut.num_streams();
	size_t	numInStreams = m_streamPackIn.num_streams();

//	used for mpi-communication.
	std::vector<MPI_Request> vSendRequests(numOutStreams);
	std::vector<MPI_Request> vReceiveRequests(numInStreams);
	
//	wait until data has been received
	std::vector<MPI_Status> vReceiveStates(numInStreams);//TODO: fix spelling!
	std::vector<MPI_Status> vSendStates(numOutStreams);//TODO: fix spelling!


////////////////////////////////////////////////
//	determine buffer sizes and communicate them if required.
	std::vector<int> vBufferSizesIn(numInStreams);
	bool allBufferSizesFixed = m_bSendBuffersFixed;
	
	if(allBufferSizesFixed)
	{
	//	a map with <procId, Size>. Will be used to collect stream-sizes
		std::map<int, int> mapBuffSizes;
	//	initialise all sizes with 0
		for(ug::StreamPack::iterator iter = m_streamPackIn.begin();
			iter != m_streamPackIn.end(); ++iter)
		{
			mapBuffSizes[iter->first] = 0;
		}
		
	//	iterate over all extractors and collect the buffer sizes
		for(typename ExtractorInfoList::iterator iter = m_extractorInfos.begin();
			iter != m_extractorInfos.end(); ++iter)
		{
			ExtractorInfo& info = *iter;
			if(info.m_srcProc > -1){
			//	the extractor only has a single interface.
				int buffSize = -1;
				if(info.m_extractor)
					buffSize = info.m_extractor->get_required_buffer_size(*info.m_interface);
				else
					buffSize = info.m_rawSize;

				if(buffSize < 0){
					allBufferSizesFixed = false;
					break;
				}
				else
					mapBuffSizes[info.m_srcProc] += buffSize;
			}
			else{
				if(!collect_layout_buffer_sizes(*info.m_layout,
												*info.m_extractor,
												&mapBuffSizes,
												typename TLayout::category_tag()))
				{
					allBufferSizesFixed = false;
					break;
				}
			}
		}
		
	//	if all buffer sizes are fixed, we'll copy them to vBufferSizes.
	//	to reduce the amount of possible errors, we'll iterate over
	//	m_streamPackIn...
		if(allBufferSizesFixed){
			int counter = 0;
			for(ug::StreamPack::iterator iter = m_streamPackIn.begin();
				iter != m_streamPackIn.end(); ++iter, ++counter)
			{
				vBufferSizesIn[counter] = mapBuffSizes[iter->first];
			}
		}
	}
	
	//	if the buffer size could not be determined, we have to communicate it.
	if(!allBufferSizesFixed)
	{
//		PROFILE_BEGIN(communicateBufferSizes);
		int sizeTag = 189345;//	an arbitrary number
		int counter;

	//	shedule receives first
		counter = 0;
		for(ug::StreamPack::iterator iter = m_streamPackIn.begin();
			iter != m_streamPackIn.end(); ++iter, ++counter)
		{
			MPI_Irecv(&vBufferSizesIn[counter], sizeof(int), MPI_UNSIGNED_CHAR,	
					iter->first, sizeTag, MPI_COMM_WORLD, &vReceiveRequests[counter]);
		}

	//	send buffer sizes
		counter = 0;
		for(ug::StreamPack::iterator iter = m_streamPackOut.begin();
			iter != m_streamPackOut.end(); ++iter, ++counter)
		{
			int streamSize = (int)iter->second->size();
			//int retVal =
			MPI_Isend(&streamSize, sizeof(int), MPI_UNSIGNED_CHAR,
					iter->first, sizeTag, MPI_COMM_WORLD, &vSendRequests[counter]);
		}

	//	TODO: this can be improved:
	//		instead of waiting for all, one could wait until one has finished and directly
	//		start copying the data to the local receive buffer. Afterwards on could continue
	//		by waiting for the next one etc...
		MPI_Waitall(numInStreams, &vReceiveRequests[0], &vReceiveStates[0]);
		MPI_Waitall(numOutStreams, &vSendRequests[0], &vSendStates[0]);
//		PROFILE_END();
	}

////////////////////////////////////////////////
//	communicate data.
//	PROFILE_BEGIN(communicateData);
	int dataTag = 749345;//	an arbitrary number

//	first shedule receives
	int counter = 0;
	for(ug::StreamPack::iterator iter = m_streamPackIn.begin();
		iter != m_streamPackIn.end(); ++iter, ++counter)
	{
	//	resize the buffer
		iter->second->resize(vBufferSizesIn[counter]);
		
	//	receive the data
		MPI_Irecv(iter->second->buffer(), vBufferSizesIn[counter], MPI_UNSIGNED_CHAR,	
				iter->first, dataTag, MPI_COMM_WORLD, &vReceiveRequests[counter]);
	}

//	now send data
	counter = 0;
	for(ug::StreamPack::iterator iter = m_streamPackOut.begin();
		iter != m_streamPackOut.end(); ++iter, ++counter)
	{
		//int retVal =
		MPI_Isend(iter->second->buffer(), iter->second->size(), MPI_UNSIGNED_CHAR,
				iter->first, dataTag, MPI_COMM_WORLD, &vSendRequests[counter]);
	}

//	TODO: this can be improved:
//		instead of waiting for all, one could wait until one has finished and directly
//		start copying the data to the local receive buffer. Afterwards on could continue
//		by waiting for the next one etc...
	MPI_Waitall(numInStreams, &vReceiveRequests[0], &vReceiveStates[0]);
	MPI_Waitall(numOutStreams, &vSendRequests[0], &vSendStates[0]);
//	PROFILE_END();
	
//	call the extractors with the received data
	for(typename ExtractorInfoList::iterator iter = m_extractorInfos.begin();
		iter != m_extractorInfos.end(); ++iter)
	{
		ExtractorInfo& info = *iter;
		if(info.m_srcProc > -1)
		{
		//	extract the data for single proc
			ug::BinaryStream& stream = *m_streamPackIn.get_stream(info.m_srcProc);
		//	this can be either an interface, a void* buffer or a
		//	binary-stream.
			if(info.m_interface){
				info.m_extractor->extract(stream, *info.m_interface);
			}
			else if(info.m_buffer){
				stream.read((char*)info.m_buffer, info.m_rawSize);
			}
			else{
				assert(info.m_stream && "ERROR in ParallelCommunicator::communicate: No valid receiver specified.");
				
				int rawSize = info.m_rawSize;
				if(rawSize < 0){
				//	the raw size is encoded in the buffer in this case.
					stream.read((char*)&rawSize, sizeof(int));
				}
				
				info.m_stream->reset();
				info.m_stream->resize(rawSize);
				stream.read((char*)info.m_stream->buffer(), rawSize);
			}
		}
		else
		{
		//	extract the data for a layout
		//	notify the extractor that extraction for a layout begins.
			info.m_extractor->begin_layout_extraction(info.m_layout);
			
		//	extract the data
			extract_data(*info.m_layout,
						m_streamPackIn,
						*info.m_extractor);
			
		//	notify the extractor that extraction is complete.
			info.m_extractor->end_layout_extraction();
		}
	}

//	clean up
	m_streamPackOut.clear();
	m_extractorInfos.clear();

//	reset m_bSendBuffersFixed
	m_bSendBuffersFixed = true;
	
//	done
	return true;
}

}//	end of namespace pcl

#endif
