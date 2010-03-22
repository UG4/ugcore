//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d18

#ifndef __H__PCL__SCL_COMMUNICATOR_IMPL__
#define __H__PCL__SCL_COMMUNICATOR_IMPL__

#include <iostream>
#include "cl_base.h"
#include "scl_communicator.h"

namespace pcl
{
////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerialCommunicator<TLayout>::
send_data(int targetProc, Interface& interface,
			  ICommunicationPolicy<TLayout>& commPol)
{
	ug::BinaryStream& stream = *m_streamPackMap[interface.get_local_src_id()].
									get_stream(targetProc);
	commPol.collect(stream, interface);
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerialCommunicator<TLayout>::
send_data(Layout& layout, ICommunicationPolicy<TLayout>& commPol)
{
//	through the the category_tag we're able to find the correct send method.
	send_data(layout, commPol, typename TLayout::category_tag());
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerialCommunicator<TLayout>::
send_data(Layout& layout,
		  ICommunicationPolicy<TLayout>& commPol,
		  const layout_tags::single_level_layout_tag&)
{
	typename Layout::iterator iter = layout.begin();
	typename Layout::iterator end = layout.end();

	for(; iter != end; ++iter)
	{
		ug::BinaryStream& stream = *m_streamPackMap[layout.get_local_src_id()].
											get_stream(layout.proc_id(iter));
		commPol.collect(stream, layout.interface(iter));
	}
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerialCommunicator<TLayout>::
send_data(Layout& layout,
		  ICommunicationPolicy<TLayout>& commPol,
		  const layout_tags::multi_level_layout_tag&)
{
	for(size_t i = 0; i < layout.num_levels(); ++i)
	{
		typename Layout::iterator iter = layout.begin(i);
		typename Layout::iterator end = layout.end(i);
		
//todo: this behaviour is not well documented. it is not even clear if it is a good behaviour!
		int localSrcID = layout.get_local_src_id();
		if(localSrcID == -1)
			localSrcID = i;
			
		for(; iter != end; ++iter)
		{
			ug::BinaryStream& stream = *m_streamPackMap[localSrcID].
											get_stream(layout.proc_id(iter));
			commPol.collect(stream, layout.interface(iter));
		}
	}
}
		  
////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerialCommunicator<TLayout>::
receive_data(int srcProc, Interface& interface,
			ICommunicationPolicy<TLayout>& commPol)
{
	m_extractorInfos.push_back(ExtractorInfo(&commPol, srcProc, &interface, NULL));
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerialCommunicator<TLayout>::
receive_data(Layout& layout, ICommunicationPolicy<TLayout>& commPol)
{
	m_extractorInfos.push_back(ExtractorInfo(&commPol, -1, NULL, &layout));
}

template <class TLayout>
template <class TLayoutMap>
void SerialCommunicator<TLayout>::
exchange_data(TLayoutMap& layoutMap,
				typename TLayoutMap::Key keyFrom,
				typename TLayoutMap::Key keyTo,
				ICommunicationPolicy<TLayout>& commPol)
{
	if(layoutMap.template has_layout<Type>(keyFrom))
		send_data(layoutMap.template get_layout<Type>(keyFrom), commPol);
		
	if(layoutMap.template has_layout<Type>(keyTo))
		receive_data(layoutMap.template get_layout<Type>(keyTo), commPol);
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerialCommunicator<TLayout>::
extract_data(TLayout& layout, ug::StreamPack& streamPack, CommPol& extractor)
{
	extract_data(layout, streamPack,
				extractor,
				typename TLayout::category_tag());
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerialCommunicator<TLayout>::
extract_data(TLayout& layout, ug::StreamPack& streamPack, CommPol& extractor,
				const layout_tags::single_level_layout_tag&)
{
	extractor.begin_level_extraction(0);
//	extract data for the layouts interfaces
	for(typename Layout::iterator li = layout.begin();
		li != layout.end(); ++li)
	{
		extractor.extract(*streamPack.get_stream(layout.get_local_src_id()),
						layout.interface(li));
	}
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerialCommunicator<TLayout>::
extract_data(TLayout& layout, ug::StreamPack& streamPack, CommPol& extractor,
				const layout_tags::multi_level_layout_tag&)
{
//	extract data for the layouts interfaces
	for(size_t i = 0; i < layout.num_levels(); ++i)
	{
//todo: this behaviour is not well documented. it is not even clear if it is a good behaviour!
		int localSrcID = layout.get_local_src_id();
		if(localSrcID == -1)
			localSrcID = i;
			
		extractor.begin_level_extraction(i);
		for(typename Layout::iterator li = layout.begin(i);
			li != layout.end(i); ++li)
		{
			extractor.extract(*streamPack.get_stream(localSrcID),
								layout.interface(li));
		}
	}
}
				
////////////////////////////////////////////////////////////////////////
template <class TLayout>
bool SerialCommunicator<TLayout>::
communicate()
{
//	call the extractors with the received data
	for(typename ExtractorInfoList::iterator iter = m_extractorInfos.begin();
		iter != m_extractorInfos.end(); ++iter)
	{
		ExtractorInfo& info = *iter;
		if(info.m_srcProc > -1)
		{
		//	extract the data for single proc
			info.m_extractor->extract(*m_streamPackMap[info.m_interface->get_local_src_id()].
											get_stream(info.m_srcProc),
										*info.m_interface);
		}
		else
		{
		//	extract the data for a layout
		//	notify the extractor that extraction for a layout begins.
			info.m_extractor->begin_layout_extraction(info.m_layout);
			
		//	extract the data
			extract_data(*info.m_layout,
						m_streamPackMap[info.m_srcProc],
						*info.m_extractor);
			
		//	notify the extractor that extraction is complete.
			info.m_extractor->end_layout_extraction();
		}
	}

//	clean up
	m_streamPackMap.clear();
	m_extractorInfos.clear();
}

}//	end of namespace pcl

#endif
