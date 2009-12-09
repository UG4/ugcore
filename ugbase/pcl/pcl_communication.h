//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

#ifndef __H__PCL__PCL_COMMUNICATION__
#define __H__PCL__PCL_COMMUNICATION__

#include <iostream>
#include <map>
#include "common/util/binary_stream.h"
#include "common/util/stream_pack.h"
#include "pcl_base.h"

namespace pcl
{

////////////////////////////////////////////////////////////////////////
//	ICollector
///	interface for collectors.
template <class TElementGroup>
class ICollector
{
	public:
	//	typedefs
		typedef typename group_traits<TElementGroup>::Interface Interface;

	public:
	///	writes the data which is contained in the interface to the buffer.
		virtual bool
		collect(std::ostream& buff, Interface& interface) = 0;
};

////////////////////////////////////////////////////////////////////////
//	IExtractor
///	interface for extractors.
template <class TElementGroup>
class IExtractor
{
	public:
	//	typedefs
		typedef typename group_traits<TElementGroup>::Layout Layout;
		typedef typename group_traits<TElementGroup>::Interface Interface;

	public:
	///	signals the beginning of a layout extraction.
	/**	the default implementation returns true and does nothing else.*/
		virtual bool
		begin_layout_extraction(Layout* pLayout)	{return true;}

	///	signals the end of a layout extraction
	/**	the default implementation returns true and does nothing else.*/
		virtual bool
		end_layout_extraction()						{return true;}

	///	extracts the data from the buffer and assigns it to the interface-elements.
	/**	If this method is called between calls to begin_layout_extraction and
		end_layout_extraction, the interface that is passed to this method
		belongs to the layout.
	*/
		virtual bool
		extract(std::istream& buff, Interface& interface) = 0;
};

////////////////////////////////////////////////////////////////////////
//	Communicator
///	Performs communication between processes.
template <class TElementGroup>
class Communicator
{
	public:
	//	typedefs
		typedef typename group_traits<TElementGroup>::Layout Layout;
		typedef typename group_traits<TElementGroup>::Interface Interface;
		typedef IExtractor<TElementGroup> Extractor;

	public:
	///	calls the collector with the binary stream that is associated with targetProc.
		void collect_data(int targetProc,
						  Interface& interface,
						  ICollector<TElementGroup>& collector);

	///	calls the collector with the binary stream that is associated with targetProc.
	/**	The interfaces of the layout are passed to the collector on after the other.*/
		void collect_data(Layout& layout,
						  ICollector<TElementGroup>& collector);

	///	registers an extractor to receive data on communicate.
	/**	make sure that your instance of the extractor exists until
		communicate has benn executed.*/
		void await_data(int srcProc,
						Interface& interface,
						IExtractor<TElementGroup>& extractor);

	///	registers an extractor to receive data on communicate.
	/**	make sure that your instance of the extractor exists until
		communicate has benn executed.*/
		void await_data(Layout& layout,
						IExtractor<TElementGroup>& extractor);

	///	sends and receives data the collected data.
	/**	The collected data will be send to the associated processes.
	 *	The extract routines of the extractors which were registered
	 *	through Communicator::await_data will be called with the received
	 *	data. After all received data is processed, the extractors are
	 *	freed. Make sure that you will keep your extractors in memory
	 *	until this point.*/
		bool communicate();

	protected:
	///	holds information that will be passed to the extract routines.
	/**	if srcProc == -1, the layout will be used for extraction.
	 *	if srcProc >= 0, the srcProc and the interface will be used.*/
		struct ExtractorInfo
		{
			ExtractorInfo()			{}
			ExtractorInfo(Extractor* pExtractor, int srcProc,
						Interface* pInterface, Layout* pLayout) : m_extractor(pExtractor), m_srcProc(srcProc),
																m_interface(pInterface), m_layout(pLayout)	{}

			Extractor*	m_extractor;
			int			m_srcProc;
			Interface*	m_interface;
			Layout*		m_layout;
		};

	///	A list that holds information about extractors.
		typedef std::list<ExtractorInfo> ExtractorInfoList;

	protected:
	///	holds the streams that are used to send data
		ug::StreamPack		m_streamPackOut;
	///	holds the streams that are used to receive data
		ug::StreamPack		m_streamPackIn;
	///	holds information about the extractors that are awaiting data.
		ExtractorInfoList	m_extractorInfos;
};

}//	end of namespace pcl

////////////////////////////////////////
//	include implementation
#include "pcl_communication_impl.hpp"

#endif
