//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

#ifndef __H__PCL__PCL_COMMUNICATOR__
#define __H__PCL__PCL_COMMUNICATOR__

#include <iostream>
#include <map>
#include "common/util/binary_stream.h"
#include "common/util/stream_pack.h"
#include "cl_base.h"

namespace pcl
{

////////////////////////////////////////////////////////////////////////
//	ParallelCommunicator
///	Performs communication between interfaces on different processes.
template <class TLayout>
class ParallelCommunicator
{
	public:
	//	typedefs
		typedef TLayout 					Layout;
		typedef typename Layout::Interface	Interface;
		typedef typename Layout::Type		Type;

	protected:
		typedef ICommunicationPolicy<Layout>	CommPol;
		
	public:
	////////////////////////////////
	//	COLLECT
	///	collects data that will be send during communicate.
	/**	Calls ICommunicationPolicy<TLayout>::collect with the specified
	 *	interface and the binary stream that is associated with the
	 *	specified target process.
	 *	Note that data will not be send until communicate has been called.
	 *	\sa receive_data, exchange_data*/
		void send_data(int targetProc,
						  Interface& interface,
						  ICommunicationPolicy<TLayout>& commPol);

	///	collects data that will be send during communicate.
	/**	Calls ICommunicationPolicy<TLayout>::collect with the specified
	 *	layout and the binary stream that is associated with the
	 *	layouts target processes.
	 *	Note that data will not be send until communicate has been called.
	 *	\sa receive_data, exchange_data*/
		void send_data(Layout& layout,
						  ICommunicationPolicy<TLayout>& commPol);

	////////////////////////////////
	//	AWAIT
	///	registers a communication-policy to receive data on communicate.
	/**	make sure that your instance of the communication-policy
		exists until communicate has benn executed.*/
		void receive_data(int srcProc,
						Interface& interface,
						ICommunicationPolicy<TLayout>& commPol);

	///	registers an communication-policy to receive data on communicate.
	/**	make sure that your instance of the communication-policy
		exists until communicate has benn executed.*/
		void receive_data(Layout& layout,
						ICommunicationPolicy<TLayout>& commPol);

	////////////////////////////////
	//	EXCHANGE
	///	internally calls send_data and receive_data with the specified layouts.
	/**	Note that data is not communicated until communicate has been called.
	 *
	 *	Make sure that the Layout- and the Interface-type of TLayoutMap
	 *	are compatible with the layout of the communicator.
	 *	Layouts are queried for TLayout::Type of the communicators TLayout-type.
	 *
	 *	This method is particularily useful if you categorize layouts on a
	 *	process. If you separate your layouts into master and slave layouts,
	 *	you could use this method e.g. to copy data from all master-layouts
	 *	to all slave-layouts of a type with a single call.*/
		template <class TLayoutMap>
		void exchange_data(TLayoutMap& layoutMap,
							typename TLayoutMap::Key keyFrom,
							typename TLayoutMap::Key keyTo,
							ICommunicationPolicy<TLayout>& commPol);
			
	///	sends and receives the collected data.
	/**	The collected data will be send to the associated processes.
	 *	The extract routines of the communication-policies which were registered
	 *	through Communicator::await_data will be called with the received
	 *	data. After all received data is processed, the communication-policies are
	 *	released. Make sure that you will keep your communication-policies
	 *	in memory until this point.*/
		bool communicate();
		
	protected:
	///	helper to collect data from single-level-layouts
		void send_data(Layout& layout,
				  ICommunicationPolicy<TLayout>& commPol,
				  const layout_tags::single_level_layout_tag&);

	///	helper to collect data from multi-level-layouts				  
		void send_data(Layout& layout,
				  ICommunicationPolicy<TLayout>& commPol,
				  const layout_tags::multi_level_layout_tag&);
	
	///	prepare stream-pack-in
		void prepare_receiver_stream_pack(ug::StreamPack& streamPack,
											TLayout& layout);
	/// specialization of stream-pack preparation for single-level-layouts
		void prepare_receiver_stream_pack(ug::StreamPack& streamPack,
										TLayout& layout,
										const layout_tags::single_level_layout_tag&);
	/// specialization of stream-pack preparation for multi-level-layouts
		void prepare_receiver_stream_pack(ug::StreamPack& streamPack,
										TLayout& layout,
										const layout_tags::multi_level_layout_tag&);

	///	extract data from stream-pack
		void extract_data(TLayout& layout, ug::StreamPack& streamPack,
						CommPol& extractor);
		
		void extract_data(TLayout& layout, ug::StreamPack& streamPack,
						CommPol& extractor,
						const layout_tags::single_level_layout_tag&);
		
		void extract_data(TLayout& layout, ug::StreamPack& streamPack,
						CommPol& extractor,
						const layout_tags::multi_level_layout_tag&);
		
	protected:		
	///	holds information that will be passed to the extract routines.
	/**	if srcProc == -1, the layout will be used for extraction.
	 *	if srcProc >= 0, the srcProc and the interface will be used.*/
		struct ExtractorInfo
		{
			ExtractorInfo()			{}
			ExtractorInfo(CommPol* pExtractor, int srcProc,
						Interface* pInterface, Layout* pLayout) : m_extractor(pExtractor), m_srcProc(srcProc),
																m_interface(pInterface), m_layout(pLayout)	{}

			CommPol*	m_extractor;
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
#include "pcl_communicator_impl.hpp"

#endif
