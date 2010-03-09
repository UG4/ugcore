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
//	ICommunicationPolicy
///	specializations are responsible to pack and unpack interface data during communication.
template <class TLayout>
class ICommunicationPolicy
{
	public:
		typedef TLayout						Layout;
		typedef typename Layout::Interface 	Interface;
		
	////////////////////////////////
	//	COLLECT
	///	should write data which is associated with the interface elements to the buffer.
		virtual bool
		collect(std::ostream& buff, Interface& interface) = 0;
		
	////////////////////////////////
	//	EXTRACT
	///	signals the beginning of a layout extraction.
	/**	the default implementation returns true and does nothing else.*/
		virtual bool
		begin_layout_extraction(Layout* pLayout)	{return true;}

	///	signals the end of a layout extraction
	/**	the default implementation returns true and does nothing else.*/
		virtual bool
		end_layout_extraction()						{return true;}

	///	extract data from the buffer and assigns it to the interface-elements.
	/**	If this method is called between calls to begin_layout_extraction and
		end_layout_extraction, the interface that is passed to this method
		belongs to the layout.*/
		virtual bool
		extract(std::istream& buff, Interface& interface) = 0;
};

////////////////////////////////////////////////////////////////////////
//	Communicator
///	Performs communication between processes.
template <class TLayout>
class Communicator
{
	public:
	//	typedefs
		typedef TLayout 					Layout;
		typedef typename Layout::Interface	Interface;

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
	///	internally calls collect_data and await data with the specified interface.
	/**	Note that data is not communicated until communicate has been called.*/
	/*
		void exchange_data(Interface& interface,
							ICommunicationPolicy<TLayout>& commPol);
	*/

	///	internally calls collect_data and await data with the specified layout.
	/**	Note that data is not communicated until communicate has been called.*/
	/*
		void exchange_data(Layout& layout,
							ICommunicationPolicy<TLayout>& commPol);
	*/
	
	///	sends and receives the collected data.
	/**	The collected data will be send to the associated processes.
	 *	The extract routines of the communication-policies which were registered
	 *	through Communicator::await_data will be called with the received
	 *	data. After all received data is processed, the communication-policies are
	 *	released. Make sure that you will keep your communication-policies
	 *	in memory until this point.*/
		bool communicate();

	protected:
		typedef ICommunicationPolicy<Layout>	CommPol;
		
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
#include "pcl_communication_impl.hpp"

#endif
