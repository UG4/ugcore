/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__PCL__PCL_INTERFACE_COMMUNICATOR__
#define __H__PCL__PCL_INTERFACE_COMMUNICATOR__

#include <map>
#include <set>
#include <cassert>

#include "common/util/binary_buffer.h"

#include "pcl_communication_structs.h"
#include "pcl_process_communicator.h"

namespace pcl
{

/// \addtogroup pcl
/// \{

////////////////////////////////////////////////////////////////////////
//	There are two types of communicators:
//	- InterfaceCommunicator: Data is exchanged between elements of
//				interfaces / layouts. CommunicationPolicies are used to
//				collect and to extract interface-data.
//				Benefits are low communication overhead and ease of use.
//
//	- ProcessCommunicator: Arbitrary data is exchanged between all processes
//				that are linked by layouts.
//				Benefits are high flexibility and independency of
//				interface element order.
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
//	InterfaceCommunicator
///	Performs communication between interfaces on different processes.
template <typename TLayout>
class InterfaceCommunicator
{
	public:
		using Layout = TLayout;
		using Interface = typename Layout::Interface;
		using Type = typename Layout::Type;

	protected:
		using CommPol = ICommunicationPolicy<Layout>;
		
	public:
		InterfaceCommunicator();
		
	////////////////////////////////
	//	SEND
	
	///	sends raw data to a target-proc.
	/**	Schedules the data in pBuff to be sent to the target-proc
	 *	pBuff can be reused or cleared directly after the call returns.
	 *
	 *	Data sent with this method can be received using receive_raw.
	 *
	 *	Please note that this method should only be used if custom data
	 *	should be sent in a block with data that is communicated through
	 *	interfaces, since an additional copy-operation at the target process
	 *	has to be performed.
	 *	If you're only interested in sending raw data, you should take a
	 *	look into pcl::ProcessCommunicator::send.
	 */
		void send_raw(int targetProc, const void* pBuff, int bufferSize,
					  bool bSizeKnownAtTarget = false);

	///	collects data that will be sent during communicate.
	/**	Calls ICommunicationPolicy<TLayout>::collect with the specified
	 *	interface and the binary stream that is associated with the
	 *	specified target process.
	 *	Note that data will not be sent until communicate has been called.
	 *	\sa receive_data, exchange_data*/
		void send_data(int targetProc,
					   const Interface& interface,
					   ICommunicationPolicy<TLayout>& commPol);

	///	collects data that will be sent during communicate.
	/**	Calls ICommunicationPolicy<TLayout>::collect with the specified
	 *	layout and the binary stream that is associated with the
	 *	layouts target processes.
	 *	Note that data will not be sent until communicate has been called.
	 *	\sa receive_data, exchange_data*/
		void send_data(const Layout& layout,
					   ICommunicationPolicy<TLayout>& commPol);

	////////////////////////////////
	//	RECEIVE
	
	///	registers a binary-stream to receive data from a source-proc.
	/**	Receives data that has been sent with send_raw.
	 *	If send_raw was called with bSizeKnownAtTarget == true, an
	 *	exact bufferSize has to be specified. If bSizeKnownAtTarget was
	 *	set to false, bufferSize has to be set to -1.
	 *
	 *	Make sure that binStreamOut exists until communicate has been
	 *	executed.
	 *
	 *	Please note that this method should only be used if custom data
	 *	should be sent in a block with data that is communicated through
	 *	interfaces, since an additional copy-operation has to be performed.
	 *	If you're only interested in sending raw data, you should take a
	 *	look into pcl::ProcessCommunicator::receive.
	 */
		void receive_raw(int srcProc, ug::BinaryBuffer& bufOut,
						 int bufSize = -1);

	///	registers a buffer to receive data from a source-proc.
	/**	Receives data that has been sent with send_raw.
	 *	This method may only be used if send_raw was called with
	 *	bSizeKnownAtTarget == true. Call receive_raw with a binary-
	 *	stream instead, if buffer-sizes are not known.
	 *
	 *	Make sure that the buffer points to valid memory until
	 *	communicate has been executed.
	 *
	 *	Please note that this method should only be used if custom data
	 *	should be sent in a block with data that is communicated through
	 *	interfaces, since an additional copy-operation has to be performed.
	 *	If you're only interested in sending raw data, you should take a
	 *	look into pcl::ProcessCommunicator::receive.
	 */
		void receive_raw(int srcProc, void* bufOut, int bufSize);
		
	///	registers a communication-policy to receive data on communicate.
	/**	Receives have to be registered before communicate is executed.
		make sure that your instance of the communication-policy
		exists until communicate has benn executed.*/
		void receive_data(int srcProc,
						  const Interface& interface,
						  ICommunicationPolicy<TLayout>& commPol);

	///	registers an communication-policy to receive data on communicate.
	/**	Receives have to be registered before communicate is executed.
		make sure that your instance of the communication-policy
		exists until communicate has benn executed.*/
		void receive_data(const Layout& layout,
						  ICommunicationPolicy<TLayout>& commPol);

	////////////////////////////////
	//	EXCHANGE
	///	internally calls send_data and receive_data with the specified layouts.
	/**	Note that data is not communicated until communicate has been called.
	 *
	 *	TLayout has to feature the following typedefs and methods:
	 *	\code
	 *	// The type of the key with which a layout can be identified.
	 *	Key
	 *
	 *	// returns true, if the layout exists.
	 *	template <typename TType>
	 *	bool has_layout(const TLayoutMap::Key& key);
	 *	
	 *	// returns the layout that is associated with the given key.
	 *	template <typename TType>
	 *	TLayout& get_layout(const TLayoutMap::Key& key);
	 *	\endcode
	 *
	 *	The methods will only be called with type InterfaceCommunicator::Type. 
	 *
	 *	This method is particularly useful if you categorize layouts on a
	 *	process. If you separate your layouts into master and slave layouts,
	 *	you could use this method e.g. to copy data from all master-layouts
	 *	to all slave-layouts of a type with a single call.*/
		template <typename TLayoutMap>
		void exchange_data(const TLayoutMap& layoutMap,
						   const typename TLayoutMap::Key& keyFrom,
						   const typename TLayoutMap::Key& keyTo,
						   ICommunicationPolicy<TLayout>& commPol);
	
	///	sends and receives the collected data.
	/**	The collected data will be sent to the associated processes.
	 *	The extract routines of the communication-policies which were registered
	 *	through Communicator::receive_data will be called with the received
	 *	data. After all received data is processed, the communication-policies are
	 *	released. Make sure that you will keep your communication-policies
	 *	in memory until this point.
	 *	\note	Calling communicate() is effectively the same as calling
	 *			communicate_and_resume() directly followed by wait().
	 *	\param tag	For internal communications the provided tag is used. The
	 *				default value is fine in most cases and normally only has to be
	 *				adjusted if one performs multiple communications at the same time
	 *				(e.g. with different communicators).*/
		bool communicate(int tag = 749345);
		

	///	collects data and communicates it with other processes without waiting for receive
	/**	Data will be collected from the communication-policies registered in send(...)
	 * and communicated to other processes. The method however won't wait until the
	 * data has arrived at its target processes and thus can't extract the received
	 * data. This has to be done manually by calling wait().
	 * \note	Instead of using communicate_and_resume() and wait() you could
	 * 			simply call communicate. Separating communication and wait however
	 * 			gives you the benefit of being able to continue with other calculations
	 * 			while communication is performed.
	 * \note	A call to communicate_and_resume() has to be followed by a call to wait().
	 * 			You may not call communicate_and_resume() twice on a single communicator
	 * 			without calling wait() in between.
	 * \param tag	For internal communications the provided tag is used. The
	 *				default value is fine in most cases and normally only has to be
	 *				adjusted if one performs multiple communications at the same time
	 *				(e.g. with different communicators).*/
		bool communicate_and_resume(int tag = 749345);

	///	waits for the data communicated by communicate_and_resume() and extracts it
	/**	The extract routines of the communication-policies which were registered
	 *	through Communicator::receive_data will be called with the received
	 *	data. After all received data is processed, the communication-policies are
	 *	released. Make sure that you will keep your communication-policies
	 *	in memory until this point.*/
		void wait();
	

	///	enables debugging of communication. This has a severe effect on performance!
	/**	communication debugging will execute some code during communicate(), which
	 * checks whether matching sends and receives have been scheduled with matching
	 * buffer sizes.
	 *
	 * If not all processes participate during communication, you have to specify
	 * a process-communicator (involvedProcs), which includes all and only the procs,
	 * which will call communicate later on. If no process communicator is specified,
	 * then PCD_WORLD is used.
	 *
	 * Note that communication debugging introduces additional communication which
	 * considerable slows down performance. You should only use it temporarily, if
	 * you encounter problems during communication.
	 *
	 * Don't forget to call disable_communication_debugging(), when you're done
	 * with debugging.*/
		void enable_communication_debugging(const ProcessCommunicator& involvedProcs =
												ProcessCommunicator(PCD_WORLD));
	 
	///	disables debugging of communication
		void disable_communication_debugging();
	
	///	returns true if communication debugging is enabled
		[[nodiscard]] bool communication_debugging_enabled() const;
	 
	protected:
		using BufferMap = std::map<int, ug::BinaryBuffer>;

	protected:
	///	helper to collect data from single-level-layouts
		void send_data(const Layout& layout,
				  ICommunicationPolicy<TLayout>& commPol,
				  const layout_tags::single_level_layout_tag&);

	///	helper to collect data from multi-level-layouts				  
		void send_data(const Layout& layout,
				  	   ICommunicationPolicy<TLayout>& commPol,
				  	   const layout_tags::multi_level_layout_tag&);
	
	///	prepare stream-pack-in
		void prepare_receiver_buffer_map(BufferMap& bufMap,
										 std::set<int>& curProcs,
										 const TLayout& layout);
	/// specialization of stream-pack preparation for single-level-layouts
		void prepare_receiver_buffer_map(BufferMap& bufMap,
										 std::set<int>& curProcs,
										 const TLayout& layout,
										 const layout_tags::single_level_layout_tag&);
	/// specialization of stream-pack preparation for multi-level-layouts
		void prepare_receiver_buffer_map(BufferMap& bufMap,
										 std::set<int>& curProcs,
										 const TLayout& layout,
										 const layout_tags::multi_level_layout_tag&);

	///	collects buffer sizes for a given layout and stores them in a map
	/**	The given map holds pairs of procID, bufferSize
	 *	If buffer-sizes can't be determined, false is returned.
	 *	if pMmapBuffSizesOut == nullptr, the method will simply determine
	 *	whether all buffer-sizes can be calculated.*/
	 	bool collect_layout_buffer_sizes(const TLayout& layout,
										 ICommunicationPolicy<TLayout>& commPol,
										 std::map<int, int>* pMapBuffSizesOut,
										 const layout_tags::single_level_layout_tag&);

	///	collects buffer sizes for a given layout and stores them in a map
	/**	The given map holds pairs of procID, bufferSize
	 *	If buffer-sizes can't be determined, false is returned.
	 *	if pMmapBuffSizesOut == nullptr, the method will simply determine
	 *	whether all buffer-sizes can be calculated.*/
	 	bool collect_layout_buffer_sizes(const TLayout& layout,
										 ICommunicationPolicy<TLayout>& commPol,
										 std::map<int, int>* pMapBuffSizesOut,
										 const layout_tags::multi_level_layout_tag&);
	
	///	extract data from stream-pack
		void extract_data(const TLayout& layout, BufferMap& bufMap,
						  CommPol& extractor);
		
		void extract_data(const TLayout& layout, BufferMap& bufMap,
						  CommPol& extractor,
						  const layout_tags::single_level_layout_tag&);
		
		void extract_data(const TLayout& layout, BufferMap& bufMap,
						  CommPol& extractor,
						  const layout_tags::multi_level_layout_tag&);
		
	protected:		
	///	holds information that will be passed to the extract routines.
	/**	if srcProc == -1, the layout will be used for extraction.
	 *	if srcProc >= 0, either the buffer, the binaryStream or the
	 *	interface will be used for extraction, depending on which is
	 *	not nullptr.
	 */
		struct ExtractorInfo
		{
			ExtractorInfo() = default;
			ExtractorInfo(int srcProc, CommPol* pExtractor,
						const Interface* pInterface, const Layout* pLayout,
						void* buffer, ug::BinaryBuffer* binBuffer, int rawSize) :
				m_srcProc(srcProc), m_extractor(pExtractor),
				m_interface(pInterface), m_layout(pLayout),
				m_buffer(buffer), m_binBuffer(binBuffer), m_rawSize(rawSize)
			{
				assert((srcProc == -1) || ((srcProc >= 0) && (srcProc < pcl::NumProcs())));
			}

			int					m_srcProc;
			CommPol*			m_extractor;			
			const Interface*	m_interface;
			const Layout*		m_layout;
			void*				m_buffer;
			ug::BinaryBuffer*	m_binBuffer;
			int					m_rawSize;
		};

	///	A list that holds information about extractors.
		using ExtractorInfoList = std::list<ExtractorInfo>;

	protected:
	///	holds the buffers that are used to send data
		BufferMap m_bufMapOut;
	///	stores out-procs for the next communication step
		std::set<int> m_curOutProcs;

	///	holds the buffers that are used to receive data
		BufferMap m_bufMapIn;
	///	stores in-procs for the next communication step
		std::set<int> m_curInProcs;

	///	holds information about the extractors that are awaiting data.
		ExtractorInfoList m_extractorInfos;
		
	///	used by communicate, communicate_and_resume and wait, to check whether communication is done.
		std::vector<MPI_Request> m_vSendRequests;
		std::vector<MPI_Request> m_vReceiveRequests;

	///	This procComm holds the processes that shall participate during communication-debugging.
		ProcessCommunicator	m_debugProcComm;
		
	///	This is the tag for the currently performed communication.
	/**	Set to -1 if no communication is currently performed.*/
		int	m_curComTag;

	///	true if the communication shall be debugged.
		bool m_bDebugCommunication;
		
	///	holds info whether all send-buffers are of predetermined fixed size.
	/**	reset to true after each communication-step.*/
		bool m_bSendBuffersFixed;
};

// end group pcl
/// \}

}//	end of namespace pcl

////////////////////////////////////////
//	include implementation
#include "pcl_interface_communicator_impl.hpp"

#endif
