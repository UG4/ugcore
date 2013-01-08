//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

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
template <class TLayout>
class InterfaceCommunicator
{
	public:
	//	typedefs
		typedef TLayout 					Layout;
		typedef typename Layout::Interface	Interface;
		typedef typename Layout::Type		Type;

	protected:
		typedef ICommunicationPolicy<Layout>	CommPol;
		
	public:
		InterfaceCommunicator();
		
	////////////////////////////////
	//	SEND
	
	///	sends raw data to a target-proc.
	/**	Shedules the data in pBuff to be sent to the target-proc
	 *	pBuff can be reused or cleared directly after the call returns.
	 *
	 *	Data sent with this method can be received using receive_raw.
	 *
	 *	Please note that this method should only be used if custom data
	 *	should be send in a block with data that is communicated through
	 *	interfaces, since an additional copy-operation at the target process
	 *	has to be performed.
	 *	If you're only interested in sending raw data, you should take a
	 *	look into pcl::ProcessCommunicator::send.
	 */
		void send_raw(int targetProc, const void* pBuff, int bufferSize,
					   bool bSizeKnownAtTarget = false);

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
	 *	should be send in a block with data that is communicated through
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
	 *	should be send in a block with data that is communicated through
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
						Interface& interface,
						ICommunicationPolicy<TLayout>& commPol);

	///	registers an communication-policy to receive data on communicate.
	/**	Receives have to be registered before communicate is executed.
		make sure that your instance of the communication-policy
		exists until communicate has benn executed.*/
		void receive_data(Layout& layout,
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
	 *	template <class TType>
	 *	bool has_layout(const TLayoutMap::Key& key);
	 *	
	 *	// returns the layout that is associated with the given key.
	 *	template <class TType>
	 *	TLayout& get_layout(const TLayoutMap::Key& key);
	 *	\endcode
	 *
	 *	The methods will only be called with type InterfaceCommunicator::Type. 
	 *
	 *	This method is particularily useful if you categorize layouts on a
	 *	process. If you separate your layouts into master and slave layouts,
	 *	you could use this method e.g. to copy data from all master-layouts
	 *	to all slave-layouts of a type with a single call.*/
		template <class TLayoutMap>
		void exchange_data(TLayoutMap& layoutMap,
							const typename TLayoutMap::Key& keyFrom,
							const typename TLayoutMap::Key& keyTo,
							ICommunicationPolicy<TLayout>& commPol);
	
	///	sends and receives the collected data.
	/**	The collected data will be send to the associated processes.
	 *	The extract routines of the communication-policies which were registered
	 *	through Communicator::await_data will be called with the received
	 *	data. After all received data is processed, the communication-policies are
	 *	released. Make sure that you will keep your communication-policies
	 *	in memory until this point.*/
		bool communicate();
		
	
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
		bool communication_debugging_enabled();
	 
	protected:
		typedef std::map<int, ug::BinaryBuffer>	BufferMap;

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
		void prepare_receiver_buffer_map(BufferMap& bufMap,
											std::set<int>& curProcs,
											TLayout& layout);
	/// specialization of stream-pack preparation for single-level-layouts
		void prepare_receiver_buffer_map(BufferMap& streamPack,
										std::set<int>& curProcs,
										TLayout& layout,
										const layout_tags::single_level_layout_tag&);
	/// specialization of stream-pack preparation for multi-level-layouts
		void prepare_receiver_buffer_map(BufferMap& streamPack,
										std::set<int>& curProcs,
										TLayout& layout,
										const layout_tags::multi_level_layout_tag&);

	///	collects buffer sizes for a given layout and stores them in a map
	/**	The given map holds pairs of procID, bufferSize
	 *	If buffer-sizes can't be determined, false is returned.
	 *	if pMmapBuffSizesOut == NULL, the method will simply determine
	 *	whether all buffersizes can be calculated.*/
	 	bool collect_layout_buffer_sizes(TLayout& layout,
										ICommunicationPolicy<TLayout>& commPol,
										std::map<int, int>* pMapBuffSizesOut,
										const layout_tags::single_level_layout_tag&);

	///	collects buffer sizes for a given layout and stores them in a map
	/**	The given map holds pairs of procID, bufferSize
	 *	If buffer-sizes can't be determined, false is returned.
	 *	if pMmapBuffSizesOut == NULL, the method will simply determine
	 *	whether all buffersizes can be calculated.*/
	 	bool collect_layout_buffer_sizes(TLayout& layout,
										ICommunicationPolicy<TLayout>& commPol,
										std::map<int, int>* pMapBuffSizesOut,
										const layout_tags::multi_level_layout_tag&);	
	
	///	extract data from stream-pack
		void extract_data(TLayout& layout, BufferMap& bufMap,
						CommPol& extractor);
		
		void extract_data(TLayout& layout, BufferMap& bufMap,
						CommPol& extractor,
						const layout_tags::single_level_layout_tag&);
		
		void extract_data(TLayout& layout, BufferMap& bufMap,
						CommPol& extractor,
						const layout_tags::multi_level_layout_tag&);
		
	protected:		
	///	holds information that will be passed to the extract routines.
	/**	if srcProc == -1, the layout will be used for extraction.
	 *	if srcProc >= 0, either the buffer, the binaryStream or the
	 *	interace will be used for extraction, depending on which is
	 *	not NULL.
	 */
		struct ExtractorInfo
		{
			ExtractorInfo()			{}
			ExtractorInfo(int srcProc, CommPol* pExtractor,
						Interface* pInterface, Layout* pLayout,
						void* buffer, ug::BinaryBuffer* binBuffer, int rawSize) :
				m_srcProc(srcProc), m_extractor(pExtractor),
				m_interface(pInterface), m_layout(pLayout),
				m_buffer(buffer), m_binBuffer(binBuffer), m_rawSize(rawSize)
			{
				assert((srcProc == -1) || ((srcProc >= 0) && (srcProc < pcl::GetNumProcesses())));
			}

			int					m_srcProc;
			CommPol*			m_extractor;			
			Interface*			m_interface;
			Layout*				m_layout;
			void*				m_buffer;
			ug::BinaryBuffer*	m_binBuffer;
			int					m_rawSize;
		};

	///	A list that holds information about extractors.
		typedef std::list<ExtractorInfo> ExtractorInfoList;

	protected:
	///	holds the buffers that are used to send data
		BufferMap		m_bufMapOut;
	///	stores out-procs for the next communication step
		std::set<int>	m_curOutProcs;

	///	holds the buffers that are used to receive data
		BufferMap		m_bufMapIn;
	///	stores in-procs for the next communication step
		std::set<int>	m_curInProcs;

	///	holds information about the extractors that are awaiting data.
		ExtractorInfoList	m_extractorInfos;
		
	///	This procComm holds the processes that shall participate during communication-debugging.
		ProcessCommunicator	m_debugProcComm;
		
	///	true if the communication shall be debugged.
		bool m_bDebugCommunication;
		
	///	holds info whether all send-buffers are of predetermined fixed size.
	/**	reset to true after each communication-step.*/
		bool m_bSendBuffersFixed;
};

}//	end of namespace pcl

////////////////////////////////////////
//	include implementation
#include "pcl_interface_communicator_impl.hpp"

#endif
