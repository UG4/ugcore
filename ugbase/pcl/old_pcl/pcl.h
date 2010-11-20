//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m05 d08

#ifndef __H__PCL__
#define __H__PCL__

#include <vector>
#include <list>
#include <iostream>
#include "common/util/hash.h"
#include "common/util/section_container.h"
#include "common/util/binary_stream.h"

//	Don't rely on mpi being included.
//	It is only included to allow us to define some constants.
//	This include will most likely be removed in future versions.
#include "mpi.h"

////////////////////////////////////////////////////////////////////////
///	this allows us to print messages to the users terminal
/**
 * if an output-processor is specified through \sa pcl::SetOutputProcRank
 * only messages from that processor will be printed to the screen.
 */
#define PCLLOG(msg) if((pcl::GetProcRank() == pcl::GetOutputProcRank()) || pcl::GetOutputProcRank() == -1) {std::cout << msg; std::cout.flush();}

namespace pcl
{

typedef int ProcID;

//	be careful: used as section indices, too.
enum NodeType
{
	NT_NORMAL = 0,
	NT_MASTER = 1,
	NT_SLAVE = 2,
};

const int NT_NUM_TYPES = 3;

////////////////////////////////////////////////////////////////////////
//	operations for get_data
enum DataOperation
{
	DO_FIRST = 1,
	DO_ADD,
	DO_MIN,
	DO_MAX
};

//	ReduceOperation
#define RO_MAX MPI_MAX
#define RO_MIN MPI_MIN
#define RO_SUM MPI_SUM
#define RO_PROD MPI_PROD
#define RO_LAND MPI_LAND
#define RO_BAND MPI_BAND
#define RO_LOR MPI_LOR
#define RO_BOR MPI_BOR
#define RO_LXOR MPI_LXOR
#define RO_BXOR MPI_BXOR
#define RO_MAXLOC MPI_MAXLOC
#define RO_MINLOC MPI_MINLOC

typedef MPI_Op ReduceOperation;

//	ReduceOperation
#define DT_BYTE MPI_BYTE,
#define DT_PACKED MPI_PACKED
#define DT_CHAR MPI_CHAR
#define DT_SHORT MPI_SHORT
#define DT_INT MPI_INT
#define DT_LONG MPI_LONG
#define DT_FLOAT MPI_FLOAT
#define DT_DOUBLE MPI_DOUBLE
#define DT_LONG_DOUBLE MPI_LONG_DOUBLE
#define DT_UNSIGNED_CHAR MPI_UNSIGNED_CHAR

typedef MPI_Datatype DataType;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

///	call this method before any other pcl-operations.
void Init(int argc, char* argv[]);
///	call this method right before quitting your application
void Finalize();

///	returns the number of processes
int GetNumProcesses();
///	returns the rank of the process
ProcID GetProcRank();
///	returns the rank of the process for which the output shall be printed to the screen.
ProcID GetOutputProcRank();
///	sets the rank of the process for which the output shall be printed.
void SetOutputProcRank(ProcID rank);

///	returns the time in seconds
double Time();

///	sends data to another process. data may be received using \sa ReceiveData or \sa CollectData.
void SendData(ProcID destProc, void* pBuffer, int bufferSize, int tag);

///	receives the data that was send with \sa SendData or \sa DistributeData.
void ReceiveData(void* pBuffOut, ProcID srcProc, int bufferSize, int tag);

///	collect the data send with send_data from proc firstSendProc to numSendProcs excluding destProc.
void CollectData(ProcID thisProcID, int firstSendProc, int numSendProcs,
					void* pBuffer, int bufferSizePerProc, int tag);

///	sends the data in the different sections of the buffer to the specified processes.
/**
 * pBufferSegSizes has to contain numRecProcs elements.
 * Buffer-Segments are send to the processes
 * firstRecProc, ..., firstRecProc + numRecProcs.
 * If thisProcID lies in this range, then the buffer-segments are
 * sent to firstRecProc, ..., firstRecProc + numRecProcs + 1,
 * excluding thisProcID.
 */
void DistributeData(ProcID thisProcID, int firstRecProc, int numRecProcs,
					void* pBuffer, int* pBufferSegSizes, int tag);

///	sends the data in the different sections of the buffer to the specified processes.
/**
 * pRecProcMap and pBufferSegSizes have to contain numRecProcs elements.
 * Entries in pRecProcMap specify the target-process of the i-th buffer
 * segment.
 */
void DistributeData(ProcID thisProcID, int* pRecProcMap, int numRecProcs,
					void* pBuffer, int* pBufferSegSizes, int tag);

///	reduces the data to a single buffer using the specified ReduceOperation and distributes the result to all processes.
/**
 * Both, sendBuf and recBuf have to hold count elements of the specified type.
 * \param op has to be one of the
 */
void AllReduce(void* sendBuf, void* recBuf, int count, DataType type,
				ReduceOperation op);
 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//
//TODO: update docu
/**
 * The CommunicationGroup is a class that represents an abstract
 * concept for cluster-communication.
 * Nodes that can be added by the user are organized in interfaces
 * over which data can be exchanged with other processes.
 * The type of the nodes can be specified by the user through the
 * template argument.
 * Please note that you may only add, remove, replace, ... nodes
 * between calls to begin_node_managment() and end_node_managment().
 * The same holds true for calls to set_node_type() and
 * add_node_to_interface().
 *
 * All nodes have a NodeType. The NodeType is set per node, not per
 * interface. The NodeType is important for communication, since
 * data can be exchanged between master and slave nodes or between
 * normal nodes or between all interface-nodes at once.
 * Make sure that for each master node on process A, that is contained
 * in the interface to process B, an associated slave node on process B
 * exists, which is contained in the interface to process A.
 * Both nodes have to have the same interface index!
 *
 * There are two ways to exchange data between nodes.
 * The fast one is the data exchange through send- and receive buffers.
 * You can create send- and receive-buffers through
 * create_send_buffer() and create_receive_buffer().
 * If you want to exchange data in the next communication step
 * you'll have to schedule those buffers through
 * schedule_send_buffer() and schedule_receive_buffer().
 * Make sure that you always schedule matching buffers on both
 * sides of an interface.
 *
 * The other way to exchange data is using schedule_personal_data().
 * This method allows you to exchange data directly between on node
 * and its associated nodes on other processes.
 *
 * Note that no communication requests are processesd until
 * communicate() is called.
 *
 * The class uses hashes internally to associate nodes and indices.
 * be sure to overload the hash-function
 * template <typename TKey> unsigned long hash_key(const TKey& key)
 * for TNode before you declare a specialised version of the
 * CommunicationGroup anywhere.
 */
template <class TNode>
class CommunicationGroup
{
	protected:
		class NodeInfo;
		class DataBuffer;
		typedef unsigned char byte;
		
	public:
		typedef TNode	node_type;
		typedef std::list<NodeInfo*>  NodeInfoList;
		typedef NodeInfo*	HNODE;
		typedef ug::Hash<HNODE, TNode>	NodeHandleHash;
		typedef ug::SectionContainer<NodeInfo*, NodeInfoList>	NodeHandleContainer;	//	all nodes of the distribution-group are stored in are stored in a section-container.
		typedef typename NodeHandleContainer::iterator	NodeHandleIterator;
		typedef std::vector<NodeInfo*> NodeHandleVec;	//	used in interfaces
		
	/**	The DataBufferInfo is used to identify a DataBuffer.*/
		struct DataBufferInfo
		{
			friend class CommunicationGroup<TNode>;
			
			protected:
				bool isInterfaceBuffer;
				bool isSendBuffer;
				unsigned int	bufferSetID;	// A bufferSetID of 0 marks the buffer as invalid.
				unsigned int	bufferIndex;
		};
		
		typedef DataBufferInfo DBHANDLE;


	public:
		CommunicationGroup();///< be sure to call set_node_hash_size if you use the default constructor.
		CommunicationGroup(uint nodeHashSize);
		CommunicationGroup(ProcID procID, uint nodeHashSize);
		CommunicationGroup(const CommunicationGroup& dg);///< The implemented copy-constructor only copies initial values (no nodes).
		~CommunicationGroup();
		
		void clear();

	/**	warning - set_node_hash_size automatically calls clear in the actual implementation.
		Be sure to call this method before starting adding any nodes to the group.*/
	//TODO: allow hash-resize without data-loss.
		void set_node_hash_size(uint nodeHashSize);
		
		inline void set_proc_id(ProcID procID)	{m_procID = procID;}///< WARNING: a call to this method will invalidate all existing interface-relations of other groups to this group. Use with care!
		inline ProcID get_proc_id()				{return m_procID;}

	////////////////////////////////
	//	node-managment
	/**	call this method before you add, erase or replace any nodes to the group or
		to interfaces. This method has to be called before the change of node-types, too.*/
		void begin_node_managment();
	/**	returns true if the group is in node-managment-mode.*/
		inline bool node_managment_is_active()	{return m_bNodeManagmentMode;}
	/**	add a node to the group. Call this method between
		begin_node_managment() and end_node_managment().*/
		HNODE add_node(const TNode& node);///< returns node_infos_end() if the operation did not succeed.
	/**	removes a node from the group. Call this method between
		begin_node_managment() and end_node_managment().
		Interfaces are not updated until end_node_managment() is called.*/
		void erase_node(HNODE hNode);
	/**	removes nodeOld and adds nodeNew without invalidating the handle.
		 Call this method between begin_node_managment() and end_node_managment().*/		
		HNODE replace_node(const TNode& nodeOld, const TNode& nodeNew);
	/**	changes the node type of a node. This affects all interfaces in which the node
		is contained. Interfaces are not updated until end_node_managment() however.
		Call this method between begin_node_managment() and end_node_managment().*/
		void set_node_type(HNODE hNode, NodeType nt);
		
	/**	Be careful: No checks are performed whether a node already exists in the interface!
	  *	returns the index of the node in the interface.*/
		int add_node_to_interface(HNODE hNode,
									ProcID connectedProcID);
												
	/**	Be careful: No checks are performed whether a node already exists in the interface!
	  *	please note, that the NodeType of a node is the same for the whole distribution-group.
	  *	Changing it here will thus affect other interfaces as well.
	  * This is why a change from master to slave or from slave to master will be fetched with
	  * an assertion in this method. If you really want to change the node-type from master to slave
	  *	or vice versa, then please use the set_node_type(...) method instead.
	  *	returns the index of the node in the interface.*/
		int add_node_to_interface(HNODE hNode,
									ProcID connectedProcID,
									NodeType nt);
	
	/**	ends node-managment mode and updates interfaces.*/
		void end_node_managment();
		
	////////////////////////////////
	//	node access
		inline int num_nodes()						{return m_nodeHandles.num_elements();}
		HNODE get_handle(const TNode& node);
		inline const TNode& get_node(HNODE hNode)							{return hNode->node;}
		inline HNODE invalid_handle()										{return NULL;}
		inline bool handle_is_invalid(HNODE hNode)							{return hNode == invalid_handle();}
		inline bool handle_is_valid(HNODE hNode)							{return hNode != invalid_handle();}
		inline NodeType get_node_type(HNODE hNode)							{return hNode->nt;}
		inline const std::vector<int>& get_interface_indices(HNODE hNode)	{return hNode->vInterfaceIndices;}
		inline const std::vector<int>& get_interface_entry_indices(HNODE hNode)	{return hNode->vInterfaceEntryIndices;}
		inline int num_containing_interfaces(HNODE hNode)					{return hNode->vInterfaceIndices.size();}

	////////////////////////////////
	//	iterators
	//  iteraters may point to invalid handles during node-managment-mode!
		inline NodeHandleIterator handles_begin()	{return m_nodeHandles.begin();}
		inline NodeHandleIterator handles_end()		{return m_nodeHandles.end();}
		inline NodeHandleIterator normal_handles_begin()	{return m_nodeHandles.section_begin(NT_NORMAL);}
		inline NodeHandleIterator normal_handles_end()		{return m_nodeHandles.section_end(NT_NORMAL);}
		inline NodeHandleIterator master_handles_begin()	{return m_nodeHandles.section_begin(NT_MASTER);}
		inline NodeHandleIterator master_handles_end()		{return m_nodeHandles.section_end(NT_MASTER);}
		inline NodeHandleIterator slave_handles_begin()		{return m_nodeHandles.section_begin(NT_SLAVE);}
		inline NodeHandleIterator slave_handles_end()		{return m_nodeHandles.section_end(NT_SLAVE);}
		inline NodeHandleIterator boundary_handles_begin()	{return m_nodeHandles.section_begin(NT_MASTER);}
		inline NodeHandleIterator boundary_handles_end()	{return m_nodeHandles.section_end(NT_SLAVE);}

	////////////////////////////////
	//	interface managment
		inline int num_interfaces()											{return m_vInterfaces.size();}
		inline ProcID get_connected_proc_of_interface(int interfaceIndex)	{return m_vInterfaces[interfaceIndex]->connectedProcID;}
		int get_interface_index_by_proc_id(ProcID procID);///<	returns -1 if no interface with the specified process exists.
		inline HNODE get_interface_node(int interfaceIndex, int nodeIndex)	{return m_vInterfaces[interfaceIndex]->vNodeHandles[nodeIndex];}
		inline int num_nodes_in_interface(int interfaceIndex)				{return m_vInterfaces[interfaceIndex]->vNodeHandles.size();}
		inline int num_normal_nodes_in_interface(int interfaceIndex)		{return m_vInterfaces[interfaceIndex]->numNodesOfType[NT_NORMAL];}
		inline int num_master_nodes_in_interface(int interfaceIndex)		{return m_vInterfaces[interfaceIndex]->numNodesOfType[NT_MASTER];}
		inline int num_slave_nodes_in_interface(int interfaceIndex)			{return m_vInterfaces[interfaceIndex]->numNodesOfType[NT_SLAVE];}

	////////////////////////////////
	//	data managment
	///	this method can be used to check if a buffer-handle points to a valid buffer.
		inline bool dbhandle_is_valid(const DBHANDLE& hDB)	{return hDB.bufferSetID != 0;}
		
	/**	This method creates a send data-buffer in all interfaces at once.
	  *	Data can be set using set_data. This data will be stored in all buffers in which the specified node is contained.
	  *	send-buffers are invalidated as soon as new nodes are added or old nodes are removed from any interface.
	  *	elementSize has to be specified in bytes. (use sizeof(...))*/
		inline DBHANDLE create_send_buffer(unsigned int elementSize, int dataID)
			{return create_interface_buffers(elementSize, dataID, true);}

	/**	This method deletes all send buffers. Should be called each time before a the buffers
	  *	for the next communication-step are set up - at least if something changed in the group.*/
		void delete_send_buffers();
		
	/**	This data will be stored in all interfaces in which the specified node is contained.
	  *	hDB has to be the handle of a send-buffer!*/
		bool set_data(const DBHANDLE& hDB, HNODE hNode, void* pData);
		
	/**	This method creates a receive data-buffer in all interfaces at once.
	  *	Data can be queried using get_data. This data will be stored in all buffers in which the specified node is contained.
	  *	receive-buffers are invalidated as soon as new nodes are added or old nodes are removed from any interface.
	  *	elementSize has to be specified in bytes. (use sizeof(...))*/
		inline DBHANDLE create_receive_buffer(unsigned int elementSize, int dataID)
			{return create_interface_buffers(elementSize, dataID, false);}

	/**	This method deletes all receive buffers. Should be called each time before the buffers
	  *	for the next communication-step are set up - at least if something changed in the group.*/
		void delete_receive_buffers();
		
	/**	You can query a buffer for data using this method.
	  * You have to specify the interface from which you the data should be taken.
	  *	This interface has to contain hNode!
	  *	This method is not the fastest.
	  *	Whenever possible you should use the alternate template-implementations of get_data
	  *	get_data_combined(...) and get_data_vec(...).*/
		//bool get_data(const DBHANDLE& hDB, int interfaceIndex, HNODE hNode, void* pDataOut);

	/**	This template-method allows to combine the data associated with a node from all
	  *	interfaces, that contain the node. You may specify an operation that defines
	  *	how the result is calculated from the different datas.
	  *	Be sure that the template parameter matches the original data-type of the data.
	  *	In order to average data, you can use the data operation DO_ADD and divide the
	  *	result by the number of interfaces in which hNode is contained
	  *	(use num_containing_interfaces)*/
		template <class TDataType>
		TDataType get_data(const DBHANDLE& hDB, HNODE hNode, DataOperation dop);

	/**	This template-method returns all data that is associated with a node in std::vector.
	  *	Be sure that the template parameter matches the original data-type of the data.*/
/*		template <class TDataType>
		bool get_data_vec(const DBHANDLE& hDB, HNODE hNode, std::vector<TDataType>& vDataOut);
*/

	/**	schedule a send buffer for communication if you want its data to be included in the next communication step.
	  *	Only data of scheduled buffers will be transfered.
	  *	buffers may be scheduled multiple times with varying srcNode-types.
	  * The order in which buffers are scheduled is important.
	  *	It has to match the order of the receive buffers.*/
		inline void schedule_send_buffer(const DBHANDLE& hDB, NodeType srcNodes)
			{schedule_buffer(hDB, true, srcNodes);}

	/**	Use this method if you want to transfer data from all boundary-nodes,
	  *	regardless of their type.*/
		inline void schedule_send_buffer(const DBHANDLE& hDB)
			{schedule_buffer(hDB, true, -1);}
		
	/**	schedule a receive buffer for communication if you want it to receive data in the next communication step.
	  *	buffers may be scheduled multiple times with varying srcNode-types.
	  * The order in which buffers are scheduled is important.
	  *	It has to match the order of the send buffers.*/
		inline void schedule_receive_buffer(const DBHANDLE& hDB, NodeType srcNodes)
			{schedule_buffer(hDB, false, srcNodes);}

	/**	Use this method if you want to receive data from all boundary-nodes,
	  *	regardless of their type.
	  *	Be sure to call the corresponding method when scheduling the corresponding send buffer.*/
		inline void schedule_receive_buffer(const DBHANDLE& hDB)
			{schedule_buffer(hDB, false, -1);}
		
	/**
	 * usig personal-data you may transmit a message from a node to all associated
	 * nodes on other processes.
	 * This method is suited if you want to transmit data only for a tiny fraction
	 * of the interface-nodes. Otherwise the use of a send-buffer is recommended.
	 * (\sa schedule_send_buffer).
	 * Note that the message will not be transmitted until \sa communicate() has
	 * been executed.*/
		void schedule_personal_data(HNODE hNode, const void* pData, int dataSize, int tag);
		
	/** More convenient than \sa schedule_personal_data - but has the same effect.*/
		template <class TData>
		inline void schedule_personal_data(HNODE hNode, const TData& data, int tag)
			{schedule_personal_data(hNode, &data, sizeof(TData), tag);}

	/**	This method performs the actual communication between the processes.
	  *	Data in scheduled send-buffers will be transfered and scheduled 
	  * receive-buffers will be filled.*/
		bool communicate();
		
	/**	returns true if personal data has been received.*/
		bool personal_data_received();
	/**	fills pDataOut and HNodeOut with the first received
	 * data-block for the given tag. Returns false if there is no such data.*/
		bool get_first_personal_data(HNODE& hNodeOut, void* pDataOut,
									int dataSize, int tag);

	/** More convenient than \sa get_first_personal_data - but has the same effect.*/
		template <class TData>
		inline bool get_first_personal_data(HNODE& hNodeOut, TData& dataOut, int tag)
			{return get_first_personal_data(hNodeOut, &dataOut, sizeof(TData), tag);}
			
	/** fills the HNodeOut, pDataOut with the next received
	 * personal data-block for the given tag. Call this method repeatedly
	 * until false is returned to receive all personal data of the given tag.
	 * Before the first call to this method you have to call
	 * \sa get_first_personal_data.*/ 
		bool get_next_personal_data(HNODE& hNodeOut, void* pDataOut,
									int dataSize, int tag);

	/** More convenient than \sa get_next_personal_data - but has the same effect.*/
		template <class TData>
		inline bool get_next_personal_data(HNODE& hNodeOut, TData& dataOut, int tag)
			{return get_next_personal_data(hNodeOut, &dataOut, sizeof(TData), tag);}
												
/*
		void clear_node_data();
	///	pData has to have as many entries as there are nodes.
		void add_node_data_field(void* pData, uint dataSize, DataTypeID dtID);
	///	pData has to have as many entries as there are nodes in the specified interface.
		void add_interface_data_field(int interfaceIndex, void* pData, uint dataSize, DataTypeID dtID);
	///	personal messages may be send to a single interface element on another processor
		void add_personal_data(int interfaceIndex, int nodeInterfaceID, void* pData, uint dataSize, DataTypeID dtID);
	///	send to all connected slaves or to the master-node
		void add_personal_data(const TNode& node, void* pData, uint dataSize, DataTypeID dtID);

		int num_node_data_fields();
		uint get_node_data_field_size(int nodeDataFieldIndex);
		DataTypeID get_node_data_field_type_id(int nodeDataFieldIndex);
		const void* get_node_data_field(int nodeDataFieldIndex);

		int num_interface_data_fields(int interfaceIndex);
		uint get_interface_data_field_size(int interfaceIndex, int interfaceDataFieldIndex);
		DataTypeID get_interface_data_field_type_id(int interfaceIndex, int interfaceDataFieldIndex);
		const void* get_interface_data_field(int interfaceIndex, int interfaceDataFieldIndex);

		int num_personal_datas(int interfaceIndex);
		uint get_personal_data_size(int interfaceIndex, int dataIndex);
		DataTypeID get_personal_data_type_id(int interfaceIndex, int dataIndex);
		const void* get_personal_data(int interfaceIndex, int dataIndex);


	////////////////////////////////
	//	communication packets
	//	use the following methods to send the whole group to another processor
		uint get_packet_size();
		bool pack(void** packHere);
		bool unpack(void* packet, uint packetSize);

	//	you may use this to communicate interface data
		uint get_interface_packet_size(int interfaceIndex);
		bool pack_interface_data(void** packHere, int interfaceIndex);
		bool unpack_interface_data(int interfaceIndex, void* packet, uint packetSize);
*/
	protected:
		class Interface;
		
		inline Interface* get_interface(int i)	{return m_vInterfaces[i];}
		
	/**	This method creates the data-buffers - either a send-buffer or a receive-buffer.*/
		DBHANDLE create_interface_buffers(unsigned int elementSize, int dataID, bool sendBuffers);
		
	/**	This method schedules buffers.*/
		void schedule_buffer(const DBHANDLE& hDB, bool sendBuffer, int nodeTypeID);
		
	///	returns the associated send buffer
		inline DataBuffer* get_send_buffer(const DBHANDLE& hDB, int interfaceIndex)
			{return m_vInterfaces[interfaceIndex]->vSendBuffers[hDB.bufferIndex];}

	///	returns the associated receive buffer
		inline DataBuffer* get_receive_buffer(const DBHANDLE& hDB, int interfaceIndex)
			{return m_vInterfaces[interfaceIndex]->vReceiveBuffers[hDB.bufferIndex];}
		
	protected:
	/** The NodeInfo stores information that is associated with a given node:
	  *	node:	The node itself.
	  *	nt:		the type of the node.
	  *	refIter:	an iterator to the entry of the NodeInfo in the groups section container.
  	  *	vInterfaceIndices:	indices to the interfaces in which the node is contained.
	  * vInterfaceRefIters:	indices to the entry at which the node is stored in the corresponding interface (use in combination with vInterfaceIndices).
	  */
		struct NodeInfo
		{
			friend class CommunicationGroup<TNode>;

			protected:
				NodeInfo()	{};
			//TODO: create protected copy-constructor.
				~NodeInfo()	{};
				
			protected:
				TNode				node;
				NodeType			nt;
				NodeHandleIterator	refIter;//	iterator to this NodeInfo-entry in the section container.
				std::vector<int>	vInterfaceIndices;///< indices to the interfaces in which the node is referenced.
				std::vector<int>	vInterfaceEntryIndices;///<	indices to the entry at which the node is stored in the corresponding interface (use in combination with vInterfaceIndices).
		};

	/**	A DataBuffer is used to store data that shall be send on the next communicate statement.
	  */
		struct DataBuffer
		{
			unsigned int		setID;	///< each set of data-buffers has a different id.
			unsigned int		elementSize;
			int					dataID;
			std::vector<byte>	vData;
		};
		
		typedef std::vector<DataBuffer*>	DataBufferVec;
		
	/**	The Interface struct is used internally to store the nodes that were assigned to an interface.
	  *	Since every interface represents the connection to an other process, the id of the connected
	  * process has to be stored as well.
	  */
		struct Interface
		{
			Interface();
						
			//void erase_node(HNODE hNode);

			ProcID			connectedProcID;
			NodeHandleVec	vNodeHandles;
			unsigned int	numNodesOfType[NT_NUM_TYPES];
			DataBufferVec	vSendBuffers;
			DataBufferVec	vReceiveBuffers;
			ug::BinaryStream	personalDataSendStream;
			ug::BinaryStream	personalDataReceiveStream;
			
		//TODO: don't use vectors for Send- and ReceiveData but byte* with new and delete in order to avoid unnecessary copying of data (minor performance-boost, since called seldom).
			std::vector<byte>	vSendData;///< this buffer is used to avoid reallocation in each communication step.
			std::vector<byte>	vReceiveData;///< this buffer is used to avoid reallocation in each communication step.
		};
		
		typedef std::pair<int, int>		PersonalDataNextPos; ///<first: interface-index, second: stream-pos.
		typedef std::vector<Interface*> InterfaceVec;
		typedef ug::Hash<int, ProcID>		InterfaceIndexHash;
		typedef ug::Hash<PersonalDataNextPos, int>	PersonalDataNextPosHash;
		
	/**	This struct is used to store, which buffers are scheduled.*/
		struct ScheduledBuffer
		{
			ScheduledBuffer()	{}
			ScheduledBuffer(const DBHANDLE& hDB, int nodeType) : dbHandle(hDB), nodeTypeID(nodeType)	{}
			
			DBHANDLE	dbHandle;
			int			nodeTypeID;//-1 for all nodes. Corresponds to NodeType else.
		};
		
		typedef std::vector<ScheduledBuffer>	ScheduledBufferVec;
		
	protected:
		ProcID				m_procID;
		NodeHandleContainer	m_nodeHandles;
		NodeHandleHash		m_hshNodeHandles;///< allows access to handles via node-object.
		InterfaceVec		m_vInterfaces;
		InterfaceIndexHash	m_hshInterfaceIndices;///< allows fast access to an interface given a procID.
	
		bool m_bNodeManagmentMode;
		bool m_bNodesErased;
	//	for communication.
		ScheduledBufferVec	m_vScheduledSendBuffers;
		ScheduledBufferVec	m_vScheduledReceiveBuffers;
		PersonalDataNextPosHash	m_hshPersonalDataNextPos;///< used for get_first_personal_data and get_next_personal_data.
};

}

//	include implementation
#include "pcl_impl.hpp"

#endif
