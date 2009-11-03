//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m06 d05

#ifndef __H__SDD_IMPL__
#define __H__SDD_IMPL__

#include <cassert>
#include <string.h>
#include "mpi.h"
#include "pcl.h"
#include "common.h"

namespace pcl
{

const uint DEFAULT_NODE_HASH_SIZE = 499;
const uint DEFAULT_INTERFACE_HASH_SIZE = 19;
const uint DEFAULT_PERSONAL_DATA_NEXT_HASH_SIZE = 19;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//	CommunicationGroup implementation
////////////////////////////////////////////////////////////////////////////////////////////////
template <class TNode>
CommunicationGroup<TNode>::
CommunicationGroup() : m_hshNodeHandles(DEFAULT_NODE_HASH_SIZE),
					  m_hshInterfaceIndices(DEFAULT_INTERFACE_HASH_SIZE),
					  m_hshPersonalDataNextPos(DEFAULT_PERSONAL_DATA_NEXT_HASH_SIZE)
{
	m_bNodeManagmentMode = false;
	m_bNodesErased = false;
}

template <class TNode>
CommunicationGroup<TNode>::
CommunicationGroup(uint nodeHashSize) : 
				m_hshNodeHandles(nodeHashSize),
				m_hshInterfaceIndices(DEFAULT_INTERFACE_HASH_SIZE),
				m_hshPersonalDataNextPos(DEFAULT_PERSONAL_DATA_NEXT_HASH_SIZE)
{
	m_bNodeManagmentMode = false;
	m_bNodesErased = false;

}

template <class TNode>
CommunicationGroup<TNode>::
CommunicationGroup(ProcID procID, uint nodeHashSize) : m_procID(procID), 
				m_hshNodeHandles(nodeHashSize),
				m_hshInterfaceIndices(DEFAULT_INTERFACE_HASH_SIZE),
				m_hshPersonalDataNextPos(DEFAULT_PERSONAL_DATA_NEXT_HASH_SIZE)
{
	m_bNodeManagmentMode = false;
	m_bNodesErased = false;
}

template <class TNode>
CommunicationGroup<TNode>::
CommunicationGroup(const CommunicationGroup<TNode>& dg) :
				m_hshNodeHandles(dg.m_hshNodeHandles.get_hash_size()),
				m_hshInterfaceIndices(DEFAULT_INTERFACE_HASH_SIZE),
				m_hshPersonalDataNextPos(DEFAULT_PERSONAL_DATA_NEXT_HASH_SIZE)
{
	m_bNodeManagmentMode = false;
	m_bNodesErased = false;
}

template <class TNode>
CommunicationGroup<TNode>::
~CommunicationGroup()
{
	clear();
}

template <class TNode>
void
CommunicationGroup<TNode>::
clear()
{
//	clear the hash
	m_hshNodeHandles.clear();
	
//	delete all nodes
	for(NodeHandleIterator iter = m_nodeHandles.begin();
		iter != m_nodeHandles.end(); ++iter)
	{
		delete *iter;
	}
	m_nodeHandles.clear();
	
//	clear send and receive buffers.
	delete_send_buffers();
	delete_receive_buffers();
	
//	delete all interfaces
	for(unsigned int i = 0; i < m_vInterfaces.size(); ++i)
	{
		delete m_vInterfaces[i];
	}
	m_vInterfaces.clear();
	
//	clear the interface hash
	m_hshInterfaceIndices.clear();
}

template <class TNode>
void
CommunicationGroup<TNode>::
set_node_hash_size(uint nodeHashSize)
{
	clear();
	m_hshNodeHandles.set_hash_size(nodeHashSize);
}

////////////////////////////////////////////////////////////////////////
//	node handling
template <class TNode>
void CommunicationGroup<TNode>::
begin_node_managment()
{
	m_bNodeManagmentMode = true;
	m_bNodesErased = false;
}

template <class TNode>
typename CommunicationGroup<TNode>::HNODE
CommunicationGroup<TNode>::add_node(const TNode& node)
{
//	create a new NodeInfo entry and insert it to m_nodeHandles.
	NodeInfo* pInfo = new NodeInfo();
	pInfo->node = node;
	pInfo->nt = NT_NORMAL;
	pInfo->refIter = m_nodeHandles.insert(pInfo, NT_NORMAL);
//	store the handle in the hash
	m_hshNodeHandles.add(pInfo, node);
	return pInfo;
}

template <class TNode>
void CommunicationGroup<TNode>::
erase_node(HNODE hNode)
{
//	erase the node from the handle-list and from the hash.
	m_nodeHandles.erase(hNode->refIter, hNode->nt);
	m_hshNodeHandles.erase(hNode->node);

//	if we are in node-managment-mode, we will mark the interface
//	node entry as invalid, but won't remove it yet.
	if(node_managment_is_active())
	{
	//	mark as invalid in interface-list.;
		for(uint i = 0; i < (uint)hNode->vInterfaceIndices.size(); ++i)
		{
			int interfaceInd = hNode->vInterfaceIndices[i];
			int entryInd = hNode->vInterfaceEntryIndices[i];
			get_interface(interfaceInd)->vNodeHandles[entryInd] = NULL;
		}
		m_bNodesErased = true;
	}
	else
	{
	//	since we are not in node-managment mode, we have to remove the node
	//	from all associated interfaces directly.
		const std::vector<int>& vInts = get_interface_indices(hNode);
		for(int i = 0; i < (uint)vInts.size(); ++i)
		{
			int nodeInterfaceIndex = vInts[i];
			Interface& interface = *get_interface(nodeInterfaceIndex);
		//	erase the node-entry
			interface.vNodeHandles.erase(
					interface.vNodeHandles.begin() + hNode->vInterfaceEntryIndices[i]);
		//	update the entry indices of all following nodes.
			for(int j = hNode->vInterfaceEntryIndices[i];
					j < interface.vNodeHandles.size(); ++j)
			{
				HNODE tNode = interface.vNodeHandles[j];
			//	update entry index
				for(int k = 0; k < tNode->vInterfaceIndices.size(); ++k)
				{
					if(tNode->vInterfaceIndices[k] == nodeInterfaceIndex)
					{
					//	update the corresponding entry index
						tNode->vInterfaceEntryIndices[k]--;
						break;
					}
				}
			}
		}
	}
	
	delete hNode;
}

template <class TNode>
typename CommunicationGroup<TNode>::HNODE
CommunicationGroup<TNode>::
replace_node(const TNode& nodeOld, const TNode& nodeNew)
{
//	get the handle of the old node
	HNODE hNode = get_handle(nodeOld);
	
	if(handle_is_valid(hNode))
	{
	//	remove it from the hash
		m_hshNodeHandles.erase(nodeOld);
	//	set the new node
		hNode->node = nodeNew;
	//	add the handle to the hash - this time with nodeNew as key
		m_hshNodeHandles.add(hNode, nodeNew);
	}
	
//	done
	return hNode;
}

template <class TNode>
void CommunicationGroup<TNode>::
set_node_type(HNODE hNode, NodeType nt)
{	
	NodeType oldType = hNode->nt;
	
//	with the node-type, the section in which a node is stored will also change.
//	remove it from the old one first and store it in the new one then.
	m_nodeHandles.erase(hNode->refIter, oldType);
	hNode->refIter = m_nodeHandles.insert(hNode, nt);
	
//	update interfaces too!
	for(uint i = 0; i < hNode->vInterfaceIndices.size(); ++i)
	{
		int ii = hNode->vInterfaceIndices[i];
		
		m_vInterfaces[ii]->numNodesOfType[oldType]--;
		m_vInterfaces[ii]->numNodesOfType[nt]++;
	}
	
//	store the new type
	hNode->nt = nt;
}

////////////////////////////////////////////////////////////////////////
//	interface access
////////////////////////////////
template <class TNode>
int
CommunicationGroup<TNode>::
add_node_to_interface(HNODE hNode, ProcID connectedProcID)
{
	return add_node_to_interface(hNode, connectedProcID, hNode->nt);
}

template <class TNode>
int
CommunicationGroup<TNode>::
add_node_to_interface(HNODE hNode, ProcID connectedProcID, NodeType nt)
{
//	update the node type
	if(hNode->nt != nt)
	{
	//	transitions from master to slave and from slave to master should not be done using this method!
		if(hNode->nt == NT_MASTER && nt == NT_SLAVE)
		{
			assert(!"Master node may not be transformed to slave node in this mehtod!");
		}

		if(hNode->nt == NT_SLAVE && nt == NT_MASTER)
		{
			assert(!"Master node may not be transformed to slave node in this mehtod!");
		}
		
	//	change the node-type (this method will automatically update all refIters)
		set_node_type(hNode, nt);
	}
	
//	get the interface index.
	int interfaceIndex = get_interface_index_by_proc_id(connectedProcID);
	Interface* pInt = NULL;
	
//	check whether we have to create a new interface
	if(interfaceIndex == -1)
	{
	//	yes we do.
		pInt = new Interface();
		pInt->connectedProcID = connectedProcID;
	//	push the interface to the interface vector
		m_vInterfaces.push_back(pInt);
	//	add the interface index to the hash
		interfaceIndex = m_vInterfaces.size() - 1;
		m_hshInterfaceIndices.add(interfaceIndex, connectedProcID);
	}
	else
	{
	//	the interface already exists. Store the pointer.
		pInt = m_vInterfaces[interfaceIndex];
	}

//	add the node to the interface
	unsigned int entryIndex = pInt->vNodeHandles.size();
	pInt->vNodeHandles.push_back(hNode);
	hNode->vInterfaceIndices.push_back(interfaceIndex);
	hNode->vInterfaceEntryIndices.push_back(entryIndex);
	pInt->numNodesOfType[nt]++;

	return entryIndex;
}

template <class TNode>
void CommunicationGroup<TNode>::
end_node_managment()
{
	if(!m_bNodeManagmentMode)
		return;
		
	m_bNodeManagmentMode = false;
	
	if(m_bNodesErased)
	{
		m_bNodesErased = false;
	//	iterate through all interfaces and make sure that their
	//	node-lists only contain valid handles.
	//	update nodes interface-entry-indices if required.
		for(int i = 0; i < num_interfaces(); ++i)
		{
			Interface& interface = *get_interface(i);
			int insertTo = 0;
			for(int j = 0; j < interface.vNodeHandles.size(); ++j)
			{
				HNODE hNode = interface.vNodeHandles[j];
				if(hNode)
				{
					if(insertTo != j)
					{
					//	copy the node to close the gap.
						interface.vNodeHandles[insertTo] = interface.vNodeHandles[j];
					//	update entry index
						for(int k = 0; k < hNode->vInterfaceIndices.size(); ++k)
						{
							if(hNode->vInterfaceIndices[k] == i)
							{
							//	update the corresponding entry index
								hNode->vInterfaceEntryIndices[k] = insertTo;
								break;
							}
						}
					}
					insertTo++;
				}
			}
		//	the interface and associated entry-indices in associated nodes
		//	have been updated.
		//	resize the interface-vec to the new size.
			interface.vNodeHandles.resize(insertTo);
		}
	}

//TODO: update interface-buffers.
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template <class TNode>
typename CommunicationGroup<TNode>::HNODE
CommunicationGroup<TNode>::get_handle(const TNode& node)
{
	if(m_hshNodeHandles.has_entries(node))
		return m_hshNodeHandles.first(node);
	return invalid_handle();
}

template <class TNode>
int
CommunicationGroup<TNode>::
get_interface_index_by_proc_id(ProcID procID)
{
	if(m_hshInterfaceIndices.has_entries(procID))
		return m_hshInterfaceIndices.first(procID);
	return -1;
}

//protected
template <class TNode>
typename CommunicationGroup<TNode>::DBHANDLE
CommunicationGroup<TNode>::
create_interface_buffers(unsigned int elementSize, int dataID, bool sendBuffers)
{
//	this static variable defines the buffers id - will be increased with each set of new buffers.
	static unsigned int setID = 1;
	
	DBHANDLE hDB;
	
//	check if the arguments and the set-up is valid.
	if(num_interfaces() == 0)
	{
		hDB.bufferSetID = 0;
		return hDB;
	}
	if(elementSize == 0)
	{
		hDB.bufferSetID = 0;
		return hDB;
	}
	
//	create a new buffer in each interface
	for(int i = 0; i < num_interfaces(); ++i)
	{
		DataBuffer* pDB = new DataBuffer;
		pDB->setID = setID;
		pDB->dataID = dataID;
		pDB->elementSize = elementSize;
		pDB->vData.resize(elementSize * num_nodes_in_interface(i));//	send-buffers are always created with full size, since this makes indexing a lot easier.
		
	//	store the buffer in the interface
		if(sendBuffers)
			m_vInterfaces[i]->vSendBuffers.push_back(pDB);
		else
			m_vInterfaces[i]->vReceiveBuffers.push_back(pDB);
	}
	
//	set up the handle
	hDB.bufferSetID = setID;
	hDB.isInterfaceBuffer = true;
	hDB.isSendBuffer = sendBuffers;
	if(sendBuffers)
		hDB.bufferIndex = m_vInterfaces[0]->vSendBuffers.size() - 1;
	else
		hDB.bufferIndex = m_vInterfaces[0]->vReceiveBuffers.size() - 1;
			
	setID++;
	
	return hDB;
}

template <class TNode>
void
CommunicationGroup<TNode>::
delete_send_buffers()
{
//	iterate through all interfaces
	for(int i = 0; i < num_interfaces(); ++i)
	{
	//	iterate through all send-buffers and delete each.
		for(uint j = 0; j < m_vInterfaces[i]->vSendBuffers.size(); ++j)
			delete m_vInterfaces[i]->vSendBuffers[j];
			
	//	clear the buffer-list
		m_vInterfaces[i]->vSendBuffers.clear();
	}
}

template <class TNode>
void
CommunicationGroup<TNode>::
delete_receive_buffers()
{
//	iterate through all interfaces
	for(int i = 0; i < num_interfaces(); ++i)
	{
	//	iterate through all send-buffers and delete each.
		for(uint j = 0; j < m_vInterfaces[i]->vReceiveBuffers.size(); ++j)
			delete m_vInterfaces[i]->vReceiveBuffers[j];
			
	//	clear the buffer-list
		m_vInterfaces[i]->vReceiveBuffers.clear();
	}
}

template <class TNode>
bool
CommunicationGroup<TNode>::
set_data(const DBHANDLE& hDB, HNODE hNode, void* pData)
{
//	the handle has to be valid
	assert(dbhandle_is_valid(hDB) && "invalid handle!");

//	data of send-buffers may be set only!
	if(!hDB.isSendBuffer)
		return false;
		
//	we have to distinguish between interface- and non-interface-buffers.
	if(hDB.isInterfaceBuffer)
	{
	//	iterate through the interfaces that are connected with the node.
		for(int ii = 0; ii < hNode->vInterfaceIndices.size(); ++ii)
		{
			int interfaceIndex = hNode->vInterfaceIndices[ii];
		//	get the buffer
			DataBuffer* pDB = get_send_buffer(hDB, interfaceIndex);

		//	assert that the buffers setID and the setID of the handle correspond
			assert(pDB->setID == hDB.bufferSetID && "setIDs of buffer and handle do not correspond. Old handle?");
			
		//	copy the data
			int destDataIndex = hNode->vInterfaceEntryIndices[ii] * pDB->elementSize;
			memcpy(&pDB->vData[destDataIndex], pData, pDB->elementSize);
		}
	}
	else
	{
	//TODO:	add support for non-interface buffers.
		return false;//only interface buffers are supported in the moment.
	}

//	done
	return true;
}
/*
template <class TNode>
bool
CommunicationGroup<TNode>::
get_data(const DBHANDLE& hDB, int interfaceIndex, HNODE hNode, void* pDataOut)
{
//	the handle has to be valid
	assert(handle_is_valid(hDB) && "invalid handle!");

//	data of receive-buffers may be queried only!
	if(hDB.isSendBuffer)
		return false;
		
//	we have to distinguish between interface- and non-interface-buffers.
	if(hDB.isInterfaceBuffer)
	{
		{
		//	get the buffer
			DataBuffer* pDB = get_receive_buffer(hDB, interfaceIndex);
		
		//	assert that the buffers setID and the setID of the handle correspond
			assert(pDB->setID == hDB.bufferSetID && "setIDs of buffer and handle do not correspond. Old handle?");

		//	get the entry index
			int ii = -1;
			for(int i = 0; i < hNode->vInterfaceIndices.size(); ++i)
			{
				if(hNode->vInterfaceIndices[i] == interfaceIndex)
				{
					ii = i;
					break;
				}
			}
			
			if(ii < 0)
				return false;
							
		//	copy the data
			int srcDataIndex = hNode->vInterfaceEntryIndices[ii] * pDB->elementSize;
			memcpy(pDataOut, &pDB->vData[srcDataIndex], pDB->elementSize);
		}
	}
	else
	{
	//TODO:	add support for non-interface buffers.
		return false;//only interface buffers are supported in the moment.
	}

//	done
	return true;
}
*/
template <class TNode>
void
CommunicationGroup<TNode>::
schedule_buffer(const DBHANDLE& hDB, bool sendBuffer, int nodeTypeID)
{
//	the handle has to be valid
	assert(dbhandle_is_valid(hDB) && "invalid handle!");
	
	if(sendBuffer)
		m_vScheduledSendBuffers.push_back(ScheduledBuffer(hDB, nodeTypeID));
	else
		m_vScheduledReceiveBuffers.push_back(ScheduledBuffer(hDB, nodeTypeID));
}

template <class TNode>
void
CommunicationGroup<TNode>::
schedule_personal_data(HNODE hNode, const void* pData, int dataSize, int tag)
{
//	pack the data to the interfaces personal-data-stream.
//	iterate through the interfaces
	const std::vector<int>& interfaceInds = get_interface_indices(hNode);

	for(int i = 0; i < interfaceInds.size(); ++i)
	{
	//	get the data-stream
		Interface* pInt = get_interface(interfaceInds[i]);
		ug::BinaryStream& dataStream = pInt->personalDataSendStream;
	//	get the index of the node in the interface
		int nodeInterfaceIndex = hNode->vInterfaceEntryIndices[i];

	//	write personal data to the stream
	//	format: tag, nodeInterfaceIndex, pData
		dataStream.write((char*)&tag, sizeof(int));
		dataStream.write((char*)&nodeInterfaceIndex, sizeof(int));
		dataStream.write((const char*)pData, dataSize);
	}
}

template <class TNode>
bool
CommunicationGroup<TNode>::
communicate()
{
	if(num_interfaces() == 0)
		return false;
	
////////////////////////////////////////////////
//	set up send buffers

//	collect the send-data for each interface in a single big buffer.
//	iterate through all scheduled send-buffers and calculate the
//	buffer sizes for each interface.
//	each buffer has an additional int value that defines the size of the personal data that
//	will be send in an extra send-receive-step.
	std::vector<unsigned int> vSendBufferSizes(num_interfaces(), sizeof(int));
	
	for(unsigned int i = 0; i < m_vScheduledSendBuffers.size(); ++i)
	{
		const DBHANDLE& hDB = m_vScheduledSendBuffers[i].dbHandle;
		int nt = m_vScheduledSendBuffers[i].nodeTypeID;
		
	//	get the element-size from a representative buffer of the handles buffer-set
	//	(the element size is the same for all buffers in all interfaces that belong to one handle.
		int elementSize = get_send_buffer(hDB, 0)->elementSize;
		
	//	iterate through all interfaces and sum the sizes of the buffers.
		for(int j = 0; j < num_interfaces(); ++j)
		{
			switch(nt)
			{
			case -1:		vSendBufferSizes[j] += num_nodes_in_interface(j) * elementSize; break;
			case NT_NORMAL:	vSendBufferSizes[j] += num_normal_nodes_in_interface(j) * elementSize; break;
			case NT_MASTER:	vSendBufferSizes[j] += num_master_nodes_in_interface(j) * elementSize; break;
			case NT_SLAVE:	vSendBufferSizes[j] += num_slave_nodes_in_interface(j) * elementSize; break;
			default:
				assert(!"bad node type id in scheduled send-buffer!");
				break;
			}
		}
	}
	
//	adjust the buffers
//	the buffers are stored in the respective interfaces.
//	in each communication-step one has to check whether they are big enough.
//	if not they will be resized.
//	for performance-reasons they won't be shrinked again.
	std::vector<byte*> vSendBuffers(num_interfaces());
	
	for(int i = 0; i < num_interfaces(); ++i)
	{
		Interface* pI = get_interface(i);
		if(vSendBufferSizes[i] > 0)
		{
			if(vSendBufferSizes[i] > pI->vSendData.size())
				pI->vSendData.resize(vSendBufferSizes[i]);
			vSendBuffers[i] = &(pI->vSendData.front());
		}
		else
			vSendBuffers[i] = NULL;
	}
	

////////////////////////////////////////////////
//	fill send buffers

//	copy the data from scheduled send buffers to the combined send buffers.
	std::vector<unsigned int> vSendBufferFillStatus(num_interfaces(), 0);
	
	for(unsigned int i = 0; i < m_vScheduledSendBuffers.size(); ++i)
	{
		const DBHANDLE& hDB = m_vScheduledSendBuffers[i].dbHandle;
		int nt = m_vScheduledSendBuffers[i].nodeTypeID;
				
	//	iterate through all interfaces and collect the data
		for(int j = 0; j < num_interfaces(); ++j)
		{
		//	the buffer from which we will take the data.
			DataBuffer* pDB = get_send_buffer(hDB, j);
			
		//	we have to perform differnt operations - depending on the node-type.
			if(nt == -1)
			{
			//	copy the whole buffer
				void* pDest = vSendBuffers[j] + vSendBufferFillStatus[j];
				void* pSrc = &(pDB->vData[0]);
				memcpy(pDest, pSrc, pDB->vData.size());
				vSendBufferFillStatus[j] += pDB->vData.size();
			}
			else
			{
			//	check whether there are nodes of this type.
			//	If the buffer consists only of nodes of this type, then the buffer may be copied at once.
				Interface* pI = m_vInterfaces[j];
				unsigned int numNodes = pI->numNodesOfType[nt];
				if(numNodes > 0)
				{
					if(numNodes == num_nodes_in_interface(j))
					{
					//	the buffer may be copied at once
						void* pDest = vSendBuffers[j] + vSendBufferFillStatus[j];
						void* pSrc = &(pDB->vData[0]);
						memcpy(pDest, pSrc, pDB->vData.size());
						vSendBufferFillStatus[j] += pDB->vData.size();
					}
					else
					{
					//	iterate over all nodes in the interface.
					//	if the type matches the scheduled-buffers node type, then copy its data.
						for(int k = 0; k < num_nodes_in_interface(j); ++k)
						{
							HNODE hNode = get_interface_node(j, k);
							if(get_node_type(hNode) == nt)
							{
							//	copy the data
								void* pDest = vSendBuffers[j] + vSendBufferFillStatus[j];
								void* pSrc = &(pDB->vData[k*pDB->elementSize]);
//PLOG(m_procID, "copying to send-buffer (I: " << j << "N: " << k << "): " << (*(float*)pSrc) << "\n");
								memcpy(pDest, pSrc, pDB->elementSize);
								vSendBufferFillStatus[j] += pDB->elementSize;
							}
						}
					}
				}
			}
		}
	}
//	iterate over the interfaces again and add the personal-data-size for each
//	interface to the send-buffer
	for(int i = 0; i < num_interfaces(); ++i)
	{
		int* pDest = (int*)(vSendBuffers[i] + vSendBufferFillStatus[i]);
		*pDest = get_interface(i)->personalDataSendStream.size();
		vSendBufferFillStatus[i] += sizeof(int);
	}
	
//	the send buffers are now complete.


////////////////////////////////////////////////
//	set up receive buffers

//	each buffer holds an additional integer that is used to receive the size of
//	the personal-data-stream that is received in the next communication-step.
	std::vector<unsigned int> vReceiveBufferSizes(num_interfaces(), sizeof(int));
	
	for(unsigned int i = 0; i < m_vScheduledReceiveBuffers.size(); ++i)
	{
		const DBHANDLE& hDB = m_vScheduledReceiveBuffers[i].dbHandle;
		int nt = m_vScheduledReceiveBuffers[i].nodeTypeID;
		
	//	get the element-size from a representative buffer of the handles buffer-set
	//	(the element size is the same for all buffers in all interfaces that belong to one handle.
		int elementSize = get_receive_buffer(hDB, 0)->elementSize;
		
	//	iterate through all interfaces and sum the sizes of the buffers.
		for(int j = 0; j < num_interfaces(); ++j)
		{
			switch(nt)
			{
			case -1:		vReceiveBufferSizes[j] += num_nodes_in_interface(j) * elementSize; break;
			case NT_NORMAL:	vReceiveBufferSizes[j] += num_normal_nodes_in_interface(j) * elementSize; break;
			case NT_MASTER:	vReceiveBufferSizes[j] += num_master_nodes_in_interface(j) * elementSize; break;
			case NT_SLAVE:	vReceiveBufferSizes[j] += num_slave_nodes_in_interface(j) * elementSize; break;
			default:
				assert(!"bad node type id in scheduled receive-buffer!");
				break;
			}
		}
	}
	
//	adjust the buffers
//	the buffers are stored in the respective interfaces.
//	in each communication-step one has to check whether they are big enough.
//	if not they will be resized.
//	for performance-reasons they won't be shrinked again.
	std::vector<byte*> vReceiveBuffers(num_interfaces());
	
	for(int i = 0; i < num_interfaces(); ++i)
	{
		Interface* pI = get_interface(i);
		if(vReceiveBufferSizes[i] > 0)
		{
			if(vReceiveBufferSizes[i] > pI->vReceiveData.size())
				pI->vReceiveData.resize(vReceiveBufferSizes[i]);
			vReceiveBuffers[i] = &(pI->vReceiveData.front());
		}
		else
			vReceiveBuffers[i] = NULL;
	}

// everything is ready for communication now


////////////////////////////////////////////////
//	perform communication
	int tag = 189345;//	an arbitrary number
	
//	first send all send-buffers (non-blocking)
	std::vector<MPI_Request> vSendRequests(num_interfaces());
	
	for(int i = 0; i < num_interfaces(); ++i)
	{
		if(vSendBufferSizes[i] > 0)
		{
			int retVal = MPI_Isend(vSendBuffers[i], (int)vSendBufferSizes[i], MPI_UNSIGNED_CHAR,
					get_connected_proc_of_interface(i), tag, MPI_COMM_WORLD, &vSendRequests[i]);
		}
		else
			vSendRequests[i] = MPI_REQUEST_NULL;
	}
	
//	receive data
	std::vector<MPI_Request> vReceiveRequests(num_interfaces());
	
	for(int i = 0; i < num_interfaces(); ++i)
	{
		if(vReceiveBufferSizes[i] > 0)
			MPI_Irecv(vReceiveBuffers[i], vReceiveBufferSizes[i], MPI_UNSIGNED_CHAR,	
					get_connected_proc_of_interface(i), tag, MPI_COMM_WORLD, &vReceiveRequests[i]);
		else
			vReceiveRequests[i] = MPI_REQUEST_NULL;
	}
	
//	wait until data has been received
	std::vector<MPI_Status> vReceiveStatii(num_interfaces());//TODO: fix spelling!
	
//	TODO: this can be improved:
//		instead of waiting for all, one could wait until one has finished and directly
//		start copying the data to the local receive buffer. Afterwards on could continue
//		by waiting for the next one etc...
	MPI_Waitall(num_interfaces(), &vReceiveRequests[0], &vReceiveStatii[0]);
	
////////////////////////////////////////////////
//	copy the data from the received buffers to the scheduled-receive-buffers.
	std::vector<unsigned int> vReceiveBufferPositions(num_interfaces(), 0);
	
	for(unsigned int i = 0; i < m_vScheduledReceiveBuffers.size(); ++i)
	{
		const DBHANDLE& hDB = m_vScheduledReceiveBuffers[i].dbHandle;
		int nt = m_vScheduledReceiveBuffers[i].nodeTypeID;
				
	//	iterate through all interfaces
		for(int j = 0; j < num_interfaces(); ++j)
		{
		//	the buffer to which we will copy data
			DataBuffer* pDB = get_receive_buffer(hDB, j);
			
		//	we have to perform differnt operations - depending on the node-type.
			if(nt == -1)
			{
			//	copy to the whole buffer
				void* pDest = &(pDB->vData[0]);
				void* pSrc = vReceiveBuffers[j] + vReceiveBufferPositions[j];
				memcpy(pDest, pSrc, pDB->vData.size());
				vReceiveBufferPositions[j] += pDB->vData.size();
			}
			else
			{
			//	check whether there are nodes of this type.
			//	If the buffer consists only of nodes of this type, then the buffer may be copied at once.
				Interface* pI = m_vInterfaces[j];
				unsigned int numNodes = pI->numNodesOfType[nt];
				if(numNodes > 0)
				{
					if(numNodes == num_nodes_in_interface(j))
					{
					//	the buffer may be copied at once
						void* pDest = &(pDB->vData[0]);
						void* pSrc = vReceiveBuffers[j] + vReceiveBufferPositions[j];
						memcpy(pDest, pSrc, pDB->vData.size());
						vReceiveBufferPositions[j] += pDB->vData.size();
					}
					else
					{
					//	iterate over all nodes in the interface.
					//	if the type matches the scheduled-buffers node type, then copy its data.
						for(int k = 0; k < num_nodes_in_interface(j); ++k)
						{
							HNODE hNode = get_interface_node(j, k);
							if(get_node_type(hNode) == nt)
							{
							//	copy the data
								void* pDest = &(pDB->vData[k*pDB->elementSize]);
								void* pSrc = vReceiveBuffers[j] + vReceiveBufferPositions[j];
								memcpy(pDest, pSrc, pDB->elementSize);
								vReceiveBufferPositions[j] += pDB->elementSize;
							}
						}
					}
				}
			}
		}
	}

//	resize and reset the personalDataReceiveStream
	for(int i = 0; i < num_interfaces(); ++i)
	{
		int newSize = *(int*)(vReceiveBuffers[i] + vReceiveBufferPositions[i]);
		get_interface(i)->personalDataReceiveStream.resize(newSize);
		get_interface(i)->personalDataReceiveStream.reset();
	}
	
//	 wait until all sends are finished
	std::vector<MPI_Status> vSendStatii(num_interfaces());//TODO: fix spelling!
	MPI_Waitall(num_interfaces(), &vSendRequests[0], &vSendStatii[0]);
		
//	clear the schedule-list
	m_vScheduledSendBuffers.clear();
	m_vScheduledReceiveBuffers.clear();
	
////////////////////////////////////////////////
//	now start personal-data-communication.
//	first send personalDataSendStream
	for(int i = 0; i < num_interfaces(); ++i)
	{
		ug::BinaryStream& binStream = get_interface(i)->personalDataSendStream;
		if(binStream.size() > 0)
		{
			int retVal = MPI_Isend(binStream.buffer(), binStream.size(), MPI_UNSIGNED_CHAR,
					get_connected_proc_of_interface(i), tag, MPI_COMM_WORLD, &vSendRequests[i]);
		}
		else
			vSendRequests[i] = MPI_REQUEST_NULL;
	}
	
//	receive data
	for(int i = 0; i < num_interfaces(); ++i)
	{
		ug::BinaryStream& binStream = get_interface(i)->personalDataReceiveStream;
		if(binStream.size() > 0)
		{
			MPI_Irecv(binStream.buffer(), binStream.size(), MPI_UNSIGNED_CHAR,	
					get_connected_proc_of_interface(i), tag, MPI_COMM_WORLD, &vReceiveRequests[i]);
		}
		else
			vReceiveRequests[i] = MPI_REQUEST_NULL;
	}
	
//	wait until data has been received
	MPI_Waitall(num_interfaces(), &vReceiveRequests[0], &vReceiveStatii[0]);

//	wait until all sends are finished
	MPI_Waitall(num_interfaces(), &vSendRequests[0], &vSendStatii[0]);

//	clear the personalDataSendStream
	for(int i = 0; i < num_interfaces(); ++i)
		get_interface(i)->personalDataSendStream.clear();
	
//	done. FINALLY!
	return true;
}

template <class TNode>
bool
CommunicationGroup<TNode>::
personal_data_received()
{
//	check if any personalDataReceiveStream contains data
	for(int i = 0; i < num_interfaces(); ++i)
	{
		if(get_interface(i)->personalDataReceiveStream.size() > 0)
			return true;
	}
//	none contains any data
	return false;
}

template <class TNode>
bool
CommunicationGroup<TNode>::
get_first_personal_data(HNODE& hNodeOut, void* pDataOut, int dataSize, int tag)
{
	if(!m_hshPersonalDataNextPos.has_entries(tag))
		m_hshPersonalDataNextPos.add(PersonalDataNextPos(0, 0), tag);
	else
		m_hshPersonalDataNextPos.first(tag) = PersonalDataNextPos(0, 0);
	
	return get_next_personal_data(hNodeOut, pDataOut, dataSize, tag);
}
						
template <class TNode>
bool
CommunicationGroup<TNode>::
get_next_personal_data(HNODE& hNodeOut, void* pDataOut, int dataSize, int tag)
{
	assert(m_hshPersonalDataNextPos.has_entries(tag) &&
			"call get_first_personal_data before calling this method.");

	PersonalDataNextPos& nextPos = m_hshPersonalDataNextPos.first(tag);
	
//	iterate through the interfaces until the next valid entry is found.
	while(nextPos.first < num_interfaces())	
	{
		Interface* pInterface = get_interface(nextPos.first);
	//	get the stream and the buffer of the active interface
		ug::BinaryStream& stream = pInterface->personalDataReceiveStream;
		char* pBuff = stream.buffer() + nextPos.second;

	//	iterate through the stream and find the next entry with the matching tag
		while(nextPos.second < stream.size())
		{
		//	format: {tag, nodeInterfaceIndex, pData}
		//	get the tag
			int streamTag = *(int*)pBuff;
			pBuff += sizeof(int);
		//	if the tag matches, then fill the data
			if(streamTag == tag)
			{
				int nodeInterfaceIndex = *(int*)pBuff;
				pBuff += sizeof(int);
				memcpy(pDataOut, pBuff, dataSize);
			//	assign hNode
				hNodeOut = pInterface->vNodeHandles[nodeInterfaceIndex];
			//	increase nextPos
				nextPos.second += 2*sizeof(int)+dataSize;
			//	done
				return true;
			} 
		//	increase nextPos
			nextPos.second += 2*sizeof(int)+dataSize;
		}
		nextPos.second = 0;
		nextPos.first++;
	}
	return false;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	interface implementation
////////////////////////////////////////////////////////////////////////
template <class TNode>
CommunicationGroup<TNode>::
Interface::
Interface()
{
//	initialize all numbers
	for(unsigned int i = 0; i < NT_NUM_TYPES; ++i)
		numNodesOfType[i] = 0;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	template class-methods
////////////////////////////////////////////////////////////////////////
template <class TNode>
template <class TDataType>
TDataType
CommunicationGroup<TNode>::
get_data(const DBHANDLE& hDB, HNODE hNode, DataOperation dop)
{
//	the handle has to be valid
	assert(dbhandle_is_valid(hDB) && "invalid handle!");

//	hNode has to be an interface node.
	assert((hNode->vInterfaceIndices.size() > 0) && "node is not an interface node!");
	
//	data of receive-buffers may be queried only!
	assert(!hDB.isSendBuffer && "get_data may only be performed on receive buffers!");
		
//	we have to distinguish between interface- and non-interface-buffers.
	if(hDB.isInterfaceBuffer)
	{			
		switch(dop)
		{
			case DO_FIRST:
			{
				DataBuffer* pDB = get_receive_buffer(hDB, hNode->vInterfaceIndices[0]);
				return *((TDataType*)&pDB->vData[hNode->vInterfaceEntryIndices[0] * pDB->elementSize]);
			}
			case DO_ADD:
			{
			//	sum the datas that belong to this node from the different interfaces.
				TDataType retVal = 0;
				for(int i = 0; i < hNode->vInterfaceIndices.size(); ++i)
				{
				//	get the buffer
					DataBuffer* pDB = get_receive_buffer(hDB, hNode->vInterfaceIndices[i]);
				//	sum it up
					retVal += *((TDataType*)&pDB->vData[hNode->vInterfaceEntryIndices[i] * pDB->elementSize]);
				}
				return retVal;
			}
			
			case DO_MIN:
			{
			//	find the minimum
				DataBuffer* pDB = get_receive_buffer(hDB, hNode->vInterfaceIndices[0]);
				TDataType retVal = *((TDataType*)&pDB->vData[hNode->vInterfaceEntryIndices[0] * pDB->elementSize]);
				for(int i = 1; i < hNode->vInterfaceIndices.size(); ++i)
				{
				//	get the buffer
					DataBuffer* pDB = get_receive_buffer(hDB, hNode->vInterfaceIndices[i]);
					retVal = min(retVal, *((TDataType*)&pDB->vData[hNode->vInterfaceEntryIndices[i] * pDB->elementSize]));
				}
				return retVal;
			}

			case DO_MAX:
			{
			//	find the maximum
				DataBuffer* pDB = get_receive_buffer(hDB, hNode->vInterfaceIndices[0]);
				TDataType retVal = *((TDataType*)&pDB->vData[hNode->vInterfaceEntryIndices[0] * pDB->elementSize]);
				for(int i = 1; i < hNode->vInterfaceIndices.size(); ++i)
				{
				//	get the buffer
					DataBuffer* pDB = get_receive_buffer(hDB, hNode->vInterfaceIndices[i]);
					retVal = max(retVal, *((TDataType*)&pDB->vData[hNode->vInterfaceEntryIndices[i] * pDB->elementSize]));
				}
				return retVal;
			}
		}
	}
	
	assert(!"program shouldn't reach this point. Forgot to implement some data-operations?");
	return TDataType();
}

}//	end of namespace

#endif
