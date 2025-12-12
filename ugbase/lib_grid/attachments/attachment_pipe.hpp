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

#ifndef __UTIL__ATTACHMENT_PIPE__IMPL__
#define __UTIL__ATTACHMENT_PIPE__IMPL__

#include "attachment_pipe.h"
#include "common/profiler/profiler.h"

namespace ug {

template <typename TElem, typename TElemHandler>
AttachmentPipe<TElem, TElemHandler>::
AttachmentPipe() :
	m_attachmentEntryIteratorHash(13),
	m_numElements(0),
	m_numDataEntries(0),
	m_containerSize(0),
	m_pHandler(nullptr)
{
}

template <typename TElem, typename TElemHandler>
AttachmentPipe<TElem, TElemHandler>::
AttachmentPipe(typename atraits::ElemHandlerPtr pHandler) :
	m_attachmentEntryIteratorHash(13),
	m_numElements(0),
	m_numDataEntries(0),
	m_containerSize(0),
	m_pHandler(pHandler)
{
}

template <typename TElem, typename TElemHandler>
AttachmentPipe<TElem, TElemHandler>::
~AttachmentPipe()
{
//	clear attachment entries
	for(auto iter = m_attachmentEntryContainer.begin();
	    iter != m_attachmentEntryContainer.end(); ++iter)
	{
		delete (*iter).m_pAttachment;
		delete (*iter).m_pContainer;
	}
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
clear()
{
	clear_elements();
	clear_attachments();
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
clear_elements()
{
	m_stackFreeEntries = UINTStack();
	resize_attachment_containers(0);
	m_numElements = 0;
	m_numDataEntries = 0;
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
clear_attachments()
{
	m_attachmentEntryContainer.clear();
	m_attachmentEntryIteratorHash.clear();
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
reset_values(size_t dataIndex)
{
	for(auto iter = m_attachmentEntryContainer.begin();
	    iter != m_attachmentEntryContainer.end(); ++iter)
	{
		(*iter).m_pContainer->reset_entry(dataIndex);
	}
}



template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
reserve(size_t numElems)
{
	if(numElems > get_container_size())
		resize_attachment_containers(numElems);
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
register_element(TElem elem)
{
	size_t newInd;
//	check if there is a free entry

	if(!m_stackFreeEntries.empty())
	{
	//	take it
		newInd = m_stackFreeEntries.top();
		m_stackFreeEntries.pop();

	//	initialize attachments with default constructor
		reset_values(newInd);
	}
	else
	{
	//	add a new entry
		newInd = m_numDataEntries++;
	//	make sure that the data containers are big enough.
		grow_attachment_containers(m_numDataEntries);
	}

//	assign new attachment index to the element
	atraits::set_data_index(m_pHandler, elem, newInd);

//	increase element count
	m_numElements++;
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
unregister_element(const TElem& elem)
{
//	get the attachment index
	size_t ind = atraits::get_data_index(m_pHandler, elem);
//	store the index in the stack of free entries
	m_stackFreeEntries.push(ind);
//	decrease element count
	m_numElements--;
}


template <typename TElem, typename TElemHandler>
template <typename TAttachment>
void
AttachmentPipe<TElem, TElemHandler>::
attach(TAttachment& attachment,
		const typename TAttachment::ValueType& defaultValue,
		uint options)
{
//	make sure, that the attachment is not already attached
	if(!has_attachment(attachment))
	{
	//	create a new element and insert it into the entryContainer.
	//	resize the new data-container to the size of the element container.
	//	store the returned iterator in the iteratorHash to allow fast access.
		IAttachment* pClonedAttachment = attachment.clone();
		IAttachmentDataContainer* pClonedContainer = attachment.create_container(defaultValue);
		pClonedContainer->resize(get_container_size());
		auto iter = m_attachmentEntryContainer.insert(
				m_attachmentEntryContainer.end(),
						AttachmentEntry(pClonedAttachment,
						pClonedContainer, options));
		m_attachmentEntryIteratorHash.insert(pClonedAttachment->id(), iter);
	}
	else
	{
//TODO: do something!
	}
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
attach(IAttachment& attachment, uint options)
{
//	make sure, that the attachment is not already attached
	if(!has_attachment(attachment))
	{
	//	create a new element and insert it into the entryContainer.
	//	resize the new data-container to the size of the element container.
	//	store the returned iterator in the iteratorHash to allow fast access.
		IAttachment* pClonedAttachment = attachment.clone();
		IAttachmentDataContainer* pClonedContainer = pClonedAttachment->create_container();
		pClonedContainer->resize(get_container_size());
		auto iter = m_attachmentEntryContainer.insert(m_attachmentEntryContainer.end(), AttachmentEntry(pClonedAttachment, pClonedContainer, options));
		m_attachmentEntryIteratorHash.insert(pClonedAttachment->id(), iter);
	}
	else
	{
//TODO: do something!
	}
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
detach(IAttachment& attachment)
{
	if(has_attachment(attachment))
	{
		auto iter = m_attachmentEntryIteratorHash.get_entry(attachment.id());
		delete (iter->m_pAttachment);
		(iter->m_pAttachment) = nullptr;
		delete(iter->m_pContainer);
		(iter->m_pContainer) = nullptr;
		m_attachmentEntryContainer.erase(iter);
		m_attachmentEntryIteratorHash.erase(attachment.id());
	}
}

template <typename TElem, typename TElemHandler>
bool
AttachmentPipe<TElem, TElemHandler>::
has_attachment(IAttachment& attachment) const
{
	return m_attachmentEntryIteratorHash.has_entry(attachment.id());
}


template <typename TElem, typename TElemHandler>
template <typename TAttachment>
typename TAttachment::ValueType*
AttachmentPipe<TElem, TElemHandler>::
get_data_array(TAttachment& attachment)
{
	if(has_attachment(attachment))
		return get_data_container(attachment)->get_ptr();

	return nullptr;
}

template <typename TElem, typename TElemHandler>
IAttachmentDataContainer*
AttachmentPipe<TElem, TElemHandler>::
get_data_container(IAttachment& attachment) const
{
	AttachmentEntryIterator iter;
	if(m_attachmentEntryIteratorHash.get_entry(iter, attachment.id()))
		return iter->m_pContainer;
	return nullptr;
}

template <typename TElem, typename TElemHandler>
template <typename TAttachment>
typename TAttachment::ContainerType*
AttachmentPipe<TElem, TElemHandler>::
get_data_container(TAttachment& attachment)
{
	AttachmentEntryIterator iter;
	if(m_attachmentEntryIteratorHash.get_entry(iter, attachment.id()))
		return static_cast<typename TAttachment::ContainerType*>(iter->m_pContainer);
	return nullptr;
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
defragment()
{
	if(!is_fragmented())
		return;

//	if num_elements == 0, then simply resize all data-containers to 0.
	if(num_elements() == 0)
	{
	//	just clear the attachment containers.
		for(auto iter = m_attachmentEntryContainer.begin();
		    iter != m_attachmentEntryContainer.end(); ++iter)
		{
			iter->m_pContainer->resize(0);
		}
		m_stackFreeEntries = UINTStack();
		m_numDataEntries = 0;
	}
	else
	{
	//	calculate the fragmentation array. It has to be of the same size as the fragmented data containers.
		std::vector<uint> vNewIndices(num_data_entries(), ATTACHMENT_CONSTANTS::INVALID_ATTACHMENT_INDEX);

	//	iterate through the elements and calculate the new index of each
		size_t counter = 0;
		typename atraits::element_iterator iter = atraits::elements_begin(m_pHandler);
		typename atraits::element_iterator end = atraits::elements_end(m_pHandler);

		for(; iter != end; ++iter){
			vNewIndices[atraits::get_data_index(m_pHandler, (*iter))] = counter;
			atraits::set_data_index(m_pHandler, (*iter), counter);
			++counter;
		}

	//	after defragmentation there are no free indices.
		m_stackFreeEntries = UINTStack();
		m_numDataEntries = counter;

	//	now iterate through the attached data-containers and defragment each one.
		{
			for(auto iter = m_attachmentEntryContainer.begin(); iter != m_attachmentEntryContainer.end(); ++iter)
			{
				iter->m_pContainer->defragment(&vNewIndices.front(), num_elements());
			}
		}
	}
}

/*
template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
defragment()
{
	if(!is_fragmented())
		return;

//	if num_elements == 0, then simply resize all data-containers to 0.
	if(num_elements() == 0)
	{
	//	just clear the attachment containers.
		for(AttachmentEntryIterator iter = m_attachmentEntryContainer.begin();
			iter != m_attachmentEntryContainer.end(); iter++)
		{
			(*iter).m_pContainer->resize(0);
		}
	}
	else
	{
	//	calculate the fragmentation array. It has to be of the same size as the fragmented data containers.
		std::vector<uint> vNewIndices(num_data_entries(), INVALID_ATTACHMENT_INDEX);

	//	iterate through the elements and calculate the new index of each
		size_t counter = 0;
		{
			typename ElemEntryVec::iterator iter;
			typename ElemEntryVec::iterator copyHere = m_vEntries.begin();
			for(iter = m_vEntries.begin(); iter != m_vEntries.end(); ++iter)
			{
				if(!atraits::entry_is_invalid(*iter))
				{
					vNewIndices[atraits::get_data_index(m_pHandler, (*iter))] = counter;
					atraits::set_data_index(m_pHandler, (*iter), counter);
					++counter;

				//	copy entries
					*copyHere = *iter;
					++copyHere;
				}
			}

		//	resize m_vEntries
			m_vEntries.resize(counter);
		}

	//	now iterate through the attached data-containers and defragment each one.
		{
			for(AttachmentEntryIterator iter = m_attachmentEntryContainer.begin();
						iter != m_attachmentEntryContainer.end(); iter++)
			{
				(*iter).m_pContainer->defragment(&vNewIndices.front(), num_elements());
			}
		}
	}
}
*/
template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
resize_attachment_containers(size_t newSize)
{
//UG_LOG("resizing attachment containers\n");
	//PROFILE_BEGIN(AttachmentResize);
	m_containerSize = newSize;
	for(auto iter = m_attachmentEntryContainer.begin();
	    iter != m_attachmentEntryContainer.end(); ++iter)
	{
		iter->m_pContainer->resize(newSize);
	}
	//PROFILE_END();
//UG_LOG("done\n");
}

template <typename TElem, typename TElemHandler>
void
AttachmentPipe<TElem, TElemHandler>::
grow_attachment_containers(size_t newMinSize)
{
	//if(!m_attachmentEntryContainer.empty())
	{
	//	if the container-size is smaller than newMinSize, we will increase it
	//	by a factor of 2 (at least)
		size_t actSize = get_container_size();

		if(actSize < newMinSize)
		{
			size_t newSize = actSize * 2;
			if(newSize < newMinSize)
				newSize = newMinSize;

			resize_attachment_containers(newSize);
		}
	}
}

template <typename TElem, typename TElemHandler>
size_t
AttachmentPipe<TElem, TElemHandler>::
get_container_size()
{
	return m_containerSize;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	AttachmentAccessor
template <typename TElem, typename TAttachment, typename TElemHandler>
AttachmentAccessor<TElem, TAttachment, TElemHandler>::
AttachmentAccessor() : m_pContainer(nullptr), m_pHandler(nullptr)
{
}

template <typename TElem, typename TAttachment, typename TElemHandler>
AttachmentAccessor<TElem, TAttachment, TElemHandler>::
AttachmentAccessor(const AttachmentAccessor& aa)
{
	m_pContainer = aa.m_pContainer;
	m_pHandler = aa.m_pHandler;
}

template <typename TElem, typename TAttachment, typename TElemHandler>
AttachmentAccessor<TElem, TAttachment, TElemHandler>::
AttachmentAccessor(AttachmentPipe<TElem, TElemHandler>& attachmentPipe, TAttachment& attachment)
{
	m_pContainer = static_cast<typename TAttachment::ContainerType*>(attachmentPipe.get_data_container(attachment));
	m_pHandler = attachmentPipe.get_elem_handler();
	assert(m_pContainer && "ERROR in AttachmentAccessor::AttachmentAccessor(attachmentPipe, attachment): attachment not attached to attachmentPipe!");
}

template <typename TElem, typename TAttachment, typename TElemHandler>
bool
AttachmentAccessor<TElem, TAttachment, TElemHandler>::
access(AttachmentPipe<TElem, TElemHandler>& attachmentPipe, TAttachment& attachment)
{
	if(!attachmentPipe.has_attachment(attachment))
		return false;

	m_pContainer = static_cast<typename TAttachment::ContainerType*>(attachmentPipe.get_data_container(attachment));
	m_pHandler = attachmentPipe.get_elem_handler();
	assert(m_pContainer && "ERROR in AttachmentAccessor::access(attachmentPipe, attachment): attachment not attached to attachmentPipe!");

	return true;
}

}//	end of namespace

#endif
