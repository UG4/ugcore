//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d06

#ifndef __UTIL__ATTACHMENT_PIPE__IMPL__
#define __UTIL__ATTACHMENT_PIPE__IMPL__

#include "attachment_pipe.h"
#include "../common.h"

#define ATRAITS attachment_traits<TElem, TElemHandlerTag>

namespace ug
{

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
clear()
{
	clear_elements();
	clear_attachments();
}

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
clear_elements()
{
	m_vEntries.resize(0);
	m_stackFreeEntries = UINTStack();
	resize_attachment_containers(0);
	m_numElements = 0;
}

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
clear_attachments()
{
	m_attachmentEntryContainer.clear();
	m_attachmentEntryIteratorHash.clear();
}

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
reset_values(uint dataIndex)
{
	for(AttachmentEntryIterator iter = m_attachmentEntryContainer.begin();
					iter != m_attachmentEntryContainer.end(); iter++)
	{
		(*iter).m_pContainer->reset_entry(dataIndex);
	}
}

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
register_element(TElem elem)
{
	uint newInd;
//	check if there is a free entry

	if(!m_stackFreeEntries.empty())
	{
	//	take it
		newInd = m_stackFreeEntries.top();
		m_stackFreeEntries.pop();
		m_vEntries[newInd] = elem;

	//	initialize attachments with default constructor
		reset_values(newInd);
	}
	else
	{
	//	add a new entry
		newInd = m_vEntries.size();
	//	push the elem to the end of the elem-vec.
		m_vEntries.push_back(elem);
	//	make sure that the data containers are big enough.
		grow_attachment_containers(m_vEntries.size());
	}

//	assign new attachment index to the element
	ATRAITS::set_data_index(elem, newInd);

//	increase element count
	m_numElements++;
}

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
unregister_element(const TElem& elem)
{
//	get the attachment index
	uint ind = ATRAITS::get_data_index(elem);
//	clear the entry
	m_vEntries[ind] = NULL;
//	store the index in the stack of free entries
	m_stackFreeEntries.push(ind);
//	decrease element count
	m_numElements--;
}


template <class TElem, class TElemHandlerTag>
template <class TAttachment>
void
AttachmentPipe<TElem, TElemHandlerTag>::
attach(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, uint options)
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
		AttachmentEntryIterator iter = m_attachmentEntryContainer.insert(m_attachmentEntryContainer.end(), AttachmentEntry(pClonedAttachment, pClonedContainer, options));
		m_attachmentEntryIteratorHash.add(iter, pClonedAttachment->id());
	}
	else
	{
		LOG("WARNING: attachment " << attachment.get_name() << " already attached." << std::endl);
	}
}

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
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
		AttachmentEntryIterator iter = m_attachmentEntryContainer.insert(m_attachmentEntryContainer.end(), AttachmentEntry(pClonedAttachment, pClonedContainer, options));
		m_attachmentEntryIteratorHash.add(iter, pClonedAttachment->id());
	}
	else
	{
		LOG("WARNING: attachment " << attachment.get_name() << " already attached." << std::endl);
	}
}

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
detach(IAttachment& attachment)
{
	if(has_attachment(attachment))
	{
		AttachmentEntryIterator iter = m_attachmentEntryIteratorHash.first(attachment.id());
		SAFE_DELETE((*iter).m_pAttachment);
		SAFE_DELETE((*iter).m_pContainer);
		m_attachmentEntryContainer.erase(iter);
		m_attachmentEntryIteratorHash.erase(attachment.id());
	}
}

template <class TElem, class TElemHandlerTag>
bool
AttachmentPipe<TElem, TElemHandlerTag>::
has_attachment(IAttachment& attachment) const
{
	return m_attachmentEntryIteratorHash.has_entries(attachment.id());
}

template <class TElem, class TElemHandlerTag>
IAttachmentDataContainer*
AttachmentPipe<TElem, TElemHandlerTag>::
get_data_container(IAttachment& attachment) const
{
	if(has_attachment(attachment))
		return (*m_attachmentEntryIteratorHash.first(attachment.id())).m_pContainer;
	return NULL;
}

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
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
		uint counter = 0;
		{
			typename ElemEntryVec::iterator iter;
			for(iter = m_vEntries.begin(); iter != m_vEntries.end(); ++iter)
			{
				if((*iter) != NULL)
				{
					vNewIndices[ATRAITS::get_data_index((*iter))] = counter;
					ATRAITS::set_data_index((*iter), counter);
					++counter;
				}
			}
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

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
resize_attachment_containers(uint newSize)
{
	for(AttachmentEntryIterator iter = m_attachmentEntryContainer.begin();
				iter != m_attachmentEntryContainer.end(); iter++)
	{
		(*iter).m_pContainer->resize(newSize);
	}
}

template <class TElem, class TElemHandlerTag>
void
AttachmentPipe<TElem, TElemHandlerTag>::
grow_attachment_containers(uint newMinSize)
{

	//if(!m_attachmentEntryContainer.empty())
	{
	//	if the container-size is smaller than newMinSize, we will increase it
	//	by a factor of 2 (at least)
		uint actSize = get_container_size();

		if(actSize < newMinSize)
		{
			uint newSize = actSize * 2;
			if(newSize < newMinSize)
				newSize = newMinSize;

			for(AttachmentEntryIterator iter = m_attachmentEntryContainer.begin();
				iter != m_attachmentEntryContainer.end(); iter++)
			{
				(*iter).m_pContainer->resize(newSize);
			}
		}
	}
}

template <class TElem, class TElemHandlerTag>
uint
AttachmentPipe<TElem, TElemHandlerTag>::
get_container_size()
{
	if(!m_attachmentEntryContainer.empty())
	{
		AttachmentEntryIterator iter = m_attachmentEntryContainer.begin();
		return (*iter).m_pContainer->size();
	}
	return m_vEntries.size();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	AttachmentAccessor
template <class TElem, class TAttachment, class ElemHandlerTag>
AttachmentAccessor<TElem, TAttachment, ElemHandlerTag>::
AttachmentAccessor() : m_pContainer(NULL)
{
}

template <class TElem, class TAttachment, class ElemHandlerTag>
AttachmentAccessor<TElem, TAttachment, ElemHandlerTag>::
AttachmentAccessor(const AttachmentAccessor& aa)
{
	m_pContainer = aa.m_pContainer;
}

template <class TElem, class TAttachment, class ElemHandlerTag>
AttachmentAccessor<TElem, TAttachment, ElemHandlerTag>::
AttachmentAccessor(const AttachmentPipe<TElem, ElemHandlerTag>& attachmentPipe, TAttachment& attachment)
{
	m_pContainer = static_cast<typename TAttachment::ContainerType*>(attachmentPipe.get_data_container(attachment));
	assert(m_pContainer && "ERROR in AttachmentAccessor::AttachmentAccessor(attachmentPipe, attachment): attachment not attached to attachmentPipe!");
}

template <class TElem, class TAttachment, class ElemHandlerTag>
void
AttachmentAccessor<TElem, TAttachment, ElemHandlerTag>::
access(const AttachmentPipe<TElem, ElemHandlerTag>& attachmentPipe, TAttachment& attachment)
{
	m_pContainer = static_cast<typename TAttachment::ContainerType*>(attachmentPipe.get_data_container(attachment));
	assert(m_pContainer && "ERROR in AttachmentAccessor::access(attachmentPipe, attachment): attachment not attached to attachmentPipe!");
}

}//	end of namespace

#endif
