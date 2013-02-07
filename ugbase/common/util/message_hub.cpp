// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 18.11.2011 (m,d,y)
 
#include "message_hub.h"
#include "common/assert.h"
#include "common/error.h"

namespace ug{

MessageHub::CallbackEntry::
CallbackEntry(const Callback& cb, CallbackId* cbId) :
	m_callback(cb),
	m_callbackId(cbId)
{
}

MessageHub::CallbackId::
CallbackId(MessageHub* hub, size_t msgTypeId,
		   CallbackEntryIterator callbackEntryIter, bool autoFree) :
	m_hub(hub),
	m_msgTypeId(msgTypeId),
	m_callbackEntryIter(callbackEntryIter),
	m_autoFree(autoFree)
{
}

MessageHub::CallbackId::~CallbackId()
{
//	Make sure that the associated message hub still exists
	if(m_hub){
		if(m_autoFree)
			m_hub->unregister_callback_impl(this);
		else{
		//	we have to set the callback-id of the associated callback-entry to
		//	NULL, to avoid memory access errors when the associated MessageHub
		//	is destroyed.
			m_callbackEntryIter->m_callbackId = NULL;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	MessageHub implementation
MessageHub::MessageHub()
{
}

MessageHub::~MessageHub()
{
//	we have to make sure to invalidate all associated callback-ids.
//	All entry-lists have to be deleted
	for(CallbackMap::iterator i_table = m_callbackMap.begin();
		i_table != m_callbackMap.end(); ++i_table)
	{
		CallbackEntryList& entryList = i_table->second;
		for(CallbackEntryList::iterator i_entry = entryList.begin();
			i_entry != entryList.end(); ++i_entry)
		{
		//	if the associated callback-id is valid, then set its hub to NULL
			if(i_entry->m_callbackId != NULL)
				i_entry->m_callbackId->m_hub = NULL;
		}
	}
}

void MessageHub::
unregister_callback(MessageHub::SPCallbackId cbId)
{
	unregister_callback_impl(cbId.get());
}

void MessageHub::
unregister_callback_impl(MessageHub::CallbackId* cbId)
{
	if(cbId->m_hub == NULL){
		throw(Error("MessageHub::unregister_callback: Invalid callback-id. "
						"The callback was probably already unregistered.",
						MSG_HUB_BAD_CALLBACK_ID));
	}

	UG_ASSERT(cbId->m_hub == this, "Wrong MessageHub");

	CallbackEntryList& callbacks = m_callbackMap[cbId->m_msgTypeId];

//	clear the entry
	callbacks.erase(cbId->m_callbackEntryIter);

//	set the associated hub to NULL, since it was just unregistered
	cbId->m_hub = NULL;
}

}// end of namespace
