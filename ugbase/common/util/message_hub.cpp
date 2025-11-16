/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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
		//	nullptr, to avoid memory access errors when the associated MessageHub
		//	is destroyed.
			m_callbackEntryIter->m_callbackId = nullptr;
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
		//	if the associated callback-id is valid, then set its hub to nullptr
			if(i_entry->m_callbackId != nullptr)
				i_entry->m_callbackId->m_hub = nullptr;
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
	if(cbId->m_hub == nullptr){
		throw(Error("MessageHub::unregister_callback: Invalid callback-id. "
						"The callback was probably already unregistered.",
						MSG_HUB_BAD_CALLBACK_ID));
	}

	UG_ASSERT(cbId->m_hub == this, "Wrong MessageHub");

	CallbackEntryList& callbacks = m_callbackMap[cbId->m_msgTypeId];

//	clear the entry
	callbacks.erase(cbId->m_callbackEntryIter);

//	set the associated hub to nullptr, since it was just unregistered
	cbId->m_hub = nullptr;
}

}// end of namespace
