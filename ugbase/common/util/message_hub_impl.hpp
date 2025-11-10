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

#ifndef __H__UG__message_hub_impl__
#define __H__UG__message_hub_impl__

namespace ug
{


template <class TMsg>
MessageHub::SPCallbackId MessageHub::
register_function_callback(void (*callback)(const TMsg&),
						   bool autoFree)
{
	typedef void (*FuncCallback)(const IMessage&);

	return register_callback_impl<TMsg>((FuncCallback)callback, autoFree);
}


template <class TMsg, class TClass>
MessageHub::SPCallbackId MessageHub::
register_class_callback(TClass* cls,
						void (TClass::*callback)(const TMsg&),
						bool autoFree)
{
	typedef void (TClass::*ClassCallback)(const IMessage&);

	return register_callback_impl<TMsg>(
					boost::bind((ClassCallback)callback, cls, _1),
					autoFree);
}


template <class TMsg>
void MessageHub::
post_message(const TMsg& msg)
{

	size_t id = GetUniqueTypeID<TMsg>();
//	call the callbacks
	CallbackEntryList& callbacks = m_callbackMap[id];
	for(CallbackEntryList::iterator iter = callbacks.begin();
		iter != callbacks.end(); ++iter)
	{
		iter->m_callback(msg);
	}
}


template <class TMsg>
MessageHub::SPCallbackId MessageHub::
register_callback_impl(boost::function<void (const IMessage&)> callback,
					   bool autoFree)
{
	size_t id = GetUniqueTypeID<TMsg>();
	CallbackEntryList& callbacks = m_callbackMap[id];
	CallbackEntryIterator cbIter = callbacks.insert(callbacks.end(),
													CallbackEntry(callback, NULL));

	CallbackId* cbId = new CallbackId(this, id, cbIter, autoFree);

	cbIter->m_callbackId = cbId;

	return SPCallbackId(cbId);
}


}//	end of namespace

#endif
