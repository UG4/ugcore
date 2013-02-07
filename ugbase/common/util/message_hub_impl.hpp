// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.11.2011 (m,d,y)

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
