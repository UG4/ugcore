// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.11.2011 (m,d,y)

#ifndef __H__UG__message_hub_impl__
#define __H__UG__message_hub_impl__

namespace ug
{

template <class TMsg>
int MessageHub::
get_message_id(const char* messageIdName)
{
//	get the type id
	size_t msgTypeId = GetUniqueTypeID<TMsg>();

//	check whether the messageID was already registered
	IdMap::iterator iter = m_idMap.find(messageIdName);
	if(iter != m_idMap.end()){
		int msgId = iter->second;
	//	compare types
		if(m_msgTypeIds[msgId] != msgTypeId){
			throw(Error("MessageHub::get_message_id: message-name was "
						"already registered with a different message type!",
						MSG_HUB_TYPE_MISMATCH));
		}
		else
			return msgId;
	}

//	it wasn't yet registered. Do this now.
	int newId = m_highestMsgId++;
	m_idMap[messageIdName] = newId;
	m_callbackTable.push_back(SPCallbackList(new CallbackList()));
	m_msgTypeIds.push_back(msgTypeId);

	return newId;
}


template <class TMsg>
MessageHub::SPCallbackId MessageHub::
register_callback(int msgId, void (*callback)(int, const TMsg*), bool autoFree)
{
	typedef void (*FuncCallback)(int, const IMessage*);

	return register_callback_impl<TMsg>(msgId, (FuncCallback)callback, autoFree);
}


template <class TMsg, class TClass>
MessageHub::SPCallbackId MessageHub::
register_callback(int msgId, TClass* cls,
				  void (TClass::*callback)(int, const TMsg*),
				  bool autoFree)
{
	typedef void (TClass::*ClassCallback)(int, const IMessage*);

	return register_callback_impl<TMsg>(msgId,
					boost::bind((ClassCallback)callback, cls, _1, _2),
					autoFree);
}


template <class TMsg>
void MessageHub::
post_message(int msgId, const TMsg* msg)
{
//	check whether the given msgId is valid
	if((msgId < 0) || (msgId >= m_highestMsgId)){
		throw(Error("MessageHub::post_message: Bad message id!",
					MSG_HUB_BAD_MESSAGE_ID));
	}

	size_t msgTypeId = GetUniqueTypeID<TMsg>();
	if(msgTypeId != m_msgTypeIds[msgId]){
		throw(Error("MessageHub::post_message: message type mismatch!",
					MSG_HUB_TYPE_MISMATCH));
	}

//	call the callbacks
	CallbackList& callbacks = *m_callbackTable[msgId].get_impl();
	for(CallbackList::iterator iter = callbacks.begin();
		iter != callbacks.end(); ++iter)
	{
		(*iter)(msgId, msg);
	}
}


template <class TMsg>
MessageHub::SPCallbackId MessageHub::
register_callback_impl(int msgId,
					   boost::function<void (int, const IMessage*)> callback,
					   bool autoFree)
{
//	check whether the given msgId is valid
	if((msgId < 0) || (msgId >= m_highestMsgId)){
		throw(Error("MessageHub::register_callback: Bad message id!",
					MSG_HUB_BAD_MESSAGE_ID));
	}

	size_t msgTypeId = GetUniqueTypeID<TMsg>();
	if(msgTypeId != m_msgTypeIds[msgId]){
		throw(Error("MessageHub::register_callback: message type mismatch!",
					MSG_HUB_TYPE_MISMATCH));
	}

//	everything is valid. Register the callback. Check whether the vector
//	contains any empty callbacks. Use such a slot, if one exists.
//todo: implement a speed up
	CallbackList& callbacks = *m_callbackTable[msgId].get_impl();
	CallbackList::iterator cbIter = callbacks.insert(callbacks.end(), callback);

	return SPCallbackId(new CallbackId(this, msgId, cbIter, autoFree));
}


}//	end of namespace

#endif
