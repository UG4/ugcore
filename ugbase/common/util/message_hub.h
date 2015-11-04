#ifndef __H__UG__message_hub__
#define __H__UG__message_hub__

#include <vector>
#include <list>
#include <string>
#include <map>
#include <boost/function.hpp>
#include <boost/function_equal.hpp>
#include <boost/bind.hpp>
#include "common/assert.h"
#include "common/util/smart_pointer.h"
#include "common/util/metaprogramming_util.h"
#include "common/error.h"
#include "common/log.h"

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

class MessageHub;
typedef SmartPtr<MessageHub> SPMessageHub;

///	Allows to register callbacks and post messages to those callbacks.
/**
 * The MessageHub provides a flexible way to register callbacks for a given
 * message-type, which can be called through the MessageHub::post_message method.
 *
 * Functions and member-functions of arbitrary classes can both be used as a
 * callback. Callback functions and methods have to have the following signature:
 *
 * void (*)(const TMsg&)
 *
 * where the message-type TMsg is derived from MessageHub::IMessage.
 * They can be registered through the register_callback methods.
 * If you're intending to register a member method of a class, then you have to use
 * additionally pass the pointer to the class instance, whose member-method shall
 * be called.
 *
 * On registration of a callback, an identifier is returned, which allows to
 * unregister the callback later on. The identifier is wrapped in a smart-pointer,
 * and supports an auto-free property. If this is property is enabled, then
 * the associated callback is automatically unregistered, as soon as the last
 * copy of the encapsulating smart-pointer is deleted.
 */
class MessageHub
{
	public:
	//	predeclarations
		class IMessage;
		class CallbackId;

	private:
	//	private type definitions
		typedef boost::function<void (const IMessage&)> Callback;

	///	The CallbackEntry holds the actual callback and the associated callback-id.
	/**	If a MessageHub is destroyed before all callbacks are unregistered, this
	 * information can be used to notify associated callback-ids, that the associated
	 * hub is gone. Important if autoFree is enabled.
	 */
		struct CallbackEntry{
			CallbackEntry(const Callback& cb, CallbackId* cbId);
			Callback	m_callback;
			CallbackId*	m_callbackId;
		};

		typedef std::list<CallbackEntry>	CallbackEntryList;
		typedef CallbackEntryList::iterator	CallbackEntryIterator;
		typedef std::map<size_t, CallbackEntryList>	CallbackMap;

	public:
	///	Error codes which give information on the error-reason
		enum ErrorIds{
			MSG_HUB_UNKNOWN_ERROR,
			MSG_HUB_TYPE_MISMATCH,
			MSG_HUB_BAD_MESSAGE_ID,
			MSG_HUB_BAD_CALLBACK_ID
		};

	///	Instances of this class are thrown if an error occurs in MessageHub
	/**	This class derives from UGError.
	 * Use MessageHub::Error::get_message_hub_error_id() to get a more information
	 * on what caused the error.
	 */
		class Error : public UGError
		{
			public:
				Error(const char* msg, ErrorIds errorId) :
					UGError(msg), m_errorId(errorId) {}

				ErrorIds get_message_hub_error_id()	{return m_errorId;}

			protected:
				ErrorIds m_errorId;
		};

	///	This is the base class of all messages, which may be passed to callbacks.
		class IMessage{
			public:
				IMessage()	{}
				virtual ~IMessage()	{}
		};

	///	The callback-id allows to deregister previously registered callbacks.
	/**	Note that the class features an autoFree mechanism, which automatically
	 * frees the associated callback. Since this class is always wrapped in a
	 * smart-pointer, the the associated callback won't be freed, until the last
	 * copy of that smart-pointer is deleted.
	 * You can disable the auto-free mechanism through the set set_auto_free method
	 * or by setting the autoFree parameter of the MessageHub::register_callback
	 * methods to false.
	 */
		class CallbackId{
			friend class MessageHub;

			public:
				~CallbackId();
				void set_auto_free(bool autoFree)	{m_autoFree = autoFree;}

			private:
				CallbackId(MessageHub* hub, size_t msgTypeId,
						   CallbackEntryIterator callbackEntryIter, bool autoFree);

				MessageHub*				m_hub;
			///	Make sure to only access the iterator while m_hub != NULL.
				size_t					m_msgTypeId;
				CallbackEntryIterator	m_callbackEntryIter;
				bool					m_autoFree;
		};

		typedef SmartPtr<CallbackId>	SPCallbackId;

	public:
		MessageHub();
		~MessageHub();

	///	registers a function callback given a message-type.
	/** The callback has to be of the type
	 *
	 * void (*FuncCallback)(const TMsg&)
	 *
	 * where TMsg is a type derived from MessageHub::IMessage.
	 *
	 * The method returns a smart-pointer to a callback-identifier.
	 * The auto-free property is disabled for the returned callback-id by default.
	 * Note that this behavior differs from the similar register_class_callback
	 * method for class-methods.*/
		template <class TMsg>
		SPCallbackId register_function_callback(void (*callback)(const TMsg&),
											   bool autoFree = false);

	///	registers a method callback given a message-type.
	/** The callback has to be of the type
	 *
	 * void (TClass::*ClassCallback)(const TMsg&)
	 *
	 * where TClass is the class whose member function is registered as callback
	 * and TMsg a type derived from MessageHub::IMessage.
	 *
	 * The method returns a smart-pointer to a callback-identifier. When the
	 * instance is deleted, the callback is unregistered by default.
	 * Note that this behavior differs from the similar register_function_callback
	 * method for function pointers.
	 * It's a good idea to store the smart-pointer as a member in the class from
	 * which you register the callback (if it is registered from a class at all).
	 * You then won't have to deal with unregistration manually.*/
		template <class TMsg, class TClass>
		SPCallbackId register_class_callback(TClass* cls,
										   void (TClass::*callback)(const TMsg&),
										   bool autoFree = true);

	///	Call this method to explicitly unregister a callback.
	/**	Note that if you're storing the callback-id in a class and if autoFree
	 * is enabled for the callback-id, then the callback is automatically
	 * unregistered, when the last instance of the smart-pointer is deleted.
	 *
	 * If you use this method, the autoFree property of the given CallbackId
	 * will be automatically disabled.
	 */
		void unregister_callback(SPCallbackId cbId);

	///	Posts a message to all callbacks which are registered for the given message tpye
		template <class TMsg>
		void post_message(const TMsg& msg);

	private:
	///	registers a callback given a message-id.
	/**	Make sure to only pass msgIds which were retrieved through get_message_id
	 * before. Also be sure to use the correct msg-type, which was registered with
	 * the given message-id.
	 *
	 * The callback has to be of the type
	 *
	 * <tt>boost::function\<void (const IMessage&)\></tt>
	 *
	 * The method returns a smart-pointer to an callback-identifier.
	 *
	 * If the message-id was not registered or if it was registered with a
	 * different type, then an instance of MessageHub::Error is thrown
	 * (derives from UGError).
	 */
		template <class TMsg>
		SPCallbackId register_callback_impl(
							   boost::function<void (const IMessage&)> callback,
							   bool autoFree);

	///	performs unregistration of the given callback
		void unregister_callback_impl(CallbackId* cbId);

	private:
		CallbackMap	m_callbackMap;///< given a msg-type-id, this map returns a list of associated callbacks
};

// end group ugbase_common_util
/// \}

}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "message_hub_impl.hpp"

#endif
