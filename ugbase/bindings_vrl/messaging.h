/*
 * File:   messaging.h
 * Author: Michael Hoffer <info@michaelhoffer.de>
 *
 * Created on 15. Oktober 2010, 11:41
 */


#include <string>
#include <deque>

#ifndef MESSAGING_H
#define	MESSAGING_H


namespace ug {
    namespace vrl {
		void soutPrintln(std::string msg);
		void serrPrintln(std::string msg);


		class MessageBuffer {
		public:
//			MessageBuffer();
//			MessageBuffer* MessageBuffer::getInstance();
			static void addMessage(std::string msg);
			static std::string getMessages();
			static std::string messageString;
		private:
			static std::deque<std::string> messages;
			
		};

    } // end vrl::
}// end ug::


//namespace ug {
//	namespace vrl {
//
//		enum MessageType {
//			/**
//			 * defines message as information
//			 */
//			INFO,
//			/**
//			 * defines message as information with singe pulse signal
//			 */
//			INFO_SINGLE,
//			/**
//			 * defines message as warning
//			 */
//			WARNING,
//			/**
//			 * defines message as warning with single pulse signal
//			 */
//			WARNING_SINGLE,
//			/**
//			 * defines message as error
//			 */
//			ERROR,
//			/**
//			 * defines message as error with single pulse signal
//			 */
//			ERROR_SINGLE,
//			/**
//			 * defines message as invisible
//			 */
//			SILENT
//		};
//
//		void displayMessage(
//				std::string title, std::string message, MessageType mType);
//
//	} // vrl::
//} // ug::
//



#endif	/* MESSAGING_H */

