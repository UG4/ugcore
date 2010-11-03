/*
 * File:   messaging.h
 * Author: Michael Hoffer <info@michaelhoffer.de>
 *
 * Created on 15. Oktober 2010, 11:41
 */


#include <string>
#include <deque>
#include <vector>

#ifndef MESSAGING_H
#define	MESSAGING_H

#define EMPHASIZE_BEGIN "\"<tt><b>"
#define EMPHASIZE_END "</tt></b>\""


namespace ug {
	namespace vrl {
		/**
		 * DO NOT USE THIS METHOD!
		 * (totally broken)
		 * @param msg
		 */
		void soutPrintln(std::string msg);
		/**
		 * DO NOT USE THIS METHOD!
		 * (totally broken)
		 * @param msg
		 */
		void serrPrintln(std::string msg);

		/**
		 * Replaces each substring of <code>target</code> string that is equal to
		 * <code>oldstr</code> with <code>newstr</code>
		 * @param target string to modify
		 * @param oldstr string to raplace
		 * @param newstr replacement string
		 * @return a copy of the specified <code>target</code> string where
		 *         all occurences of <code>oldstr</code> are replaced with
		 *         <code>newstr</code
		 */
		std::string replaceAll(
				std::string target,
				const std::string oldstr,
				const std::string newstr);

		/**
		 * Splits string around given delimiter.
         * @param str string to split
         * @param delimiter delimiter
		 * @result vector of strings computed by splitting this
		 *        string around matches of the given delimiter
         */
		std::vector<std::string> split(
				const std::string& str,
				const char delimiter);

		class MessageBuffer {
		public:
			/**
			 * Adds a message to this message buffer.
			 * @param msg message to add
			 */
			static void addMessage(std::string msg);
			/**
			 * Returns a string containing all messages. Occurences of
			 * <code>\n</code> are replaced with <code><br></code>.
			 * @return all messages as HTML compatible string
			 */
			static std::string getMessages();

		private:
			static std::deque<std::string> messages;
			static std::string messageString;

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

