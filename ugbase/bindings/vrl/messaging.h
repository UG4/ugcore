/*
 * File:   messaging.h
 * Author: Michael Hoffer <info@michaelhoffer.de>
 *
 * Created on 15. Oktober 2010, 11:41
 */


#include <string>
#include <deque>
#include <vector>
#include  <jni.h>

#include "common/ug_config.h"

#ifndef MESSAGING_H
#define	MESSAGING_H

#define EMPHASIZE_BEGIN "\"<tt><b>"
#define EMPHASIZE_END "</tt></b>\""

#define RED_BEGIN "<font color=\"red\">"
#define COLOR_END "</font>"
#define GREEN_BEGIN "<font color=\"green\">"
#define YELLOW_BEGIN "<font color=\"yellow\">"

#define VRL_CRITICAL_ERROR ">> this is a critical error in UG-VRL bindings!" \
					" Please write a bug report to " \
					EMPHASIZE_BEGIN \
					"Michael Hoffer &lt;michael.hoffer@gcsc.uni-frankfurt.de&gt;" \
					EMPHASIZE_END \
					".\n"


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
 *         <code>newstr</code>
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

/**
 * Checks whether <code>str</code> starts with <code>search</code>.
 * @param str string
 * @param search string to search
 * @return <code>true</code> if <code>str</code> starts
 * with <code>search</code>; <code>false</code> otherwise
 */
bool startsWith(std::string str, std::string search);

/**
 * Checks whether <code>str</code> contains <code>search</code>.
 * @param str string
 * @param search string to search
 * @return <code>true</code> if <code>str</code> contains
 * <code>search</code>; <code>false</code> otherwise
 */
bool contains(std::string str, std::string search);

/**
 * Returns the message string of the specified Java exception. If no exception
 * occured, an empty String will be returned.
 * @param env JVM environment
 * @param ex Java exception
 * @return the message string of the specified Java exception or an empty string
 *         if no exception occured
 */
std::string getExceptionMessageString(JNIEnv* env, jthrowable ex);

/**
 * Checks whether an exception has occured. And logs available exception
 * messages.
 * @param env JVM environment
 * @param msg custom error message (shown if exception occurs)
 * @param throwCPPException defines whether to throw a C++ exception (default true)
 * @return <code>true</code> if no exception has occured; <code>false</code> otherwise
 */
bool checkException(JNIEnv* env, std::string msg = "", bool throwCPPException = true);

class UG_API MessageBuffer {
public:
	/**
	 * Adds a message to this message buffer.
	 * @param msg message to add
	 */
	static void addMessage(std::string msg);
	/**
	 * Returns a string containing all messages. Occurences of
	 * <code>&#92;n</code> are replaced with <code>&lt;br&gt;</code>.
	 * @return all messages as HTML compatible string
//	 */

private:
//	static std::deque<std::string> messages;
//	static std::string messageString;
//	static unsigned int maxQueueSize;

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

