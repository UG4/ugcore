/*
 * Copyright (c) 2010-2012:  Steinbeis Forschungszentrum (STZ Ölbronn)
 * Author: Michael Hoffer
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

