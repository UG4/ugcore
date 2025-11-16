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

#include <deque>

#include "messaging.h"
#include "bindings_vrl.h"
#include "type_converter.h"
#include "threading.h"
//#include "type_converter.h"
//#include "canvas.h"

namespace ug {
namespace vrl {

void MessageBuffer::addMessage(std::string msg) {

	JNIEnv* env = threading::attachThread(getJavaVM());
	jclass clazz = env->FindClass("edu/gcsc/vrl/ug/UG");
	jmethodID addMessage =
			env->GetStaticMethodID(clazz, "addMessage", "(Ljava/lang/String;)V");

	msg = replaceAll(msg, "\n", "<br>");

	env->CallStaticObjectMethod(clazz, addMessage, stringC2J(env, msg.c_str()));
}

std::vector<std::string> split(const std::string& str, const char delimiter) {
	std::vector<std::string> result;

	result.clear();
	std::stringstream tokenstream;
	tokenstream << str;
	std::string token;

	while (std::getline(tokenstream, token, delimiter)) {
		result.push_back(token);
	}

	return result;
}

std::string replaceAll(
		std::string target,
		const std::string oldstr,
		const std::string newstr) {

	// no substitution necessary
	if (oldstr == newstr) {
		return target;
	}

	for (size_t x = target.find(oldstr); x != std::string::npos; x = target.find(oldstr, x + newstr.size())) {
		target.erase(x, oldstr.length());
		target.insert(x, newstr);
	}

	return target;
}

bool startsWith(std::string str, std::string search) {
	return str.find(search) == 0;
}

bool contains(std::string str, std::string search) {
	return str.find(search) !=std::string::npos;
}

std::string getExceptionMessageString(JNIEnv* env, jthrowable exception) {
	std::string result = "";

	jclass cls;
	jmethodID getMessage;
	
	if (exception != nullptr) {
		cls = env->FindClass("java/lang/Throwable");
		getMessage = env->GetMethodID(cls, "getMessage", "()Ljava/lang/String;");
		jstring msgObj = (jstring)env->CallObjectMethod(exception,getMessage);
		
		if (msgObj != nullptr) {
			result = stringJ2C(env,msgObj);
		}
	}

	return result;
}

bool checkException(JNIEnv* env, std::string msg, bool throwCPPException) {
	jthrowable ex = env->ExceptionOccurred();
	env->ExceptionClear();
	std::string exMsg = getExceptionMessageString(env,ex);

	if (exMsg!="") {
		UG_LOG(msg << " (See Java exception for details)"<< std::endl);
		exMsg = msg + "<font color=\"red\">" + exMsg + "</font>\n";
		
		if (throwCPPException) {
			UG_THROW(exMsg);
		} else {
			UG_LOG(exMsg << std::endl);
		}
	}

	return ex == nullptr;
}


} // end vrl::
} // end ug::



