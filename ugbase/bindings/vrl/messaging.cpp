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
	
	if (exception != NULL) {
		cls = env->FindClass("java/lang/Throwable");
		getMessage = env->GetMethodID(cls, "getMessage", "()Ljava/lang/String;");
		jstring msgObj = (jstring)env->CallObjectMethod(exception,getMessage);
		
		if (msgObj != NULL) {
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

	return ex == NULL;
}


} // end vrl::
} // end ug::



