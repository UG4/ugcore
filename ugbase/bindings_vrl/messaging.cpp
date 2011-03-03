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
	jclass clazz = env->FindClass("edu/gcsc/vrl/ug4/UG4");
	jmethodID addMessage =
			env->GetStaticMethodID(clazz, "addMessage","(Ljava/lang/String;)V");

	msg = replaceAll(msg, "\n", "<br>");

	env->CallStaticObjectMethod(clazz, addMessage, stringC2J(env,msg.c_str()));
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

} // end vrl::
} // end ug::



