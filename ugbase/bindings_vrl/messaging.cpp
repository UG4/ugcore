#include <deque>

#include "messaging.h"
#include "bindings_vrl.h"
#include "type_converter.h"
//#include "type_converter.h"
//#include "canvas.h"

namespace ug {
	namespace vrl {

		std::deque<std::string> MessageBuffer::messages;

		void MessageBuffer::addMessage(std::string msg) {
			messages.push_front(msg);
			messages.resize(100, std::string());
		}

		std::string MessageBuffer::getMessages() {
			std::string result = "";

			for (std::deque<std::string>::iterator i = messages.begin(); i != messages.end(); i++) {
				result += *i;
			}

			result = replaceAll(result,"\n","<br>");

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

			for (size_t x = target.find(oldstr); x != std::string::npos; x = target.find(oldstr)) {
				target.erase(x, oldstr.length());
				target.insert(x, newstr);
			}

			return target;
		}

		void soutPrintln(std::string message) {

			//			std::string msg = message;
			//
			//			JNIEnv* env = getJNIEnv();
			//
			//			if (env->ExceptionCheck()) {
			//				env->ExceptionDescribe();
			//
			//			}
			//
			//			jclass systemClass = env->FindClass("java/lang/System");
			//
			//			std::cout << "ERROR: 1" << msg << std::endl;
			//
			//			if (env->ExceptionCheck()) {
			//				env->ExceptionDescribe();
			//
			//			}
			//
			//			jfieldID fieldID = env->GetStaticFieldID(
			//					systemClass, "out", "Ljava/io/PrintStream;");
			//
			//			std::cout << "ERROR: 2" << msg << std::endl;
			//
			//			if (env->ExceptionCheck()) {
			//				env->ExceptionDescribe();
			//			}
			//
			//			jobject outStream = env->GetStaticObjectField(
			//					systemClass, fieldID);
			//
			//			std::cout << "ERROR: 3" << msg << std::endl;
			//
			//			if (env->ExceptionCheck()) {
			//				env->ExceptionDescribe();
			//
			//			}
			//
			//			jclass streamClass = env->GetObjectClass(outStream);
			//			jmethodID methodID = env->GetMethodID(
			//					streamClass, "println", "(Ljava/lang/String;)V");
			//
			//			std::cout << "ERROR: 4" << msg << std::endl;
			//
			//			if (env->ExceptionCheck()) {
			//				env->ExceptionDescribe();
			//
			//			}
			//
			//			env->CallVoidMethod(
			//					outStream, methodID, string2JObject(env, msg.c_str()));
			//
			//			std::cout << "ERROR: 5" << msg << std::endl;
		}

		void serrPrintln(std::string msg) {

			//			JNIEnv* env = getJNIEnv();
			//			jclass systemClass = env->FindClass("java/lang/System");
			//
			//			jfieldID fieldID = env->GetStaticFieldID(
			//					systemClass, "err", "Ljava/io/PrintStream;");
			//
			//			jobject outStream = env->GetStaticObjectField(
			//					systemClass, fieldID);
			//
			//			jclass streamClass = env->GetObjectClass(outStream);
			//			jmethodID methodID = env->GetMethodID(
			//					streamClass, "println", "(Ljava/lang/String;)V");
			//
			//			env->CallVoidMethod(
			//					outStream, methodID, string2JObject(env, msg));
		}

		//		void displayMessage(
		//				std::string title,
		//				std::string message,
		//				MessageType mType) {
		//
		//			VRL_DBG("INIT", 1);
		//
		//			JNIEnv *env = Canvas::getInstance()->getJNIEnv();
		//
		//			VRL_DBG("AFTER_GET_INSTANCE", 1);
		//
		//			jobject messageBox = Canvas::getInstance()->getMessageBox();
		//
		//			VRL_DBG("AFTER_GET_MESSAGE_BOX", 1);
		//
		//			jclass messageBoxClass =
		//					env->FindClass("eu/mihosoft/vrl/visual/MessageBox");
		//
		//			VRL_DBG("AFTER_GET_MSG_CLASS", 1);
		//
		//
		//
		//			jmethodID methodID = env->GetMethodID(messageBoxClass,
		//					"addMessage",
		//					"(Ljava/lang/String;Ljava/lang/String;Leu/mihosoft/vrl/visual/MessageType;)Leu/mihosoft/vrl/visual/Message;");
		//
		//			if (env->ExceptionCheck()) {
		//				VRL_DBG("EXEPTION 0:",1);
		//				env->ExceptionDescribe();
		//			}
		//
		//			VRL_DBG("AFTER_GET_METHOD_ID", 1);
		//
		//			jobject msgType = messageTypeC2J(env, mType);
		//
		//			VRL_DBG("AFTER_GET_MSGTYPE", 1);
		//
		//			if (env->ExceptionCheck()) {
		//				VRL_DBG("EXEPTION 1:",1);
		//				env->ExceptionDescribe();
		//			}
		//
		//			VRL_DBG(getClassName(env,messageBox),1);
		//
		//			jobject titleStr = stringC2J(env, title);
		//			jobject msgStr = stringC2J(env, message);
		//
		//			jobject msg = env->CallObjectMethod(
		//					messageBox,
		//					methodID,
		//					titleStr,
		//					msgStr, msgType);
		//
		//			if (env->ExceptionCheck()) {
		//				VRL_DBG("EXEPTION 2:",1);
		//				env->ExceptionDescribe();
		//			}
		//
		//			VRL_DBG("AFTER_CALL_ADD_MSG", 1);
		//		}
		//	}

		//		void displayMessage(
		//				std::string title,
		//				std::string message,
		//				MessageType mType) {
		//
		//
		//
		//			JNIEnv *env = Canvas::getInstance()->getJNIEnv();
		//
		//			if (env->MonitorEnter(Canvas::getInstance()->asJObject()) != JNI_OK) {
		//				if (env->MonitorExit(Canvas::getInstance()->asJObject()) != JNI_OK) {
		//
		//				}
		//			}
		//
		//			jclass cls = env->FindClass("edu/gcsc/vrl/ug4/UG4");
		//
		//			if (env->ExceptionCheck()) {
		//				VRL_DBG("EXEPTION 0:", 1);
		//				env->ExceptionDescribe();
		//			}
		//
		//			jmethodID methodID = env->GetStaticMethodID(cls, "getUG4", "()Ledu/gcsc/vrl/ug4/UG4;");
		//
		//			if (env->ExceptionCheck()) {
		//				VRL_DBG("EXEPTION 1:", 1);
		//				env->ExceptionDescribe();
		//			}
		//
		//			jobject ug4 = env->CallObjectMethod(cls, methodID);
		//
		//			if (env->ExceptionCheck()) {
		//				VRL_DBG("EXEPTION 2:", 1);
		//				env->ExceptionDescribe();
		//			}
		//
		//			methodID = env->GetMethodID(cls, "addMessage", "(Ljava/lang/String;)V");
		//
		//			if (env->ExceptionCheck()) {
		//				VRL_DBG("EXEPTION 3:", 1);
		//				env->ExceptionDescribe();
		//			}
		//
		//			env->CallVoidMethod(ug4, methodID);
		//
		//			if (env->ExceptionCheck()) {
		//				VRL_DBG("EXEPTION 4:", 1);
		//				env->ExceptionDescribe();
		//			}
		//
		//			if (env->MonitorExit(Canvas::getInstance()->asJObject()) != JNI_OK) {
		//
		//			}
		//		}
	} // end vrl::
} // end ug::



