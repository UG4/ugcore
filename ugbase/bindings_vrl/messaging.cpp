#include "messaging.h"
#include "type_converter.h"
#include "canvas.h"

namespace ug {
	namespace vrl {

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

		void displayMessage(
				std::string title,
				std::string message,
				MessageType mType) {



			JNIEnv *env = Canvas::getInstance()->getJNIEnv();

			if (env->MonitorEnter(Canvas::getInstance()->asJObject()) != JNI_OK) {
				if (env->MonitorExit(Canvas::getInstance()->asJObject()) != JNI_OK) {

				}
			}

			jclass cls = env->FindClass("edu/gcsc/vrl/ug4/UG4");

			if (env->ExceptionCheck()) {
				VRL_DBG("EXEPTION 0:", 1);
				env->ExceptionDescribe();
			}

			jmethodID methodID = env->GetStaticMethodID(cls, "getUG4", "()Ledu/gcsc/vrl/ug4/UG4;");

			if (env->ExceptionCheck()) {
				VRL_DBG("EXEPTION 1:", 1);
				env->ExceptionDescribe();
			}

			jobject ug4 = env->CallObjectMethod(cls, methodID);

			if (env->ExceptionCheck()) {
				VRL_DBG("EXEPTION 2:", 1);
				env->ExceptionDescribe();
			}

			methodID = env->GetMethodID(cls, "addMessage", "(Ljava/lang/String;)V");

			if (env->ExceptionCheck()) {
				VRL_DBG("EXEPTION 3:", 1);
				env->ExceptionDescribe();
			}

			env->CallVoidMethod(ug4, methodID);

			if (env->ExceptionCheck()) {
				VRL_DBG("EXEPTION 4:", 1);
				env->ExceptionDescribe();
			}

			if (env->MonitorExit(Canvas::getInstance()->asJObject()) != JNI_OK) {

			}
		}
	}
}


