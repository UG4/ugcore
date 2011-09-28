///*
// * File:   canvas.cpp
// * Author: miho
// *
// * Created on 15. Oktober 2010, 13:37
// */
//
//#include "canvas.h"
//#include "type_converter.h"
//
//namespace ug {
//	namespace vrl {
//
//		Canvas::Canvas() {
//		}
//
//		Canvas* Canvas::getInstance() {
//			static Canvas canvas;
//			return &canvas;
//		}
//
//		void Canvas::setJObject(JNIEnv *env, jobject obj) {
//			_canvasJ = obj;
//			_env = env;
//		}
//
//		jobject Canvas::asJObject() {
//			return _canvasJ;
//		}
//
//		void Canvas::addObject(jobject obj) {
//			jclass canvasClass = _env->FindClass(
//					"eu/mihosoft/vrl/reflection/VisualCanvas");
//			jmethodID methodID = _env->GetMethodID(canvasClass,"addObject",
//					"(Ljava/lang/Object;)Leu/mihosoft/vrl/reflection/VisualObject;");
//
//			if (_env->ExceptionCheck()) {
//				VRL_DBG("EXEPTION 0:", 1);
//				_env->ExceptionDescribe();
//			}
//
//			_env->CallObjectMethod(asJObject(),methodID,obj);
//
//			if (_env->ExceptionCheck()) {
//				VRL_DBG("EXEPTION 1:", 1);
//				_env->ExceptionDescribe();
//			}
//
//		}
//
//		jobject Canvas::getMessageBox() {
//
//			VRL_DBG("INIT", 1);
//
//			if (_env->ExceptionCheck()) {
//				VRL_DBG("EXEPTION 0:", 1);
//				_env->ExceptionDescribe();
//			}
//
//			jclass canvasClass = _env->FindClass(
//					"eu/mihosoft/vrl/reflection/VisualCanvas");
//
//			if (_env->ExceptionCheck()) {
//				VRL_DBG("EXEPTION 1:", 1);
//				_env->ExceptionDescribe();
//			}
//
//			VRL_DBG("AFTER_GET_CLASS", 1);
//
//			if (_env->ExceptionCheck()) {
//				VRL_DBG("EXEPTION 2:", 1);
//				_env->ExceptionDescribe();
//			}
//
//			jmethodID methodID = _env->GetMethodID(canvasClass,
//					"getMessageBox", "()Leu/mihosoft/vrl/visual/MessageBox;");
//
//			VRL_DBG("AFTER_GET_METHOD_ID", 1);
//
//			return _env->CallObjectMethod(asJObject(), methodID);
//		}
//
//		JNIEnv* Canvas::getJNIEnv() {
//			return _env;
//		}
//	} // vrl::
//} // ug::
//
//
//
//
