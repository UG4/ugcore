/*
 * Copyright (c) 2010:  Steinbeis Forschungszentrum (STZ Ölbronn)
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
