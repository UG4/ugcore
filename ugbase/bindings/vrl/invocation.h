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

#include <map>
#include <sstream>
#include "registry/registry.h"
#include <jni.h>

#include"bindings_vrl.h"

#ifndef INVOKATION_H
#define	INVOKATION_H

namespace ug {
namespace vrl {
namespace invocation {

//const std::string createMethodSignature(JNIEnv* env, const char* className,
//		const char* methodName, bool readOnly, jobjectArray params);


/**
 * Returns an exported function by its signature.
 * @param env JVM environment to operate on
 * @param reg registry to search
 * @param eCls exported class
 * @param params Java object array containing method parameter
 * @return pointer to requested method if such a method exists;
 *         <code>NULL</code> otherwise
 */
const ug::bridge::ExportedConstructor* getConstructorBySignature(
		JNIEnv *env,
		ug::bridge::Registry* reg,
		ug::bridge::IExportedClass* eCls,
		jobjectArray params);

/**
 * Returns an exported method by its signature.
 * @param env JVM environment to operate on
 * @param reg registry to search
 * @param clazz class the method to search belongs to
 * @param readOnly defines whether to search a const method
 * @param methodName method name
 * @param params Java object array containing method parameter
 * @return pointer to requested method if such a method exists;
 *         <code>NULL</code> otherwise
 */
const ug::bridge::ExportedMethod* getMethodBySignature(
		JNIEnv *env,
		ug::bridge::Registry* reg,
		const ug::bridge::IExportedClass* clazz,
		bool readOnly,
		std::string methodName,
		jobjectArray params);

/**
 * Returns an exported function by its signature.
 * @param env JVM environment to operate on
 * @param reg registry to search
 * @param functionName method name
 * @param params Java object array containing method parameter
 * @return pointer to requested method if such a method exists;
 *         <code>NULL</code> otherwise
 */
const ug::bridge::ExportedFunction* getFunctionBySignature(
		JNIEnv *env,
		ug::bridge::Registry* reg,
		std::string functionName,
		jobjectArray params);

/**
 * Returns an exported class by name.
 * @param reg registry to search
 * @param className class name
 * @return pointer to requested class if such a class exists;
 *         <code>NULL</code> otherwise
 */
const ug::bridge::IExportedClass* getExportedClassPtrByName(
		ug::bridge::Registry* reg,
		std::string className);

/**
 * Returns a classnode by name.
 * @param reg registry to search
 * @param className class name
 * @return pointer to requested classnode if such a node exists;
 *         <code>NULL</code> otherwise
 */
const ug::bridge::ClassNameNode* getClassNodePtrByName(
		ug::bridge::Registry* reg,
		std::string className);

void initClasses(ug::bridge::Registry &reg);

} // invocation::
} // vrl::
} // ug::

#endif
