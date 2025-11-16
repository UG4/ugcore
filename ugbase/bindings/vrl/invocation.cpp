/*
 * Copyright (c) 2010-2014:  Steinbeis Forschungszentrum (STZ Ölbronn)
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

#include "invocation.h"
#include "registry/class.h"
#include "type_converter.h"
#include <string>

#include "common/util/string_util.h"
#include "common/util/hash.h"

namespace ug {
namespace vrl {
namespace invocation {

//static ug::Hash<std::string, const ug::bridge::ClassNameNode*> classNameNodes;
static ug::Hash<std::string, const ug::bridge::IExportedClass*> classes;

void initClasses(ug::bridge::Registry &reg) {
	using namespace ug::bridge;

	classes = ug::Hash<std::string, const ug::bridge::IExportedClass*>(reg.num_classes()*2);
	classes.reserve(reg.num_classes());

	// only classes, no groups !
	for (unsigned int i = 0; i < reg.num_classes(); i++) {
		const ug::bridge::IExportedClass* c = &reg.get_class(i);
		classes.insert(c->name(), c);
	}
}

const ug::bridge::ExportedMethod* getMethodBySignature(
		JNIEnv *env,
		ug::bridge::Registry* reg,
		const ug::bridge::IExportedClass* clazz, bool readOnly,
		std::string methodName,
		jobjectArray params) {

	//	// create signature
	//	std::string signature = createMethodSignature(
	//			env, clazz->name(), methodName.c_str(), readOnly, params);
	//
	//	// search in map first and return result if entry exists
	//	if (methods.find(signature.c_str()) != methods.end()) {
	//		UG_LOG("FOUND:" << signature << std::endl);
	//		return methods[signature.c_str()];
	//	}

	const ug::bridge::ExportedMethod* method = nullptr;

	// we allow invocation of methods defined in parent classes
	std::vector<const ug::bridge::IExportedClass*> classList =
			getParentClasses(reg, clazz);

	// iterate over all classes of the inheritance path
	for (size_t i = 0; i < classList.size(); i++) {

		const ug::bridge::IExportedClass* cls = classList[i];

		// check whether to search const or non-const methods
		const ug::bridge::ExportedMethodGroup* methodGroup = nullptr;

		if (readOnly) {
			methodGroup = cls->get_const_exported_method_group(methodName);
		} else {
			methodGroup = cls->get_exported_method_group(methodName);
		}

		if(methodGroup == nullptr) continue;

		size_t numOverloads = methodGroup->num_overloads();

		// search without smart->raw conversion
		for(size_t k = 0; k < numOverloads; k++) {
			method = methodGroup->get_overload(k);
			if(compareParamTypes(env, params, reg, method->params_in(), false)) {
				// found correct overload
				return method;
			}
		} // end for k

		// if nothing found search with smart->raw conversion
                for(size_t k = 0; k < numOverloads; k++) {
			method = methodGroup->get_overload(k);
			if(compareParamTypes(env, params, reg, method->params_in(), true)) {
				// found correct overload
				return method;
			}
		} // end for k

	} // end for i

	return nullptr;
}

const ug::bridge::ExportedFunction* getFunctionBySignature(
		JNIEnv *env,
		ug::bridge::Registry* reg,
		std::string functionName,
		jobjectArray params) {

	//	// create signature
	//	std::string signature = createMethodSignature(
	//			env, "", "", false, params);
	//
	//	// search in map first and return result if entry exists
	//	if (methods.find(signature.c_str()) != methods.end()) {
	//		UG_LOG("FOUND:" << signature << std::endl);
	//		return functions[signature.c_str()];
	//	}

	const ug::bridge::ExportedFunction* func = nullptr;

	const ug::bridge::ExportedFunctionGroup *funcGroup = reg->get_exported_function_group(functionName);
	if(funcGroup == nullptr) return nullptr;
	size_t numOverloads = funcGroup->num_overloads();

	// search without smart->raw conversion
	for(size_t k = 0; k < numOverloads; k++) {
		func = funcGroup->get_overload(k);
		if (compareParamTypes(env, params, reg, func->params_in(), false)) {
			// found correct overload
			return func;
		}
	}

	// if nothing found search with smart->raw conversion
	for(size_t k = 0; k < numOverloads; k++) {
		func = funcGroup->get_overload(k);
		if (compareParamTypes(env, params, reg, func->params_in(), true)) {
			// found correct overload
			return func;
		}
	}
	// we did not find the correct overload
	// there is no other group having the same name
	return nullptr;
}

const ug::bridge::ExportedConstructor* getConstructorBySignature(
		JNIEnv *env,
		ug::bridge::Registry* reg,
		ug::bridge::IExportedClass* eCls,
		jobjectArray params) {

	size_t numConstructors = 0;

	numConstructors = eCls->num_constructors();

	// search without smart->raw conversion
	for (size_t i = 0; i < numConstructors; i++) {
		const ug::bridge::ExportedConstructor* constructor = 
		&eCls->get_constructor(i);

		// if the parameter types are equal
		// we found the correct constructor
		if (compareParamTypes(
				env, params, reg, constructor->params_in(), false)) {

			return constructor;
		}
	}

	// if nothing found search with smart->raw conversion
	for (size_t i = 0; i < numConstructors; i++) {
		const ug::bridge::ExportedConstructor* constructor = 
		&eCls->get_constructor(i);

		// if the parameter types are equal
		// we found the correct constructor
		if (compareParamTypes(
				env, params, reg, constructor->params_in(), true)) {

			return constructor;
		}
	}

	return nullptr;
}

const ug::bridge::IExportedClass* getExportedClassPtrByName(
		ug::bridge::Registry* reg,
		std::string className) {

	return classes.get_entry(className);
}

const ug::bridge::ClassNameNode* getClassNodePtrByName(
		ug::bridge::Registry* reg,
		std::string className) {

	if (className == "") {
		return nullptr;
	}

	return &classes.get_entry(className)->class_name_node();
}

//static std::map<const char*, const ug::bridge::ExportedMethod*> methods;
//static std::map<const char*, const ug::bridge::ExportedFunction*> functions;
//
//
//const std::string createMethodSignature(JNIEnv* env, const char* className,
//		const char* methodName, bool readOnly, jobjectArray params) {
//	std::stringstream signature;
//
//	signature << className << "::" << methodName << "::const=" << readOnly;
//	size_t size = env->GetArrayLength(params);
//
//	for (size_t i = 0; i < size; i++) {
//
//		jobject param = env->GetObjectArrayElement(params, i);
//
//		signature << "::" << paramClass2ParamType(env, param);
//	}
//
//	return signature.str();
//}

//const ug::bridge::ExportedMethod* getMethodBySignature(
//		JNIEnv *env,
//		ug::bridge::Registry* reg,
//		const ug::bridge::IExportedClass* clazz, bool readOnly,
//		std::string methodName,
//		jobjectArray params) {
//
//	//	// create signature
//	//	std::string signature = createMethodSignature(
//	//			env, clazz->name(), methodName.c_str(), readOnly, params);
//	//
//	//	// search in map first and return result if entry exists
//	//	if (methods.find(signature.c_str()) != methods.end()) {
//	//		UG_LOG("FOUND:" << signature << std::endl);
//	//		return methods[signature.c_str()];
//	//	}
//
//	// we allow invocation of methods defined in parent classes
//	std::vector<const ug::bridge::IExportedClass*> classList =
//			getParentClasses(reg, clazz);
//
//	// iterate over all classes of the inheritance path
//	for (unsigned i = 0; i < classList.size(); i++) {
//
//		const ug::bridge::IExportedClass* cls = classList[i];
//		unsigned int numMethods = 0;
//
//		// check whether to search const or non-const methods
//		if (readOnly) {
//			numMethods = cls->num_const_methods();
//		} else {
//			numMethods = cls->num_methods();
//		}
//
//		for (unsigned int j = 0; j < numMethods; j++) {
//			const ug::bridge::ExportedMethod* method = nullptr;
//
//			unsigned int numOverloads = 1;
//
//			if (readOnly) {
//				numOverloads = cls->num_const_overloads(j);
//			} else {
//				numOverloads = cls->num_overloads(j);
//			}
//
//			for (unsigned int k = 0; k < numOverloads; k++) {
//
//				// check whether to search const or non-const methods
//				if (readOnly) {
//					method = &cls->get_const_overload(j, k);
//				} else {
//					method = &cls->get_overload(j, k);
//				}
//
//				// if the method name and the parameter types are equal
//				// we found the correct method
//				if (strcmp(method->name().c_str(),
//						methodName.c_str()) == 0 &&
//						compareParamTypes(
//						env, params, reg, method->params_in())) {
//
//					//	// improve lookup time: add method signature to map
//					//	methods[signature.c_str()] = method;
//
//					return method;
//				}
//			} // for k
//		} // for j
//	} // for i
//
//	return nullptr;
//}

//const ug::bridge::ExportedFunction* getFunctionBySignature(
//		JNIEnv *env,
//		ug::bridge::Registry* reg,
//		std::string functionName,
//		jobjectArray params) {
//
//	//	// create signature
//	//	std::string signature = createMethodSignature(
//	//			env, "", "", false, params);
//	//
//	//	// search in map first and return result if entry exists
//	//	if (methods.find(signature.c_str()) != methods.end()) {
//	//		UG_LOG("FOUND:" << signature << std::endl);
//	//		return functions[signature.c_str()];
//	//	}
//
//	unsigned int numFunctions = 0;
//
//	numFunctions = reg->num_functions();
//
//	for (unsigned int i = 0; i < numFunctions; i++) {
//		const ug::bridge::ExportedFunction* func = nullptr;
//
//		unsigned int numOverloads = 1;
//
//		numOverloads = reg->num_overloads(i);
//
//		for (unsigned int k = 0; k < numOverloads; k++) {
//			func = &reg->get_overload(i, k);
//
//			// if the function name and the parameter types are equal
//			// we found the correct function
//			if (strcmp(func->name().c_str(),
//					functionName.c_str()) == 0 &&
//					compareParamTypes(
//					env, params, reg, func->params_in())) {
//
//				// improve lookup time
//				// functions[signature.c_str()] = func;
//				return func;
//			}
//		} // end for k
//	}// end for i
//
//	// no matching function found
//	return nullptr;
//}


} // invocation::
} // vrl::
} // ug::
