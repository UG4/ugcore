#include "invocation.h"
#include "registry/class.h"
#include "type_converter.h"
#include <string>

#include "common/util/string_util.h"
#include "common/util/hash.h"

namespace ug {
namespace vrl {
namespace invocation {

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

	const ug::bridge::ExportedMethod* method = NULL;

	// we allow invocation of methods defined in parent classes
	std::vector<const ug::bridge::IExportedClass*> classList =
			getParentClasses(reg, clazz);

	// iterate over all classes of the inheritance path
	for (size_t i = 0; i < classList.size(); i++) {

		const ug::bridge::IExportedClass* cls = classList[i];

		// check whether to search const or non-const methods
		const ug::bridge::ExportedMethodGroup* methodGroup = NULL;

		if (readOnly) {
			methodGroup = cls->get_const_exported_method_group(methodName);
		} else {
			methodGroup = cls->get_exported_method_group(methodName);
		}

		if(methodGroup == NULL) continue;

		size_t numOverloads = methodGroup->num_overloads();

		for(size_t k = 0; k < numOverloads; k++) {
			method = methodGroup->get_overload(k);
			if(compareParamTypes(env, params, reg, method->params_in())) {
				// found correct overload
				return method;
			}
		} // end for k
	} // end for i

	return NULL;
}
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
//			const ug::bridge::ExportedMethod* method = NULL;
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
//	return NULL;
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
//		const ug::bridge::ExportedFunction* func = NULL;
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
//	return NULL;
//}

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

	const ug::bridge::ExportedFunction* func = NULL;

	const ug::bridge::ExportedFunctionGroup *funcGroup = reg->get_exported_function_group(functionName);
	if(funcGroup == NULL) return NULL;
	size_t numOverloads = funcGroup->num_overloads();
	for(size_t k = 0; k < numOverloads; k++) {
		func = funcGroup->get_overload(k);
		if (compareParamTypes(env, params, reg, func->params_in())) {
			// found correct overload
			return func;
		}
	}
	// we did not find the correct overload
	// there is no other group having the same name
	return NULL;
}

const ug::bridge::ExportedConstructor* getConstructorBySignature(
		JNIEnv *env,
		ug::bridge::Registry* reg,
		ug::bridge::IExportedClass* eCls,
		jobjectArray params) {

	size_t numConstructors = 0;

	numConstructors = eCls->num_constructors();

	for (size_t i = 0; i < numConstructors; i++) {
		const ug::bridge::ExportedConstructor* constructor = 
		&eCls->get_constructor(i);

		// if the parameter types are equal
		// we found the correct constructor
		if (compareParamTypes(
				env, params, reg, constructor->params_in())) {

			return constructor;
		}
	}

	return NULL;
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
		return NULL;
	}

	return &classes.get_entry(className)->class_name_node();
}

} // invocation::
} // vrl::
} // ug::
