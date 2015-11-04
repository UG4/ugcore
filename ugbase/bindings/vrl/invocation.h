
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

#endif	/* INVOKATION_H */

