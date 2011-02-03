/* 
 * File:   type_converter.h
 * Author: Michael Hoffer <info@michaelhoffer.de>
 *
 * Created on 5. Oktober 2010, 14:54
 */

#ifndef TYPE_CONVERTER_H
#define	TYPE_CONVERTER_H

#include<jni.h>
#include<string>
#include<vector>
#include "ug_bridge/class.h"
#include "ug_bridge/registry.h"
#include "messaging.h"

namespace ug {
	namespace vrl {

		enum InvokationType {
			Function,
			METHOD
		};

		/**
		 * Converts a native string to a Java string.
		 * @param env JVM environment to operate on
		 * @param s string to convert
		 * @return a java string
		 */
		jstring stringC2J(JNIEnv *env, const char* s);

		//		jstring stringC2J(JNIEnv *env, std::string const& s);

		/**
		 * <p>
		 * Converts a Java string to a native string.
		 * </p>
		 * <p>
		 * <b>Note:</b> this function must not be used to convert large amounts
		 * of strings because of inefficient memory handling.
		 * </p>
		 * @param env JVM environment to operate on
		 * @param s string to convert
		 * @return a native string
		 */
		std::string stringJ2C(JNIEnv *env, jstring const& s);

		/**
		 * Converts a native string array to a Java object array.
		 * @param env JVM environment to operate on
		 * @param strings array to convert
		 * @return a java object array
		 */
		jobjectArray stringArrayC2J(JNIEnv *env, const std::string* strings);

		/**
		 * Converts a native string array to a java object array.
		 * @param env JVM environment to operate on
		 * @param strings array to convert
		 * @return a java object array
		 */
		jobjectArray stringArrayC2J(
				JNIEnv *env,
				std::vector<std::string> const& strings);

		/**
		 * <p>
		 * Converts a Java string array to a native string array.
		 * </p>
		 * <p>
		 * <b>Note:</b> this function must not be used to convert large amounts
		 * of strings because of inefficient memory handling.
		 * </p>
		 * @param env JVM environment to operate on
		 * @param array array to convert
		 * @return a native string array
		 */
		std::vector<std::string> stringArrayJ2C(
				JNIEnv *env, jobjectArray const& array);

		/**
		 * Converts a Java boolean to a native boolean.
		 * @param b Java boolean to convert
		 * @return native boolean
		 */
		inline bool boolJ2C(jboolean b) {
			return b == JNI_TRUE;
		}

		/**
		 * Converts a native boolean to a Java boolean.
		 * @param b native boolean to convert
		 * @return Java boolean
		 */
		inline jboolean boolC2J(bool b) {
			return b ? JNI_TRUE : JNI_FALSE;
		}

		/**
		 * Converts an exported function to Groovy source code.
		 * @param func function to convert
		 * @return Groovy source code
		 */
		std::string exportedFunction2Groovy(
				ug::bridge::ExportedFunction const& func);

		/**
		 * Converts an exported class to Groovy source code.
		 * @param reg registry the class is registered with
		 * @param clazz class to convert
		 * @return Groovy source code
		 */
		std::string exportedClass2Groovy(ug::bridge::Registry* reg,
				ug::bridge::IExportedClass const& clazz);

		/**
		 * Converts a jboolean to a Java object (jobject).
		 * @param env JVM environment to operate on
		 * @param value value to convert
		 * @return Java object (jobject)
		 */
		jobject boolean2JObject(JNIEnv *env, jboolean value);

		/**
		 * Converts a Java object to a native boolean value.
		 * @param env JVM environment to operate on
		 * @param obj object to convert
		 * @return a native boolean value
		 */
		jboolean jObject2Boolean(JNIEnv *env, jobject obj);

		/**
		 * Converts an integer to a Java object (jobject).
		 * @param env JVM environment to operate on
		 * @param value value to convert
		 * @return Java object (jobject)
		 */
		jobject int2JObject(JNIEnv *env, jint value);

		/**
		 * Converts a Java object (jobject) to a native int value.
		 * @param env JVM environment to operate on
		 * @param obj object to convert
		 * @return a native int value
		 */
		jint jObject2Int(JNIEnv *env, jobject obj);

		/**
		 * Converts a double to a Java object (jobject).
		 * @param env JVM environment to operate on
		 * @param value value to convert
		 * @return a Java object (jobject)
		 */
		jobject double2JObject(JNIEnv *env, jdouble value);

		/**
		 * Converts a Java object (jobject) to a native double value.
		 * @param env JVM environment to operate on
		 * @param obj object to convert
		 * @return a native double value
		 */
		jdouble jObject2Double(JNIEnv *env, jobject obj);

		/**
		 * Converts a string to a Java object (jobject).
		 * @param env JVM environment to operate on
		 * @param value value to convert
		 * @return a Java object (jobject)
		 */
		jobject string2JObject(JNIEnv *env, const char* value);

		/**
		 * Converts a Java object (jobject) to a native string.
		 * @param env JVM environment to operate on
		 * @param obj object to convert
		 * @return a native string
		 */
		std::string jObject2String(JNIEnv *env, jobject obj);

		/**
		 * Converts a native pointer to a Java object (jobject).
		 * @param env JVM environment to operate on
		 * @param value pointer to convert
		 * @return a Java object (jobject)
		 */
		jobject pointer2JObject(JNIEnv *env, void* value);

		/**
		 * Converts a Java object (jobject) to a native pointer.
		 * @param env JVM environment to operate on
		 * @param obj object to convert
		 * @return a native pointer
		 */
		void* jObject2Pointer(JNIEnv *env, jobject obj);


		/**
		 * Creates a VRL param info to customize type representation and to
		 * ensure type-safe connections.
//		 * @param paramName name of the parameter
		 * @param className name of the param class
		 * @param classNames own class name plus class names of all super
		 *        classes
		 * @param isConst defines whether this parameter shall be const
		 * @param paramInfo additional param info options
		 * @param customParamInfo additional value options
		 * @return a string containing a VRL param info (Groovy source code)
		 */
		std::string createParamInfo(const char* paramName, const char* className,
				const std::vector<const char*>* classNames, bool isConst,
				std::vector<std::string> const& paramInfo,
				std::string const& customParamInfo="");


		/**
		 * Creates a VRL method info to customize method representation.
		 * @param className name of the return value class
		 * @param classNames return value class name plus class names of all
		 *        super classes
		 * @param isConst defines whether the return value shall be const
		 * @param customInfo additional method info options
		 * @param customOptions additional return value options
		 * @return a string containing a VRL method info (Groovy source code)
		 */
		std::string createMethodInfo(const char* className,
				const std::vector<const char*>* classNames, bool isConst,
				std::string customInfo = "", std::string customOptions = "");


		/**
		 * Converts a ug::bridge::ParameterTypes value to the equivalent Java
		 * type name.
		 * @param paramName name of the parameter
		 * @param paramType param type to convert
		 * @param className name of the value class
		 * @param classNames value class name plus class names of all super
		 *        classes
		 * @param paramInfo param info
		 * @param isOutput defines whether this value is a return value
		 *                 (in this case no param info will be added)
		 * @return a String containing the Java type name
		 */
		std::string paramType2String(int paramType, const char* paramName,
				const char* className,
				const std::vector<const char*>* classNames,
				std::vector<std::string> const& paramInfo, bool isOutput = false);

		/**
		 * Converts an array of Java objects to a parameter stack.
		 * @param env JVM environment to operate on
		 * @param paramsOut converted parameter stack (return value)
		 * @param paramsTemplate template parameter stack used to get correct
		 *                       parameter type
		 * @param array object array to convert
		 */
		void jobjectArray2ParamStack(JNIEnv *env,
				ug::bridge::ParameterStack& paramsOut,
				const ug::bridge::ParameterStack& paramsTemplate,
				jobjectArray const& array);

		/**
		 * Creates an empty Java array using the specified class as element
		 * type.
		 * @param env JVM environment to operate on
		 * @param className name of the element class
		 * @return emtpy Java array
		 */
		jobjectArray createEmptyJavaArray(
				JNIEnv *env, std::string className);

		/**
		 * Creates an empty Java array using the specified class as element
		 * type.
		 * @param env JVM environment to operate on
		 * @param elementClass element class
		 * @return emtpy Java array
		 */
		jobjectArray createEmptyJavaArray(
				JNIEnv *env, jclass elementClass);

		/**
		 * Converts a parameter stack entry to a Java object.
		 * @param env JVM environment to operate on
		 * @param params parameter stack to convert
		 * @param index index of the element to convert
		 * @return a Java object (jobject)
		 */
		jobject param2JObject(JNIEnv *env,
				ug::bridge::ParameterStack& params, int index);

		/**
		 * Returns the class object of the specified Java object.
		 * @param env JVM environment to operate on
		 * @param obj Java object
		 * @return class object of the specified Java object
		 */
		jobject getClass(JNIEnv *env, jobject obj);

		/**
		 * Returns the class name of the specified Java object.
		 * @param env JVM environment to operate on
		 * @param obj Java object
		 * @return class name of the specified Java object
		 */
		std::string getClassName(JNIEnv *env, jobject obj);

		/**
		 * Returns the parameter type (ug::bridge::ParameterTypes) of the
		 * specified Java object.
		 * @param env JVM environment to operate on
		 * @param obj Java object
		 * @return parameter type (ug::bridge::ParameterTypes) of the
		 *         specified Java object
		 */
		uint paramClass2ParamType(JNIEnv *env, jobject obj);

		/**
		 * Compares the parameter types of a Java object array and a
		 * paramter stack
		 * @param env JVM environment to operate on
		 * @param params array of Java objects
		 * @param paramStack parameter stack
		 * @return <code>true</code> if parameter types are equal;
		 *         <code>false</code> otherwise
		 */
		bool compareParamTypes(JNIEnv *env, jobjectArray params,
				const ug::bridge::ParameterStack& paramStack);

		/**
		 * Returns parent classes (super classes) of an exported class.
		 * @param reg registry to search
		 * @param clazz exported class
		 * @return a vector containing all parent classes of the given class
		 */
		const std::vector<const ug::bridge::IExportedClass*> getParentClasses(
				ug::bridge::Registry* reg,
				const ug::bridge::IExportedClass* clazz);

		/**
		 * Generates Groovy source code for all non const methods of the
		 * specified class.
		 * @param result stream to use (appends return value)
		 * @param clazz exported class
		 */
		void generateMethods(std::stringstream& result,
				ug::bridge::IExportedClass* clazz);

		/**
		 * Generates Groovy source code for all const methods of the specified
		 * class.
		 * @param result string stream to use (appends return value)
		 * @param clazz exported class
		 */
		void generateConstMethods(std::stringstream& result,
				ug::bridge::IExportedClass* clazz);

		/**
		 * Generates the method header to the specified method as Groovy
		 * source code. The header contains the return value, the method name
		 * and the parameter list.
		 * Parameters are stored as array of java objects
		 * (variable name is <code>params</code>).
		 * @param result string stream to use (appends return value)
		 * @param method method to generate
		 * @param isFunction defines whether to generate a function header
		 * @param isVisual defines whether the method shall be visualized
		 * @param prefix method prefix (optional)
		 */
		void generateMethodHeader(
				std::stringstream& result,
				ug::bridge::ExportedFunctionBase const& method,
				bool isFunction = false, bool isVisual = true, std::string prefix = "");

	} // end vrl::
}// end ug::

#endif	/* TYPE_CONVERTER_H */

