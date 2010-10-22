/* 
 * File:   type_converter.h
 * Author: miho
 *
 * Created on 5. Oktober 2010, 14:54
 */

#ifndef TYPE_CONVERTER_H
#define	TYPE_CONVERTER_H

#ifdef __GNUC__
#define PRETTY_FUNCTION __PRETTY_FUNCTION__
#else
#define PRETTY_FUNCTION "function name not available (not using GCC)"
#endif

#define VRL_DEBUG_MODE 1

// TODO: crappy, please improve me
#define VRL_DBG(msg, level) if( (VRL_DEBUG_MODE) >= (level)) std::cout << "UG4-VRL: in function \"" << PRETTY_FUNCTION << "\": " << msg << std::endl


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
		jstring stringC2J(JNIEnv *env, std::string const& s);

		/**
		 * <p>
		 * Converts a Java string to a native string.
		 * </p>
		 * <p>
		 * <b>Note:</p> this function must not be used to convert large amounts
		 * of strings because of inefficient memory handling.
		 * </p>
		 * @param env JVM environment to operate on
		 * @param s string to convert
		 * @return a native string
		 */
		std::string stringJ2C(JNIEnv *env, jstring const& s);

		/**
		 * Converts a native string array to a java object array.
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
		 * Converts a java string array to a native string array.
		 * </p>
		 * <p>
		 * <b>Note:</p> this function must not be used to convert large amounts
		 * of strings because of inefficient memory handling.
		 * </p>
		 * @param env JVM environment to operate on
		 * @param array array to convert
		 * @return a native string array
		 */
		std::vector<std::string> stringArrayJ2C(
				JNIEnv *env, jobjectArray const& array);

		inline bool boolJ2C(jboolean b) {
			return b == JNI_TRUE;
		}

		inline jboolean boolC2J(bool b) {
			return b ? JNI_TRUE : JNI_FALSE;
		}

		std::string exportedFunction2Groovy(
				ug::bridge::ExportedFunction const& func);

		std::string exportedClass2Groovy(ug::bridge::Registry* reg,
				ug::bridge::IExportedClass const& clazz);

		jobject boolean2JObject(JNIEnv *env, jboolean value);
		jboolean jObject2Boolean(JNIEnv *env, jobject obj);

		jobject int2JObject(JNIEnv *env, jint value);
		jint jObject2Int(JNIEnv *env, jobject obj);

		jobject double2JObject(JNIEnv *env, jdouble value);
		jdouble jObject2Double(JNIEnv *env, jobject obj);

		jobject string2JObject(JNIEnv *env, const char* value);
		std::string jObject2String(JNIEnv *env, jobject obj);

		jobject pointer2JObject(JNIEnv *env, void* value);
		void* jObject2Pointer(JNIEnv *env, jobject obj);


		std::string createParamInfo(const char* className,
				const std::vector<const char*>* classNames, bool isConst,
				std::string customInfo = "", std::string customOptions = "");


		std::string createMethodInfo(const char* className,
				const std::vector<const char*>* classNames, bool isConst,
				std::string customInfo = "", std::string customOptions = "");


		std::string paramType2String(int paramType,
				const char* className,
				const std::vector<const char*>* classNames, std::string paramOptions, bool isOutput = false);

		int jobjectArray2ParamStack(JNIEnv *env, ug::bridge::ParameterStack& paramsOut,
				const ug::bridge::ParameterStack& paramsTemplate,
				jobjectArray const& array);

		jobject param2JObject(JNIEnv *env, ug::bridge::ParameterStack& params, int index);

		jobject getClass(JNIEnv *env, jobject obj);
		std::string getClassName(JNIEnv *env, jobject obj);

		uint paramClass2ParamType(JNIEnv *env, jobject obj);

		bool compareParamTypes(JNIEnv *env, jobjectArray params,
				ug::bridge::ParameterStack& paramStack);

		const ug::bridge::ExportedMethod* getMethodBySignature(
				JNIEnv *env,
				ug::bridge::Registry* reg,
				const ug::bridge::IExportedClass* clazz,
				bool readOnly,
				std::string methodName,
				jobjectArray params);

		const ug::bridge::IExportedClass* getExportedClassPtrByName(
				ug::bridge::Registry* reg,
				std::string className);

		const std::vector<const ug::bridge::IExportedClass*> getParentClasses(
				ug::bridge::Registry* reg, const ug::bridge::IExportedClass* clazz);

		void generateMethods(std::stringstream& result,
				ug::bridge::IExportedClass* clazz);

//		jobject messageTypeC2J(JNIEnv *env, MessageType type);

	} // end vrl::
}// end ug::

#endif	/* TYPE_CONVERTER_H */

