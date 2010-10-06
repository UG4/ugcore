/* 
 * File:   type_converter.h
 * Author: miho
 *
 * Created on 5. Oktober 2010, 14:54
 */

#ifndef TYPE_CONVERTER_H
#define	TYPE_CONVERTER_H

#include<jni.h>
#include<string>
#include<vector>
#include "ug_bridge/class.h"

namespace ug {
	namespace vrl {
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

		std::string exportedFunction2Groovy(
				ug::bridge::ExportedFunction const& func);

		jobject int2JObject(JNIEnv *env, jint value);
		jint jObject2Int(JNIEnv *env, jobject obj);

		const char* paramType2String(int paramType);

		int jobjectArray2ParamStack(JNIEnv *env, ug::bridge::ParameterStack& paramsOut,
					const ug::bridge::ParameterStack& paramsTemplate,
					jobjectArray const& array);

		jobject param2JObject(JNIEnv *env, ug::bridge::ParameterStack& params, int index);

	} // end vrl::
}// end ug::

#endif	/* TYPE_CONVERTER_H */

