#include "type_converter.h"
#include "ug_bridge/registry.h"
#include "messaging.h"
#include "invocation.h"
#include <sstream>

namespace ug {
namespace vrl {

//		jstring stringC2J(JNIEnv *env, std::string const& s) {
//			return stringC2J(env, s.c_str());
//		}

jstring stringC2J(JNIEnv *env, const char* s) {
	return env->NewStringUTF(s);
}

std::string stringJ2C(JNIEnv *env, jstring const& s) {
	const char* tmpStr = env->GetStringUTFChars(s, false);
	std::string result = (std::string)tmpStr;
	env->ReleaseStringUTFChars(s, tmpStr);
	return result;
}

jobjectArray stringArrayC2J(
		JNIEnv *env,
		const std::string* strings,
		const unsigned int array_length) {
	jclass stringClass = env->FindClass("java/lang/String");

	jobjectArray result =
			env->NewObjectArray(array_length, stringClass, 0);

	// convert array elements
	for (unsigned int i = 0; i < array_length; i++) {
		std::string s = strings[ i ];
		jstring javaString = env->NewStringUTF(s.c_str());
		env->SetObjectArrayElement(result, i, javaString);
	}

	return result;
}

jobjectArray createEmptyJavaArray(
		JNIEnv *env, std::string className) {
	jclass elementClass = env->FindClass(className.c_str());

	return createEmptyJavaArray(env, elementClass);
}

jobjectArray createEmptyJavaArray(
		JNIEnv *env, jclass elementClass) {

	return env->NewObjectArray(0, elementClass, 0);
}

jobjectArray stringArrayC2J(
		JNIEnv *env, std::vector<std::string> const& strings) {

	if (strings.size() > 0) {
		// it is safe to give a pointer to the first vector element as
		// std::vector implementation uses contiguous memory
		return stringArrayC2J(env, &strings[0], strings.size());
	} else {
		// create an empty string array
		jclass cls = env->FindClass("java/lang/String");
		return env->NewObjectArray(0, cls, 0);
	}
}

std::vector<std::string> stringArrayJ2C(
		JNIEnv *env, jobjectArray const& array) {

	std::vector<std::string> result;

	unsigned int length = env->GetArrayLength(array);

	// convert each element of the java object array to a std string
	// and add it to the result vector
	for (unsigned int i = 0; i < length; i++) {
		result.push_back(stringJ2C(env,
				(jstring) env->GetObjectArrayElement(array, i)));
	}

	return result;
}

void generateMethodHeader(
		std::stringstream& result,
		ug::bridge::ExportedFunctionBase const& method,
		bool isFunction, std::string prefix) {

	std::stringstream methodHeaderParams;
	std::stringstream paramArrayForInvokation;

	std::string name = method.name();

	const ug::bridge::ParameterStack& paramStackIn = method.params_in();
	const ug::bridge::ParameterStack& paramStackOut = method.params_out();

	size_t numParams = paramStackIn.size();

	// generate method paramters for method header
	// (including VRL param info) and for method invokation
	for (unsigned int i = 0; i < numParams; i++) {
		if (i > 0) {
			methodHeaderParams << ",\n ";
			paramArrayForInvokation << ", ";
		}

		//				UG_LOG("BEFORE param2string:in loop: \n");

		methodHeaderParams << paramType2String(
				paramStackIn.get_type(i),
				method.parameter_name(i).c_str(),
				paramStackIn.class_name(i),
				paramStackIn.class_names(i),
				method.parameter_info_vec(i)) << " p" << i;

		paramArrayForInvokation << " p" << i;
	}

	if (!isFunction) {
		// we always need the visual id to get a reference to the
		// visualization that invokes this method,
		// that is why we add a visual id request to the param list
		if (numParams > 0) {
			methodHeaderParams << ", ";
			paramArrayForInvokation << ", ";
		}
		methodHeaderParams << " VisualIDRequest id ";
		paramArrayForInvokation << " id";
	}

	//			UG_LOG("RETURN_VEC:" << method.return_info_vec().size() << std::endl);

	bool readOnly = false;

	std::string className = "";
	const std::vector<const char*>* classNames = NULL;

	// return value generation
	std::string outType;
	if (paramStackOut.size() > 0) {
		outType = paramType2String(
				paramStackOut.get_type(0),
				paramStackOut.class_name(0), // param name
				paramStackOut.class_name(0),
				paramStackOut.class_names(0),
				method.return_info_vec(), true);

		className = paramStackOut.class_name(0);
		classNames = paramStackOut.class_names(0);

		readOnly =
				paramStackOut.get_type(0) == ug::bridge::PT_CONST_POINTER;
	} else {
		outType = "void";
	}

	// generate method info including return value info
	// (equivalent to param info)
	result << createMethodInfo(
			className.c_str(),
			classNames,
			readOnly, method.options()) << "\n";

	// putting it all together
	result << "public " << outType << " " << prefix << name << " ("
			<< methodHeaderParams.str() << ") {\n"
			<< "Object[] params = ["
			<< paramArrayForInvokation.str() << "] \n";

	// If we are not generating a function it is necessary to
	// call the updatePointer method. Its purpose is to visually invoke
	// the setPointer() and getPointer() method before this method is
	// invoked. Otherwise we might operate on wrong instance
	if (!isFunction) {
		result << "updatePointer(id);\n";
	}

	if (paramStackOut.size() > 0) {
		result << "return ";
	}
}

std::string exportedFunction2Groovy(
		ug::bridge::ExportedFunction const& func) {
	std::stringstream result;

	std::string group = func.group();

	// create component info that specifies the menu group this
	// function shall be added to
	result << "@ComponentInfo(name=\"" << func.name()
			<< "\", category=\"" << group << "\")\n"
			<< "public class UG4_" << func.name()
			<< " implements Serializable {\n"
			<< "private static final long serialVersionUID=1L;\n";

	// function generation
	generateMethodHeader(
			result, func, true);
	result << "edu.gcsc.vrl.ug4.UG4.getUG4().invokeFunction( \""
			<< func.name() << "\", false, params)";

	result << "\n}\n}";

	return result.str();
}

void generateMethods(std::stringstream& result,
		const ug::bridge::IExportedClass* clazz) {
	for (unsigned int i = 0; i < clazz->num_methods(); i++) {
		const ug::bridge::ExportedMethod &method = clazz->get_method(i);

		generateMethodHeader(result, method);

		result << "edu.gcsc.vrl.ug4.UG4.getUG4().invokeMethod("
				<< "getClassName(),"
				<< " getPointer().getAddress(), false, \""
				<< method.name() << "\", params)";

		result << "\n}\n\n";
	}
}

void generateConstMethods(std::stringstream& result,
		const ug::bridge::IExportedClass* clazz) {
	for (unsigned int i = 0; i < clazz->num_const_methods(); i++) {
		const ug::bridge::ExportedMethod &method = clazz->get_const_method(i);

		generateMethodHeader(result, method, false, "const_");

		result << "edu.gcsc.vrl.ug4.UG4.getUG4().invokeMethod("
				<< "getClassName(),"
				<< " getPointer().getAddress(), true, \""
				<< method.name() << "\", params)";

		result << "\n}\n\n";
	}
}

std::string exportedClass2Groovy(ug::bridge::Registry* reg,
		ug::bridge::IExportedClass const& clazz) {
	std::stringstream result;

	std::string className = "UG4_" + std::string(clazz.name());
	std::string group = clazz.group();

	// create component info that specifies the menu group this
	// class shall be added to
	result << "@ComponentInfo(name=\"" << clazz.name()
			<< "\", category=\"" << group << "\")\n"
			<< "public class " << className
			<< " extends edu.gcsc.vrl.ug4.UGObject {\n"
			<< "private static final long serialVersionUID=1L;\n";

	result << "public " << className << "() {\n"
			<< "setClassName(\"" << clazz.name() << "\");\n}\n";

	// gather inheritance information (necessary for type-safety)
	std::vector<const ug::bridge::IExportedClass*> baseClasses =
			getParentClasses(reg, &clazz);
	// generate method implementations
	for (unsigned int i = 0; i < baseClasses.size(); i++) {
		generateMethods(result, baseClasses[i]);
		generateConstMethods(result, baseClasses[i]);
	}
	// generate methods for this pointer handling
	result << createMethodInfo(clazz.name(), clazz.class_names(),
			false, "interactive=false")
			<< "\nedu.gcsc.vrl.ug4.Pointer getPointer()"
			<< " { super.getPointer()}\n";

	std::vector<std::string> paramInfo;

	paramInfo.push_back("");
	paramInfo.push_back(""); // put nullIsValid in this option field
	paramInfo.push_back("");

	result << "@MethodInfo(interactive=false,hide=true)\n"
			<< " void setPointer("
			<< createParamInfo(clazz.name(), clazz.name(), clazz.class_names(),
			false, paramInfo, "nullIsValid=true")
			<< " edu.gcsc.vrl.ug4.Pointer p) { super.setPointer(p)}\n";

	result << "\n}";

	return result.str();
}

const std::vector<const ug::bridge::IExportedClass*> getParentClasses(
		ug::bridge::Registry* reg,
		const ug::bridge::IExportedClass* clazz) {
	std::vector<const ug::bridge::IExportedClass*> result;
	// search registered classes by name as specified in the
	// class_names vector and add them to the result vector
	for (unsigned int i = 0; i < clazz->class_names()->size(); i++) {
		const ug::bridge::IExportedClass* baseCls =
				ug::vrl::invocation::getExportedClassPtrByName(
				reg, (*clazz->class_names())[i]);
		if (baseCls != NULL) {
			result.push_back(baseCls);
		}
	}

	return result;
}

jobject int2JObject(JNIEnv *env, jint value) {
	jclass cls = env->FindClass("java/lang/Integer");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "(I)V");
	return env->NewObject(cls, methodID, value);
}

jint jObject2Int(JNIEnv *env, jobject obj) {
	jclass argClass = env->GetObjectClass(obj);
	jmethodID ajf = env->GetMethodID(argClass, "intValue", "()I");
	return env->CallIntMethod(obj, ajf);
}

jobject double2JObject(JNIEnv *env, jdouble value) {
	jclass cls = env->FindClass("java/lang/Double");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "(D)V");
	return env->NewObject(cls, methodID, value);
}

jdouble jObject2Double(JNIEnv *env, jobject obj) {
	jclass argClass = env->GetObjectClass(obj);
	jmethodID ajf = env->GetMethodID(argClass, "doubleValue", "()D");
	return env->CallDoubleMethod(obj, ajf);
}

jobject boolean2JObject(JNIEnv *env, jboolean value) {
	jclass cls = env->FindClass("java/lang/Boolean");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "(Z)V");
	return env->NewObject(cls, methodID, value);
}

jboolean jObject2Boolean(JNIEnv *env, jobject obj) {
	jclass argClass = env->GetObjectClass(obj);
	jmethodID ajf = env->GetMethodID(argClass, "booleanValue", "()Z");
	return env->CallBooleanMethod(obj, ajf);
}

void* jObject2Pointer(JNIEnv *env, jobject obj) {
	jclass argClass = env->GetObjectClass(obj);
	jmethodID ajf = env->GetMethodID(argClass, "getAddress", "()J");
	return (void*) env->CallLongMethod(obj, ajf);
}

jobject pointer2JObject(JNIEnv *env, void* value) {
	jclass cls = env->FindClass("edu/gcsc/vrl/ug4/Pointer");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "(J)V");
	return env->NewObject(cls, methodID, (jlong) value);
}

jobject string2JObject(JNIEnv *env, const char* value) {
	return env->NewStringUTF(value);
}

std::string jObject2String(JNIEnv *env, jobject obj) {
	return stringJ2C(env, (jstring) obj);
}

std::string createParamInfo(const char* paramName, const char* className,
		const std::vector<const char*>* classNames, bool isConst,
		std::vector<std::string> const& paramInfo,
		std::string const& additionalParamInfo) {

	std::string customInfo = paramInfo.at(1);
	// add escape layer to simplify syntax
	std::string customOptions = replaceAll(paramInfo.at(2), "\"", "\\\"");

	std::stringstream paramInfoStream;
	std::stringstream classNameOptionsStream;

	// if class name or class names strings are empty no param info can
	// be generated, return error message as Groovy comment instead
	if (className == NULL) {
		//return std::string("/*ERROR PARAMINFO CLASSNAME == NULL*/");
		className = "";
	}

	if (classNames == NULL) {
		//return std::string("/*ERROR PARAMINFO CLASSNAMES == NULL*/");
		classNames = new std::vector<const char*>;
	}

	// creating VRL param options to ensure type-safety
	classNameOptionsStream
			<< ", options=\"className=\\\"" << className << "\\\";"
			<< "classNames=[";

	for (unsigned int i = 0; i < classNames->size(); i++) {
		if (i > 0) {
			classNameOptionsStream << ",";
		}
		classNameOptionsStream << "\\\"" << (*classNames)[i] << "\\\"";
	}

	classNameOptionsStream << "]";

	if (isConst) {
		classNameOptionsStream << "; readOnly=true";
	} else {
		classNameOptionsStream << "; readOnly=false";
	}

	// putting it all together
	paramInfoStream
			<< "@ParamInfo( name=\""
			<< paramName << "\""
			<< classNameOptionsStream.str()
			<< "; " << customOptions << "\"";

	if (customInfo.size() > 0) {
		paramInfoStream << ", style=\"" << customInfo << "\"";

	}

	if (additionalParamInfo.size() > 0) {
		paramInfoStream << ", " << additionalParamInfo;

	}

	paramInfoStream << ") ";

	return paramInfoStream.str();
}

std::string createMethodInfo(const char* className,
		const std::vector<const char*>* classNames, bool isConst,
		std::string customInfo, std::string customOptions) {
	std::stringstream methodInfo;
	std::stringstream classNameOptions;
	//
	if (className == NULL) {
		//return std::string("/*ERROR METHODINFO CLASSNAME == NULL*/");
		className = "";
	}

	if (classNames == NULL) {
		//return std::string("/*ERROR METHODINFO CLASSNAMES == NULL*/");
		classNames = new std::vector<const char*>;
	}

	// creating VRL param options to ensure type-safety for the return
	// value
	classNameOptions
			<< ", valueOptions=\"className=\\\"" << className << "\\\";"
			<< "classNames=[";

	for (unsigned int i = 0; i < classNames->size(); i++) {
		if (i > 0) {
			classNameOptions << ",";
		}

		classNameOptions << "\\\"" << (*classNames)[i] << "\\\"";
	}

	classNameOptions << "]";

	if (isConst) {
		classNameOptions << "; readOnly=true";
	} else {
		classNameOptions << "; readOnly=false";
	}

	// putting it all together
	methodInfo
			<< "@MethodInfo( valueName=\""
			<< className << "\""
			<< classNameOptions.str() << "; " << customOptions << "\"";

	if (customInfo.size() > 0) {
		methodInfo << ", " << customInfo;
	}

	methodInfo << ") ";

	return methodInfo.str();
}

std::string paramType2String(int paramType, const char* paramName,
		const char* className,
		const std::vector<const char*>* classNames,
		std::vector<std::string> const& paramInfo, bool isOutput) {

	switch (paramType) {
		case ug::bridge::PT_BOOL:
		{
			if (isOutput) {
				return "boolean";
			} else {
				std::string result =
						createParamInfo(paramName, className,
						classNames, false, paramInfo) +
						std::string("boolean");

				return result.c_str();
			}
		}
		case ug::bridge::PT_INTEGER:
		{
			if (isOutput) {
				return "int";
			} else {
				std::string result =
						createParamInfo(paramName, className,
						classNames, false, paramInfo) +
						std::string("int");

				return result.c_str();
			}
		}
		case ug::bridge::PT_NUMBER:
		{
			if (isOutput) {
				return "double";
			} else {
				std::string result =
						createParamInfo(paramName, className,
						classNames, false, paramInfo) +
						std::string("double");

				return result.c_str();
			}
		}
		case ug::bridge::PT_STRING:
		{
			if (isOutput) {
				return "String";
			} else {
				std::string result =
						createParamInfo(paramName, "",
						new std::vector<const char*>(),
						false, paramInfo) +
						std::string("String");

				return result.c_str();
			}
		}
		case ug::bridge::PT_POINTER:
		{
			if (isOutput) {
				return "edu.gcsc.vrl.ug4.Pointer";
			} else {
				std::string result =
						createParamInfo(paramName,
						className, classNames, false, paramInfo) +
						std::string("edu.gcsc.vrl.ug4.Pointer");

				return result.c_str();
			}
		}
		case ug::bridge::PT_CONST_POINTER:
		{
			if (isOutput) {
				return "edu.gcsc.vrl.ug4.Pointer";
			} else {
				std::string result =
						createParamInfo(paramName,
						className, classNames, true, paramInfo) +
						std::string("edu.gcsc.vrl.ug4.Pointer");

				return result.c_str();
			}
		}
		default: return "Object";
	}
}

jobject getClass(JNIEnv *env, jobject obj) {
	jclass classMethodAccess = env->FindClass("java/lang/Class");

	jmethodID classNameMethodID = env->GetMethodID(classMethodAccess,
			"getClass", "()Ljava/lang/Class;");

	return env->CallObjectMethod(obj,
			classNameMethodID);
}

std::string getClassName(JNIEnv *env, jobject obj) {
	jclass classMethodAccess = env->FindClass("java/lang/Class");

	jmethodID classNameMethodID = env->GetMethodID(classMethodAccess,
			"getClass", "()Ljava/lang/Class;");

	jobject clazz = env->CallObjectMethod(obj,
			classNameMethodID);

	classNameMethodID = env->GetMethodID(classMethodAccess,
			"getName", "()Ljava/lang/String;");
	jobject resultObj = env->CallObjectMethod(clazz,
			classNameMethodID);

	return jObject2String(env, resultObj);
}

uint paramClass2ParamType(JNIEnv *env, jobject obj) {
	int result = ug::bridge::PT_UNKNOWN;

	std::string className = getClassName(env, obj);

	if (className.compare("java.lang.Boolean") == 0) {
		result = ug::bridge::PT_BOOL;
	} else if (className.compare("java.lang.Integer") == 0) {
		result = ug::bridge::PT_INTEGER;
	} else if (className.compare("java.lang.Double") == 0) {
		result = ug::bridge::PT_NUMBER;
	} else if (className.compare("java.lang.String") == 0) {
		result = ug::bridge::PT_STRING;
	} else if (className.compare("edu.gcsc.vrl.ug4.Pointer") == 0) {
		result = ug::bridge::PT_POINTER;
	}
	if (className.compare("edu.gcsc.vrl.ug4.Pointer") == 0) {
		result = ug::bridge::PT_POINTER;
	}
	// What about const pointer?
	// Answer: compare param types allows
	// non-const* to const* conversion
	// That is why we do not check that. We also use the same class for
	// const and non const pointer. Const checking is done via readOnly
	// bit of the pointer instance (Java wrapper).

	return result;
}

bool compareParamTypes(JNIEnv *env, jobjectArray params,
		ug::bridge::ParameterStack const& paramStack) {

	// iterate over all param stack elements and compare their type with
	// the corresponding elements in the specified Java array
	for (unsigned int i = 0; i < (unsigned int) paramStack.size(); i++) {
		jobject param = env->GetObjectArrayElement(params, i);
		uint paramType = paramClass2ParamType(env, param);

		// allow non-const * to const *
		if (paramType == ug::bridge::PT_POINTER &&
				paramStack.get_type(i) == ug::bridge::PT_CONST_POINTER) {
			paramType = ug::bridge::PT_CONST_POINTER;
		}

		if (paramType != paramStack.get_type(i)) {
			return false;
		}
	}

	return true;
}

void jobjectArray2ParamStack(
		JNIEnv *env, ug::bridge::ParameterStack& paramsOut,
		const ug::bridge::ParameterStack& paramsTemplate,
		jobjectArray const& array) {
	using namespace ug::bridge;

	//	iterate through the parameter list and copy the value in the
	//  associated stack entry.
	for (int i = 0; i < paramsTemplate.size(); ++i) {
		int type = paramsTemplate.get_type(i);

		jobject value = env->GetObjectArrayElement(array, i);

		switch (type) {
			case PT_BOOL:
			{
				paramsOut.push_bool(jObject2Boolean(env, value));
			}
				break;
			case PT_INTEGER:
			{
				paramsOut.push_integer(jObject2Int(env, value));
			}
				break;
			case PT_NUMBER:
			{
				paramsOut.push_number(jObject2Double(env, value));
			}
				break;
			case PT_STRING:
			{
				paramsOut.push_string(
						jObject2String(env, value).c_str(), true);
			}
				break;

			case PT_POINTER:
			{
				paramsOut.push_pointer(jObject2Pointer(env, value),
						paramsTemplate.class_names(i));
				//DON'T USE paramsTemplate here!!!
				// Use the original type string of value
				// VRL now checks the type string and does not allow
				// incompatible connections. Thus, this should not a
				// problem anymore.
			}
				break;
			case PT_CONST_POINTER:
			{
				paramsOut.push_const_pointer(
						jObject2Pointer(env, value),
						paramsTemplate.class_names(i));
			}
				break;
		}

	} // end for
}

jobject param2JObject(
		JNIEnv *env, ug::bridge::ParameterStack& params, int index) {
	using namespace ug::bridge;
	//	iterate through the parameter list and copy the value in the
	//	associated stack entry.
	int type = params.get_type(index);

	switch (type) {
		case PT_BOOL:
		{
			return boolean2JObject(env, params.to_bool(index));
		}
			break;
		case PT_INTEGER:
		{
			return int2JObject(env, params.to_integer(index));
		}
			break;
		case PT_NUMBER:
		{
			return double2JObject(env, params.to_number(index));
		}
			break;
		case PT_STRING:
		{
			return string2JObject(env, params.to_string(index));
		}
			break;
		case PT_POINTER:
		{
			return pointer2JObject(env, params.to_pointer(index));
		}
			break;
		case PT_CONST_POINTER:
		{
			return pointer2JObject(
					env, (void*) params.to_const_pointer(index));
		}
			break;
	}

	return jobject();
}
} // end vrl::
}// end ug::
