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

jobject createEmptyString(JNIEnv *env) {
	jclass cls = env->FindClass("java/lang/String");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	return env->NewObject(cls, methodID);
}

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

jobjectArray stringArrayC2J(
		JNIEnv *env,
		const char* const* strings,
		const unsigned int array_length) {
	jclass stringClass = env->FindClass("java/lang/String");

	jobjectArray result =
			env->NewObjectArray(array_length, stringClass, 0);

	// convert array elements
	for (unsigned int i = 0; i < array_length; i++) {
		const char* s = strings[ i ];
		jstring javaString = env->NewStringUTF(s);
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
		return env->NewObjectArray(0, cls, createEmptyString(env));
	}
}

jobjectArray stringArrayC2J(
		JNIEnv *env, const std::vector<const char*>* strings) {

	if (strings != NULL && strings->size() > 0) {
		// it is safe to give a pointer to the first vector element as
		// std::vector implementation uses contiguous memory
		return stringArrayC2J(env, &(*strings)[0], strings->size());
	} else {

		// create an empty string array
		jclass cls = env->FindClass("java/lang/String");
		return env->NewObjectArray(0, cls, createEmptyString(env));
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

std::string name2ClassName(std::string className) {
	// Class names start with uppercase. For functions we need to ensure that.
	// Otherwise we get name conflicts with methods.
	std::string result = std::string(className);
	result[0] = toupper(result[0]);

	return result;
}

std::string name2InterfaceName(std::string className) {
	return name2ClassName(className) + "Interface";
}

std::string name2MethodName(std::string methodName) {
	// Class names start with uppercase. For functions we need to ensure that.
	// Otherwise we get name conflicts with methods.
	std::string result = std::string(methodName);
	result[0] = tolower(result[0]);

	return result;
}

std::vector<std::string> getInterfaceNames(ug::bridge::Registry* reg,
		ug::bridge::IExportedClass const& clazz) {
	const std::vector<const ug::bridge::IExportedClass*> baseClasses =
			getParentClasses(reg, &clazz);

	std::vector<std::string> result;

	for (unsigned int i = 0; i < baseClasses.size(); i++) {
		result.push_back(name2InterfaceName(baseClasses[i]->name()));
	}

	return result;
}




//std::string getReturnValueType(const ug::bridge::ExportedMethod &method) {
//	return " java.lang.Object ";
//}
//
//void createMethodSignature(std::stringstream& result,
//		const ug::bridge::ExportedMethod &method) {
//	std::string methodName = name2MethodName(method.name());
//
//	const ug::bridge::ParameterStack& paramStackIn = method.params_in();
//
//	size_t numParams = paramStackIn.size();
//
//	std::stringstream methodHeaderParams;
//
//	// generate method paramters for method header
//	// (including VRL param info) and for method invokation
//	for (unsigned int i = 0; i < numParams; i++) {
//
//		if (i > 0) {
//			methodHeaderParams << ",\n ";
//		}
//
//		methodHeaderParams << paramType2String(
//				paramStackIn.get_type(i),
//				method.parameter_name(i).c_str(),
//				paramStackIn.class_name(i),
//				paramStackIn.class_names(i),
//				method.parameter_info_vec(i)) << " p" << i;
//	}
//
//	result << "public " << getReturnValueType(method) << methodName
//			<< "(" << methodHeaderParams.str() << ")";
//}
//
//void createMethodBody(std::stringstream& result,
//		const ug::bridge::ExportedMethod &method) {
//	std::string methodName = name2MethodName(method.name());
//
//	const ug::bridge::ParameterStack& paramStackIn = method.params_in();
//
//	size_t numParams = paramStackIn.size();
//
//	std::stringstream paramArrayForInvocation;
//
//	// generate method paramters for method header
//	// (including VRL param info) and for method invokation
//	for (unsigned int i = 0; i < numParams; i++) {
//
//		if (i > 0) {
//			paramArrayForInvocation << ",\n ";
//		}
//
//		paramArrayForInvocation << " p" << i;
//	}
//
//	result << " {\n"
//			<< "edu.gcsc.vrl.ug4.UG4.getUG4().invokeMethod("
//			<< "getClassName(),"
//			<< " this, false, \""
//			<< methodName << "\", params);\n";
//}
//
//void createMethods(std::stringstream& result,
//		ug::bridge::IExportedClass const& clazz, bool isInterface) {
//
//	for (unsigned int i = 0; i < clazz->num_methods(); i++) {
//		const ug::bridge::ExportedMethod &method = clazz->get_method(i);
//
//		createMethodSignature(result, method);
//
//		result << "edu.gcsc.vrl.ug4.UG4.getUG4().invokeMethod("
//				<< "getClassName(),"
//				<< " getPointer().getAddress(), false, \""
//				<< method.name() << "\", params);\n";
//	}
//
//}
//
//std::string createClass(ug::bridge::Registry* reg,
//		ug::bridge::IExportedClass const& clazz) {
//	std::stringstream result;
//
//	std::string group = clazz.group();
//	std::string className = name2ClassName(clazz.name());
//
//	// create component info that specifies the menu group this
//	// class shall be added to
//	result << "@ComponentInfo(name=\"" << className
//			<< "\", category=\"" << group << "\", allowRemoval=false)\n"
//			<< "public class " << className
//			<< " extends edu.gcsc.vrl.ug4.UGObject"
//			<< " implements ";
//
//	std::vector<std::string> interfaces = getInterfaceNames(reg, clazz);
//
//	for (unsigned int i = 0; i < interfaces.size(); i++) {
//		if (i > 0) {
//			result << ", ";
//		}
//		result << interfaces[i];
//	}
//
//	result << " {\n private static final long serialVersionUID=1L;\n"
//			// create constructor
//			<< "public " << className << "() {\n"
//			<< "setClassName(\"" << className << "\");\n}\n";
//
//	// gather inheritance information (necessary for type-safety)
//	std::vector<const ug::bridge::IExportedClass*> baseClasses =
//			getParentClasses(reg, &clazz);
//
//	// adding methods
//	for (unsigned int i = 0; i < baseClasses.size(); i++) {
//		createMethods(result, baseClasses[i]);
//		//		generateConstMethods(result, baseClasses[i]);
//	}
//
//
//	return result.str();
//}

void generateMethodHeader(
		std::stringstream& result,
		ug::bridge::ExportedFunctionBase const& method,
		bool isFunction, bool isVisual, std::string prefix) {

	std::stringstream methodHeaderParams;
	std::stringstream paramArrayForInvocation;
	std::stringstream pointerListCode;

	std::string name = method.name();

	const ug::bridge::ParameterStack& paramStackIn = method.params_in();
	const ug::bridge::ParameterStack& paramStackOut = method.params_out();

	size_t numParams = paramStackIn.size();

	// generate method paramters for method header
	// (including VRL param info) and for method invokation
	for (unsigned int i = 0; i < numParams; i++) {

		if (paramStackIn.get_type(i) == ug::bridge::PT_CONST_POINTER
				|| paramStackIn.get_type(i) == ug::bridge::PT_POINTER) {
			pointerListCode << "addPointer(params[" << i << "]);\n";
		}

		if (i > 0) {
			methodHeaderParams << ",\n ";
			paramArrayForInvocation << ", ";
		}

		methodHeaderParams << paramType2String(
				paramStackIn.get_type(i),
				method.parameter_name(i).c_str(),
				paramStackIn.class_name(i),
				paramStackIn.class_names(i),
				method.parameter_info_vec(i)) << " p" << i;

		paramArrayForInvocation << " p" << i;
	}

	if (!isFunction && isVisual) {
		// we always need the visual id to get a reference to the
		// visualization that invokes this method,
		// that is why we add a visual id request to the param list
		if (numParams > 0) {
			methodHeaderParams << ", ";
			//			paramArrayForInvocation << ", ";
		}
		methodHeaderParams << " VisualIDRequest id ";
		//		paramArrayForInvocation << " id";
	}

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
	if (isVisual) {
		result << createMethodInfo(
				className.c_str(),
				classNames,
				readOnly, method.options()) << "\n";
	} else {
		result << "@MethodInfo(noGUI=true)\n";
	}

	// Methods start with lowercase. For functions we need to ensure that.
	// Otherwise we get name conflicts with constructors.
	std::string methodName = std::string(name);
	methodName[0] = tolower(methodName[0]);

	// putting it all together
	result << "public " << outType << " " << prefix << methodName << " ("
			<< methodHeaderParams.str() << ") {\n"
			<< "Object[] params = ["
			<< paramArrayForInvocation.str() << "] \n";

	// If we are not generating a function it is necessary to
	// call the updatePointer method with the current visual id.
	// Its purpose is to visually invoke
	// the setPointer() and getPointer() method before this method is
	// invoked. Otherwise we might operate on wrong instance
	if (!isFunction && isVisual) {
		result << "updatePointer(id);\n";
	} else {
		result << "updatePointer(null);\n";
	}

	result << pointerListCode.str() << "\n";

	if (returnsPointer(method)) {
		result << "edu.gcsc.vrl.ug4.Pointer result = ";
	} else if (paramStackOut.size() > 0) {
		result << "return ";
	}
}

bool returnsPointer(ug::bridge::ExportedFunctionBase const& func) {
	return func.params_out().size() > 0 &&
			(func.params_out().get_type(0) == ug::bridge::PT_CONST_POINTER
			|| func.params_out().get_type(0) == ug::bridge::PT_POINTER);
}

std::string exportedFunction2Groovy(
		ug::bridge::ExportedFunction const& func) {
	std::stringstream result;

	std::string group = func.group();

	// Class names start with uppercase. For functions we need to ensure that.
	// Otherwise we get name conflicts with methods.
	std::string className = std::string(func.name());
	className[0] = toupper(className[0]);

	// create component info that specifies the menu group this
	// function shall be added to
	result << "@ComponentInfo(name=\"" << func.name()
			<< "\", category=\"" << group << "\", allowRemoval=false)\n"
			<< "public class " << className
			<< " extends edu.gcsc.vrl.ug4.UGObject {\n"
			<< "private static final long serialVersionUID=1L;\n"
			<< "public " << className << "() {\n"
			<< "setClassName(\"" << className << "\");\n}\n";

	// visual function generation
	generateMethodHeader(
			result, func, true, true);
	result << "edu.gcsc.vrl.ug4.UG4.getUG4().invokeFunction( \""
			<< func.name() << "\", false, params);\n";

	if (returnsPointer(func)) {
		result << "result.setClassName(\""
				<< func.params_out().class_name(0) << "\");\n"
				<< "addPointer(result);\n"
				<< "return result;\n";
	}

	result << "}\n}\n";

	return result.str();
}

void generateMethods(std::stringstream& result,
		const ug::bridge::IExportedClass* clazz) {
	for (unsigned int i = 0; i < clazz->num_methods(); i++) {
		const ug::bridge::ExportedMethod &method = clazz->get_method(i);

		// non-visual method generation
		generateMethodHeader(result, method, false, false);

		result << "edu.gcsc.vrl.ug4.UG4.getUG4().invokeMethod("
				<< "getClassName(),"
				<< " getPointer().getAddress(), false, \""
				<< method.name() << "\", params);\n";

		if (returnsPointer(method)) {
			result << "result.setClassName(\""
					<< method.params_out().class_name(0) << "\");\n"
					<< "addPointer(result);\n"
					<< "return result;\n";
		}

		result << "}\n\n";

		// visual method generation
		generateMethodHeader(result, method, false, true);

		result << "edu.gcsc.vrl.ug4.UG4.getUG4().invokeMethod("
				<< "getClassName(),"
				<< " getPointer().getAddress(), false, \""
				<< method.name() << "\", params);\n";

		if (returnsPointer(method)) {
			result << "result.setClassName(\""
					<< method.params_out().class_name(0) << "\");\n"
					<< "addPointer(result);\n"
					<< "return result;\n";
		}

		result << "}\n\n";
	}
}

void generateConstMethods(std::stringstream& result,
		const ug::bridge::IExportedClass* clazz) {
	for (unsigned int i = 0; i < clazz->num_const_methods(); i++) {
		const ug::bridge::ExportedMethod &method = clazz->get_const_method(i);

		// non-visual method generation
		generateMethodHeader(result, method, false, false, "const_");

		result << "edu.gcsc.vrl.ug4.UG4.getUG4().invokeMethod("
				<< "getClassName(),"
				<< " getPointer().getAddress(), true, \""
				<< method.name() << "\", params);\n";

		if (returnsPointer(method)) {
			result << "result.setClassName(\""
					<< method.params_out().class_name(0) << "\");\n"
					<< "addPointer(result);\n"
					<< "return result;\n";
		}

		result << "}\n\n";

		// visual method generation
		generateMethodHeader(result, method, false, true, "const_");

		result << "edu.gcsc.vrl.ug4.UG4.getUG4().invokeMethod("
				<< "getClassName(),"
				<< " getPointer().getAddress(), true, \""
				<< method.name() << "\", params);\n";

		if (returnsPointer(method)) {
			result << "result.setClassName(\""
					<< method.params_out().class_name(0) << "\");\n"
					<< "addPointer(result);\n"
					<< "return result;\n";
		}

		result << "}\n\n";
	}
}

std::string exportedClass2Groovy(ug::bridge::Registry* reg,
		ug::bridge::IExportedClass const& clazz) {
	std::stringstream result;

	std::string group = clazz.group();

	// Class names start with uppercase. For functions we need to ensure that.
	// Otherwise we get name conflicts with methods.
	std::string className = /*"UG4_" + */ std::string(clazz.name());
	className[0] = toupper(className[0]);

	// create component info that specifies the menu group this
	// class shall be added to
	result << "@ComponentInfo(name=\"" << clazz.name()
			<< "\", category=\"" << group << "\", allowRemoval=false)\n"
			<< "public class " << className
			<< " extends edu.gcsc.vrl.ug4.UGObject ";

	// inheritance
	std::vector<std::string> interfaces = getInterfaceNames(reg, clazz);
	for (unsigned int i = 0; i < interfaces.size(); i++) {
		if (i > 0) {
			result << ", ";
		} else {
			result << " implements ";
		}
		result << interfaces[i];
	}

	result << "{\n"
			<< "private static final long serialVersionUID=1L;\n"

			<< "public " << className << "() {\n"
			<< "setClassName(\"" << className << "\");\n}\n";

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

	std::string paramInfoString = paramInfoStream.str();

	// We do not support invokeOnChange anymore!
	paramInfoString = replaceAll(paramInfoString, "invokeOnChange=true", "");
	paramInfoString = replaceAll(paramInfoString, ";invokeOnChange=true", "");
	paramInfoString = replaceAll(paramInfoString, "invokeOnChange=true;", "");

	return paramInfoString;
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
				//				return "edu.gcsc.vrl.ug4.Pointer";
				return name2ClassName(className);
			} else {
				std::string result =
						createParamInfo(paramName,
						className, classNames, false, paramInfo) +
						name2ClassName(className);
				//				std::string("edu.gcsc.vrl.ug4.Pointer");

				return result.c_str();
			}
		}
		case ug::bridge::PT_CONST_POINTER:
		{
			if (isOutput) {
				//				return "edu.gcsc.vrl.ug4.Pointer";
				return name2ClassName(className);
			} else {
				std::string result =
						createParamInfo(paramName,
						className, classNames, true, paramInfo) +
						name2ClassName(className);
				//						std::string("edu.gcsc.vrl.ug4.Pointer");

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
				//
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

int paramType2Int(const ug::bridge::ParameterStack& params, int index) {
	using namespace ug::bridge;
	//	iterate through the parameter list and return corresponding int
	int type = params.get_type(index);

	switch (type) {
		case PT_UNKNOWN:
		{
			return 0;
		}
			break;
		case PT_BOOL:
		{
			return 1;
		}
			break;
		case PT_INTEGER:
		{
			return 2;
		}
			break;
		case PT_NUMBER:
		{
			return 3;
		}
			break;
		case PT_STRING:
		{
			return 4;
		}
			break;
		case PT_POINTER:
		{
			return 5;
		}
			break;
		case PT_CONST_POINTER:
		{
			return 6;
		}
			break;
		default:
			return 0;
			break;
	}

	return 0;
}

jobjectArray params2NativeParams(JNIEnv *env,
		const ug::bridge::ExportedFunctionBase& func) {

	jclass cls = env->FindClass("edu/gcsc/vrl/ug4/NativeParamInfo");

	jobjectArray result =
			env->NewObjectArray(func.num_parameter(), cls, 0);

	const ug::bridge::ParameterStack& params = func.params_in();

	for (unsigned int i = 0; i < func.num_parameter(); i++) {

		//		std::cout << "P->P: " << i << std::endl;

		// create instance
		jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
		jobject obj = env->NewObject(cls, methodID);

		//assign values

		jmethodID setType = env->GetMethodID(cls, "setType", "(I)V");
		jmethodID setID = env->GetMethodID(cls, "setId", "(I)V");
		jmethodID setClassName = env->GetMethodID(
				cls, "setClassName", "(Ljava/lang/String;)V");
		jmethodID setClassNames = env->GetMethodID(
				cls, "setClassNames", "([Ljava/lang/String;)V");
		//	jmethodID setHelp = env->GetMethodID(cls,
		//		"setHelp", "(Ljava/lang/String;)V");
		//	jmethodID setToolTip = env->GetMethodID(cls,
		//		"setTooltip", "(Ljava/lang/String;)V");
		jmethodID setParamInfo = env->GetMethodID(cls,
				"setParamInfo", "([Ljava/lang/String;)V");

		using namespace ug::bridge;

		// TODO Unfortunately we don't know a better way to convert an enumeration
		// from C++ to Java. Currently we just use integers :(
		int type = paramType2Int(params, i);

		//		if (type == 0 || type > 6) {
		//			std::cout << "*******METHOD: " << func.name() << std::endl;
		//			std::cout << "*******PARAM: " << type << " " << params.class_name(i) << "\n";
		//		}

		env->CallVoidMethod(obj, setType, type);
		env->CallVoidMethod(obj, setID, i);

		env->CallVoidMethod(obj, setClassName, stringC2J(env, params.class_name(i)));
		env->CallVoidMethod(obj, setClassNames, stringArrayC2J(env, params.class_names(i)));
		env->CallVoidMethod(obj, setParamInfo,
				stringArrayC2J(env, func.parameter_info_vec(i)));

		// set array element
		env->SetObjectArrayElement(result, i, obj);
	}

	return result;
}

jobject retVal2NativeParam(JNIEnv *env,
		const ug::bridge::ExportedFunctionBase& func) {

	//	std::cout << "***R->HERE0***" << std::endl;

	jclass cls = env->FindClass("edu/gcsc/vrl/ug4/NativeParamInfo");

	unsigned int i = 0; // C/C++/Java only allow one return value

	const ug::bridge::ParameterStack& params = func.params_out();

	// create instance
	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	jobject obj = env->NewObject(cls, methodID);

	//assign values

	//	std::cout << "***R->HERE1***" << std::endl;

	jmethodID setType = env->GetMethodID(cls, "setType", "(I)V");
	jmethodID setID = env->GetMethodID(cls, "setId", "(I)V");
	jmethodID setClassName = env->GetMethodID(
			cls, "setClassName", "(Ljava/lang/String;)V");
	jmethodID setClassNames = env->GetMethodID(
			cls, "setClassNames", "([Ljava/lang/String;)V");
	//	jmethodID setHelp = env->GetMethodID(cls,
	//	"setHelp", "(Ljava/lang/String;)V");
	//	jmethodID setToolTip = env->GetMethodID(cls,
	//	"setTooltip", "(Ljava/lang/String;)V");
	jmethodID setParamInfo = env->GetMethodID(cls,
			"setParamInfo", "([Ljava/lang/String;)V");


	bool returnsVoid = func.params_out().size() == 0;


	// TODO Unfortunately we don't know a better way to convert an enumeration
	// from C++ to Java. Currently we just use integers :(
	int type = -1;

	if (returnsVoid) {
		//		std::cout << "***R->VOID" << std::endl;
		type = -1; // void
	} else {
		type = paramType2Int(params, i);

		env->CallVoidMethod(obj, setClassName, stringC2J(env, params.class_name(i)));
		env->CallVoidMethod(obj, setClassNames, stringArrayC2J(env, params.class_names(i)));
	}

	env->CallVoidMethod(obj, setType, type);
	env->CallVoidMethod(obj, setID, i);

	if (!returnsVoid) {
		env->CallVoidMethod(obj, setParamInfo,
				stringArrayC2J(env, func.return_info_vec()));
	}

	//	std::cout << "***R->HERE2***" << std::endl;

	return obj;
}

jobjectArray methods2NativeMethods(JNIEnv *env,
		const ug::bridge::IExportedClass& eCls, bool constMethods) {
	jclass cls = env->FindClass("edu/gcsc/vrl/ug4/NativeMethodInfo");

	//	std::cout << "***M->HERE0***" << std::endl;

	unsigned int numMethods = 0;

	if (constMethods) {
		numMethods = eCls.num_const_methods();
	} else {
		numMethods = eCls.num_methods();
	}

	jobjectArray result =
			env->NewObjectArray(numMethods, cls, 0);

	for (unsigned int i = 0; i < numMethods; i++) {
		const ug::bridge::ExportedMethod* method;

		//		std::cout << "M->Method: " << i << std::endl;

		if (constMethods) {
			method = &eCls.get_const_method(i);
		} else {
			method = &eCls.get_method(i);
		}

		// create instance

		jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
		jobject obj = env->NewObject(cls, methodID);

		//assign values

		jmethodID setName = env->GetMethodID(cls,
				"setName", "(Ljava/lang/String;)V");
		jmethodID setHelp = env->GetMethodID(cls,
				"setHelp", "(Ljava/lang/String;)V");
		jmethodID setToolTip = env->GetMethodID(cls,
				"setToolTip", "(Ljava/lang/String;)V");
		jmethodID setOptions = env->GetMethodID(cls,
				"setOptions", "(Ljava/lang/String;)V");
		jmethodID setRetValue = env->GetMethodID(cls,
				"setReturnValue", "(Ledu/gcsc/vrl/ug4/NativeParamInfo;)V");
		jmethodID setParameters = env->GetMethodID(cls,
				"setParameters", "([Ledu/gcsc/vrl/ug4/NativeParamInfo;)V");

		using namespace ug::bridge;
		std::string name = method->name(); // TODO pre-rpocessing necessary
		env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
		env->CallVoidMethod(obj, setHelp, stringC2J(env, method->help().c_str()));
		env->CallVoidMethod(obj, setToolTip, stringC2J(env, name.c_str()));
		env->CallVoidMethod(obj, setOptions, stringC2J(env, method->options().c_str()));
		env->CallVoidMethod(obj, setRetValue,
				retVal2NativeParam(env, *method));

		//		std::cout << "***M->2***" << std::endl;

		env->CallVoidMethod(obj, setParameters,
				params2NativeParams(env, *method));

		//		std::cout << "***M->3***" << std::endl;

		// set array element
		env->SetObjectArrayElement(result, i, obj);

		//		std::cout << "***M->4***" << std::endl;
	}

	return result;
}

//jobjectArray functions2NativeFunctions(JNIEnv *env, ug::bridge::Registry* reg) {
//	jclass cls = env->FindClass("edu/gcsc/vrl/ug4/NativeFunctionInfo");
//
//	//	std::cout << "***M->HERE0***" << std::endl;
//
//	unsigned int numFunctions = reg->num_functions();
//
//	jobjectArray result =
//			env->NewObjectArray(numFunctions, cls, 0);
//
//	for (unsigned int i = 0; i < numFunctions; i++) {
//		const ug::bridge::ExportedFunction& func = reg->get_function(i);
//
//		// create instance
//
//		jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
//		jobject obj = env->NewObject(cls, methodID);
//
//		//assign values
//
//		jmethodID setName = env->GetMethodID(cls,
//				"setName", "(Ljava/lang/String;)V");
//		jmethodID setConst = env->GetMethodID(cls,
//				"setConst", "(Z)V");
//		jmethodID setCategory = env->GetMethodID(cls,
//				"setCategory", "(Ljava/lang/String;)V");
//		jmethodID setHelp = env->GetMethodID(cls,
//				"setHelp", "(Ljava/lang/String;)V");
//		jmethodID setToolTip = env->GetMethodID(cls,
//				"setToolTip", "(Ljava/lang/String;)V");
//		jmethodID setOptions = env->GetMethodID(cls,
//				"setOptions", "(Ljava/lang/String;)V");
//		jmethodID setRetValue = env->GetMethodID(cls,
//				"setReturnValue", "(Ledu/gcsc/vrl/ug4/NativeParamInfo;)V");
//		jmethodID setParameters = env->GetMethodID(cls,
//				"setParameters", "([Ledu/gcsc/vrl/ug4/NativeParamInfo;)V");
//
//		using namespace ug::bridge;
//		std::string name = func.name(); // TODO pre-rpocessing necessary
//		env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
//		env->CallVoidMethod(obj, setConst, boolC2J(false));
//		env->CallVoidMethod(obj, setCategory, stringC2J(env, func.group().c_str()));
//		env->CallVoidMethod(obj, setHelp, stringC2J(env, func.help().c_str()));
//		env->CallVoidMethod(obj, setToolTip, stringC2J(env, name.c_str()));
//		env->CallVoidMethod(obj, setOptions, stringC2J(env, func.options().c_str()));
//		env->CallVoidMethod(obj, setRetValue,
//				retVal2NativeParam(env, func));
//
//		env->CallVoidMethod(obj, setParameters,
//				params2NativeParams(env, func));
//
//		// set array element
//		env->SetObjectArrayElement(result, i, obj);
//	}
//
//	return result;
//}

jobject function2NativeFunction(JNIEnv *env, const ug::bridge::ExportedFunction& func) {
	jclass cls = env->FindClass("edu/gcsc/vrl/ug4/NativeFunctionInfo");

	// create instance

	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	jobject obj = env->NewObject(cls, methodID);

	//assign values

	jmethodID setName = env->GetMethodID(cls,
			"setName", "(Ljava/lang/String;)V");
	jmethodID setConst = env->GetMethodID(cls,
			"setConst", "(Z)V");
	jmethodID setCategory = env->GetMethodID(cls,
			"setCategory", "(Ljava/lang/String;)V");
	jmethodID setHelp = env->GetMethodID(cls,
			"setHelp", "(Ljava/lang/String;)V");
	jmethodID setToolTip = env->GetMethodID(cls,
			"setToolTip", "(Ljava/lang/String;)V");
	jmethodID setOptions = env->GetMethodID(cls,
			"setOptions", "(Ljava/lang/String;)V");
	jmethodID setRetValue = env->GetMethodID(cls,
			"setReturnValue", "(Ledu/gcsc/vrl/ug4/NativeParamInfo;)V");
	jmethodID setParameters = env->GetMethodID(cls,
			"setParameters", "([Ledu/gcsc/vrl/ug4/NativeParamInfo;)V");

	using namespace ug::bridge;
	std::string name = func.name(); // TODO pre-rpocessing necessary
	env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
	env->CallVoidMethod(obj, setConst, boolC2J(false));
	env->CallVoidMethod(obj, setCategory, stringC2J(env, func.group().c_str()));
	env->CallVoidMethod(obj, setHelp, stringC2J(env, func.help().c_str()));
	env->CallVoidMethod(obj, setToolTip, stringC2J(env, name.c_str()));
	env->CallVoidMethod(obj, setOptions, stringC2J(env, func.options().c_str()));
	env->CallVoidMethod(obj, setRetValue,
			retVal2NativeParam(env, func));

	env->CallVoidMethod(obj, setParameters,
			params2NativeParams(env, func));

	return obj;
}

jobjectArray functions2NativeGroups(JNIEnv *env, ug::bridge::Registry* reg) {
	jclass groupArrayCls =
			env->FindClass("edu/gcsc/vrl/ug4/NativeFunctionGroupInfo");

	//	std::cout << "***M->HERE0***" << std::endl;

	unsigned int numFunctions = reg->num_functions();

	// create array of functions
	jobjectArray result =
			env->NewObjectArray(numFunctions, groupArrayCls, 0);

	for (unsigned int i = 0; i < numFunctions; i++) {

		const ug::bridge::ExportedFunctionGroup& group =
				reg->get_function_group(i);

		unsigned int numOverloads = group.num_overloads();

		jclass functionsCls =
				env->FindClass("edu/gcsc/vrl/ug4/NativeFunctionInfo");

		jobjectArray functions =
				env->NewObjectArray(numOverloads, functionsCls, 0);

		for (unsigned int j = 0; j < numOverloads; j++) {
			const ug::bridge::ExportedFunction& func = *group.get_overload(j);
//			std::cout << "Func: " << func.name() << ", overload: " << j << std::endl;

			env->SetObjectArrayElement(functions, j,
					function2NativeFunction(env, func));

		} // end for j

		// create function group instance
		jmethodID methodID = env->GetMethodID(groupArrayCls, "<init>", "()V");
		jobject obj = env->NewObject(groupArrayCls, methodID);

		jmethodID setOverloads = env->GetMethodID(groupArrayCls,
				"setOverloads", "([Ledu/gcsc/vrl/ug4/NativeFunctionInfo;)V");

		env->CallVoidMethod(obj, setOverloads, functions);

		env->SetObjectArrayElement(result, i, obj);
	} // end for i

	return result;
}

jobjectArray classes2NativeClasses(JNIEnv *env, const ug::bridge::Registry* reg) {

	jclass cls = env->FindClass("edu/gcsc/vrl/ug4/NativeClassInfo");

	jobjectArray result =
			env->NewObjectArray(reg->num_classes(), cls, 0);

	for (unsigned int i = 0; i < reg->num_classes(); i++) {

		//		std::cout << "----- CLASS: " << i << " -----" << std::endl;

		const ug::bridge::IExportedClass& eCls = reg->get_class(i);

		// create instance

		jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
		jobject obj = env->NewObject(cls, methodID);

		//assign values

		jmethodID setName = env->GetMethodID(cls,
				"setName", "(Ljava/lang/String;)V");
		jmethodID setCategory = env->GetMethodID(cls,
				"setCategory", "(Ljava/lang/String;)V");
		jmethodID setClassNames = env->GetMethodID(
				cls, "setClassNames", "([Ljava/lang/String;)V");
		jmethodID setInstantiable = env->GetMethodID(cls,
				"setInstantiable", "(Z)V");
		jmethodID setMethods = env->GetMethodID(cls,
				"setMethods", "([Ledu/gcsc/vrl/ug4/NativeMethodInfo;)V");
		jmethodID setConstMethods = env->GetMethodID(cls,
				"setConstMethods", "([Ledu/gcsc/vrl/ug4/NativeMethodInfo;)V");

		//		std::cout << "***C->HERE-3***" << std::endl;

		//		using namespace ug::bridge;
		std::string name = eCls.name(); // TODO pre-rpocessing necessary
		env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));

		//		std::cout << "***C->HERE-2***" << std::endl;

		env->CallVoidMethod(obj, setCategory, stringC2J(env, eCls.group().c_str()));

		env->CallVoidMethod(obj, setClassNames, stringArrayC2J(env, eCls.class_names()));

		//		std::cout << "***C->HERE-1***" << std::endl;

		env->CallVoidMethod(obj, setInstantiable, boolC2J(eCls.is_instantiable()));
		//

		//		std::cout << "***C->HERE0***" << std::endl;

		env->CallVoidMethod(obj, setMethods, methods2NativeMethods(env, eCls, false));

		//		std::cout << "***C->HERE1***" << std::endl;

		env->CallVoidMethod(obj, setConstMethods, methods2NativeMethods(env, eCls, true));

		//		std::cout << "***C->HERE2***" << std::endl;

		// set array element
		env->SetObjectArrayElement(result, i, obj);

		//		std::cout << "***C->HERE3***" << std::endl;
	}

	//	std::cout << "***C->END" << std::endl;

	//	env->ExceptionCheck();
	//	env->ExceptionDescribe();

	return result;
}

jobject registry2NativeAPI(JNIEnv *env, ug::bridge::Registry* reg) {
	jclass cls = env->FindClass("edu/gcsc/vrl/ug4/NativeAPIInfo");

	// create instance

	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	jobject obj = env->NewObject(cls, methodID);

	//assign values

	jmethodID setClasses = env->GetMethodID(cls, "setClasses",
			"([Ledu/gcsc/vrl/ug4/NativeClassInfo;)V");
	jmethodID setFunctions = env->GetMethodID(cls, "setFunctions",
			"([Ledu/gcsc/vrl/ug4/NativeFunctionGroupInfo;)V");

	env->CallVoidMethod(obj, setClasses, classes2NativeClasses(env, reg));
	env->CallVoidMethod(obj, setFunctions, functions2NativeGroups(env, reg));

	return obj;
}


} // end vrl::
}// end ug::
