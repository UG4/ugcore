#include "type_converter.h"
#include "registry/registry.h"
#include "messaging.h"
#include "invocation.h"
#include <sstream>

namespace ug {
namespace vrl {

jobject createEmptyString(JNIEnv *env) {
	jclass cls = env->FindClass("java/lang/String");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	return env->NewObject(cls, methodID);
}

jstring stringC2J(JNIEnv *env, const char* s) {
	return env->NewStringUTF(s);
}

std::string stringJ2C(JNIEnv *env, jstring const& s) {
	const char* tmpStr = env->GetStringUTFChars(s, JNI_FALSE);
	std::string result = (std::string) tmpStr;
	env->ReleaseStringUTFChars(s, tmpStr);
	return result;
}

jobjectArray stringArrayC2J(JNIEnv *env, const std::string* strings,
		const size_t array_length) {
	jclass stringClass = env->FindClass("java/lang/String");

	jobjectArray result = env->NewObjectArray(array_length, stringClass, 0);

	// convert array elements
	for (size_t i = 0; i < array_length; i++) {
		std::string s = strings[i];
		jstring javaString = env->NewStringUTF(s.c_str());
		env->SetObjectArrayElement(result, i, javaString);
	}

	return result;
}

jobjectArray stringArrayC2J(JNIEnv *env, const char* const * strings,
		const unsigned int array_length) {
	jclass stringClass = env->FindClass("java/lang/String");

	jobjectArray result = env->NewObjectArray(array_length, stringClass, 0);

	// convert array elements
	for (size_t i = 0; i < array_length; i++) {
		const char* s = strings[i];
		jstring javaString = env->NewStringUTF(s);
		env->SetObjectArrayElement(result, i, javaString);
	}

	return result;
}

jobjectArray createEmptyJavaArray(JNIEnv *env, std::string className) {
	jclass elementClass = env->FindClass(className.c_str());

	return createEmptyJavaArray(env, elementClass);
}

jobjectArray createEmptyJavaArray(JNIEnv *env, jclass elementClass) {

	return env->NewObjectArray(0, elementClass, 0);
}

jobjectArray stringArrayC2J(JNIEnv *env,
		std::vector<std::string> const& strings) {

	if (!strings.empty()) {
		// it is safe to give a pointer to the first vector element as
		// std::vector implementation uses contiguous memory
		return stringArrayC2J(env, &strings[0], strings.size());
	} else {
		// create an empty string array
		jclass cls = env->FindClass("java/lang/String");
		return env->NewObjectArray(0, cls, createEmptyString(env));
	}
}

jobjectArray stringArrayC2J(JNIEnv *env,
		const std::vector<const char*>* strings) {

	if (strings != NULL && !strings->empty()) {
		// it is safe to give a pointer to the first vector element as
		// std::vector implementation uses contiguous memory
		return stringArrayC2J(env, &(*strings)[0], strings->size());
	} else {

		// create an empty string array
		jclass cls = env->FindClass("java/lang/String");
		return env->NewObjectArray(0, cls, createEmptyString(env));
	}
}

jobjectArray stringArrayC2J(JNIEnv *env,
		const std::vector<const char*> strings) {

	if (!strings.empty()) {
		// it is safe to give a pointer to the first vector element as
		// std::vector implementation uses contiguous memory
		return stringArrayC2J(env, &strings[0], strings.size());
	} else {

		// create an empty string array
		jclass cls = env->FindClass("java/lang/String");
		return env->NewObjectArray(0, cls, createEmptyString(env));
	}
}

std::vector<std::string> stringArrayJ2C(JNIEnv *env,
		jobjectArray const& array) {

	std::vector<std::string> result;

	size_t length = env->GetArrayLength(array);

	// convert each element of the java object array to a std string
	// and add it to the result vector
	for (size_t i = 0; i < length; i++) {
		result.push_back(
				stringJ2C(env, (jstring) env->GetObjectArrayElement(array, i)));
	}

	return result;
}

/*//	added by Christian Poliwoda
 //	christian.poliwoda@gcsc.uni-frankfurt.de
 //	y13 m05 d07
 std::vector<bool> booleanArrayJ2C(JNIEnv *env, jobjectArray const& array) {

 std::vector<bool> result;

 size_t length = env->GetArrayLength(array);

 // convert each element of the java object array to a boolean
 // and add it to the result vector
 for (size_t i = 0; i < length; i++) {
 result.push_back(
 boolJ2C((jboolean) env->GetObjectArrayElement(array, i)));
 }

 return result;
 }*/

const std::vector<const ug::bridge::IExportedClass*> getParentClasses(
		ug::bridge::Registry* reg, const ug::bridge::IExportedClass* clazz) {
	std::vector<const ug::bridge::IExportedClass*> result;
	// search registered classes by name as specified in the
	// class_names vector and add them to the result vector
	for (size_t i = 0; i < clazz->class_names()->size(); i++) {
		const ug::bridge::IExportedClass* baseCls =
				ug::vrl::invocation::getExportedClassPtrByName(reg,
						(*clazz->class_names())[i]);
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
	jmethodID methodID = env->GetMethodID(argClass, "booleanValue", "()Z");
	return env->CallBooleanMethod(obj, methodID);
}

std::vector<const char*> getBaseClassNames(
		const ug::bridge::ClassNameNode* node) {
	std::vector<const char*> classNames;

	ug::bridge::ExtractClassNameVec(classNames, *node, false);

	return classNames;
}

void* jObject2Pointer(JNIEnv *env, jobject obj) {
	jclass argClass = env->GetObjectClass(obj);
	jmethodID methodID = env->GetMethodID(argClass, "getAddress", "()J");
	long result = env->CallLongMethod(obj, methodID);

	if (result == 0) {
		jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
		env->ThrowNew(Exception, "Pointer is NULL!");
	}

	return (void*) result;
}

std::string jPointerGetName(JNIEnv *env, jobject obj) {
	jclass argClass = env->GetObjectClass(obj);
	jmethodID methodID = env->GetMethodID(argClass, "getClassName",
			"()Ljava/lang/String;");

	std::string name = stringJ2C(env,
			(jstring) env->CallObjectMethod(obj, methodID));

	//checkException(env,"ugbase/bindings/vrl/type_converter.cpp : jPointerGetName()",true);

	return name;
}

SmartPtr<void> jObject2SmartPointer(JNIEnv *env, jobject obj) {

	// christian poliwoda ug-log for debug
	UG_LOG(
			"trunk/ugbase/bindings/vrl/type_converter.cpp : jObject2SmartPointer(JNIEnv *env, jobject obj)"<<std::endl);

	jclass argClass = env->GetObjectClass(obj);

	/*UG_LOG(" jObject2SmartPointer : argClass = "<< argClass << std::endl);
	 */
	jmethodID getSmartPointer = env->GetMethodID(argClass, "getSmartPointer",
			"()[B");
	jbyteArray mem = (jbyteArray) env->CallObjectMethod(obj, getSmartPointer);
	jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

	// temporarily it is save to use the memory provided by the Java
	// byte array
	SmartPtr<void>* smartPtr = reinterpret_cast<SmartPtr<void>*>((void*) memPtr);

	// but we have to use a new instance for use outside of this function
	SmartPtr<void> result(*smartPtr);

	/*
	 UG_LOG(
	 " jObject2SmartPointer : smartPtr.get() = "<< smartPtr->get() << std::endl);
	 UG_LOG(
	 " jObject2SmartPointer : result(*smartPtr).get() = "<< result.get() << std::endl);
	 */

	env->ReleaseByteArrayElements(mem, memPtr, 0);

	return result;
}

ConstSmartPtr<void> jObject2ConstSmartPointer(JNIEnv *env, jobject obj) {

	ConstSmartPtr<void> result;

	// We allow convertion from SmartPtr to ConstSmartPtr. As these classes
	// are not in an inheritance relation we must take care to use the correct
	// methods for conversion
	if (isJSmartPointerConst(env, obj)) {

		jclass argClass = env->GetObjectClass(obj);
		jmethodID getSmartPointer = env->GetMethodID(argClass,
				"getSmartPointer", "()[B");
		jbyteArray mem = (jbyteArray) env->CallObjectMethod(obj,
				getSmartPointer);
		jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

		ConstSmartPtr<void>* smartPtr =
				reinterpret_cast<ConstSmartPtr<void>*>((void*) memPtr);
		result = ConstSmartPtr<void>(*smartPtr);

		env->ReleaseByteArrayElements(mem, memPtr, 0);

	} else {
		// conversion from SmartPtr is done by ConstSmartPtr's assignment
		// operator
		result = jObject2SmartPointer(env, obj);
	}

	return result;
}

void invalidateJSmartPointer(JNIEnv *env, jobject obj) {

	jclass argClass = env->GetObjectClass(obj);

	jmethodID getSmartPointer = env->GetMethodID(argClass, "getSmartPointer",
			"()[B");
	jbyteArray mem = (jbyteArray) env->CallObjectMethod(obj, getSmartPointer);
	jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

	// It is important to call the invalidate method directly on the
	// smart-pointer instance that has been retrieved from the java object.
	// It is not valid to use copy-constructor etc. and to work on a stack
	// instance because invalidation of the stack instance won't allow the
	// destructor to decrease the reference count. This is equivalent to
	// not invoking the invalidate() method.
	(reinterpret_cast<SmartPtr<void>*>((void*) memPtr))->invalidate();

	env->ReleaseByteArrayElements(mem, memPtr, 0);
}

void invalidateJConstSmartPointer(JNIEnv *env, jobject obj) {

	jclass argClass = env->GetObjectClass(obj);

	jmethodID getSmartPointer = env->GetMethodID(argClass, "getSmartPointer",
			"()[B");
	jbyteArray mem = (jbyteArray) env->CallObjectMethod(obj, getSmartPointer);
	jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

	// It is important to call the invalidate method directly on the
	// smart-pointer instance that has been retrieved from the java object.
	// It is not valid to use copy-constructor etc. and to work on a stack
	// instance because invalidation of the stack instance won't allow the
	// destructor to decrease the reference count. This is equivalent to
	// not invoking the invalidate() method.
	(reinterpret_cast<ConstSmartPtr<void>*>((void*) memPtr))->invalidate();

	env->ReleaseByteArrayElements(mem, memPtr, 0);
}

jobject pointer2JObject(JNIEnv *env, void* value) {

	if (value == NULL) {
		return JNULL; // exception occured
	}

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/Pointer");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "(JZ)V");
	return env->NewObject(cls, methodID, (jlong) value, boolC2J(false));
}

jobject constPointer2JObject(JNIEnv *env, void* value) {

	if (value == NULL) {
		return JNULL; // exception occured
	}

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/Pointer");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "(JZ)V");
	return env->NewObject(cls, methodID, (jlong) value, boolC2J(true));
}

jobject smartPointer2JObject(JNIEnv *env, SmartPtr<void> value) {

	jsize size = sizeof(SmartPtr<void> );
	jbyteArray mem = env->NewByteArray(size);
	jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

	if (memPtr == NULL) {
		return JNULL; // exception occured
	}

	new (memPtr) SmartPtr<void>(value);

	// release the memory to write it back to the jvm
	env->ReleaseByteArrayElements(mem, memPtr, 0);

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/SmartPointer");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "(J[BZ)V");
	jobject result = env->NewObject(cls, methodID, (jlong) value.get(), mem,
			boolC2J(false));

	return result;
}

jobject constSmartPointer2JObject(JNIEnv *env, ConstSmartPtr<void> value) {

	jsize size = sizeof(ConstSmartPtr<void> );
	jbyteArray mem = env->NewByteArray(size);
	jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

	if (memPtr == NULL) {
		return JNULL; // exception occured
	}

	new (memPtr) ConstSmartPtr<void>(value);

	// release the memory to write it back to the jvm
	env->ReleaseByteArrayElements(mem, memPtr, 0);

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/SmartPointer");
	jmethodID methodID = env->GetMethodID(cls, "<init>", "(J[BZ)V");
	jobject result = env->NewObject(cls, methodID, (jlong) value.get(), mem,
			boolC2J(true));

	return result;
}

bool isJSmartPointerConst(JNIEnv *env, jobject ptr) {
	jclass cls = env->FindClass("edu/gcsc/vrl/ug/SmartPointer");
	jmethodID methodID = env->GetMethodID(cls, "isConst", "()Z");
	jboolean result = env->CallBooleanMethod(ptr, methodID);
	return boolJ2C(result);
}

//bool isJSmartPointer(JNIEnv *env, jobject ptr) {
//	jclass cls = env->FindClass("edu/gcsc/vrl/ug/SmartPointer");
//	jclass objCls = getClass(env,ptr);
//
//	jmethodID equals =
//			env->GetMethodID(objCls, "equals", "(Ljava.lang.Object;)Z");
//
//	bool result = env->CallBooleanMethod(objCls, equals, cls);
//
//	return boolJ2C(result);
//}

jobject string2JObject(JNIEnv *env, const char* value) {
	return env->NewStringUTF(value);
}

std::string jObject2String(JNIEnv *env, jobject obj) {
	return stringJ2C(env, (jstring) obj);
}

void throwUgErrorAsJavaException(JNIEnv *env, ug::UGError error) {
	jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");

	std::stringstream ss;

	ss << "<!--NumMsg:" << error.num_msg() << ":-->";
	ss << "<div><pre>\n";
	for (size_t i = 0; i < error.num_msg(); i++) {

		ss << "<!--StartMsg:" << i << ":-->";
		std::string msg = error.get_msg(i);
		ss << "\n" << msg << "\n";
		ss << "<!--EndMsg:" << i << ":-->";

		if (msg.size() != 0 && msg[msg.size() - 1] != '\n') {
			ss << "\n";
		}

		ss << ">> File:\t";
		ss << "<!--StartFile:" << i << ":-->";
		ss << error.get_file(i);
		ss << "<!--EndFile:" << i << ":-->";
		ss << "\n";

		ss << ">> Line:\t";
		ss << "<!--StartLine:" << i << ":-->";
		ss << error.get_line(i);
		ss << "<!--EndLine:" << i << ":-->";
		ss << "\n\n";

		if (i < error.num_msg() - 1) {
			ss << "> Caused by:" << std::endl;
		}
	}

	ss << "</pre></div>";

	env->ThrowNew(Exception, ss.str().c_str());
}

jobject getClass(JNIEnv *env, jobject obj) {
	jclass classMethodAccess = env->FindClass("java/lang/Class");

	jmethodID classNameMethodID = env->GetMethodID(classMethodAccess,
			"getClass", "()Ljava/lang/Class;");

	return env->CallObjectMethod(obj, classNameMethodID);
}

std::string getClassName(JNIEnv *env, jobject obj) {
	jclass classMethodAccess = env->FindClass("java/lang/Class");

	jmethodID classNameMethodID = env->GetMethodID(classMethodAccess,
			"getClass", "()Ljava/lang/Class;");

	jobject clazz = env->CallObjectMethod(obj, classNameMethodID);

	classNameMethodID = env->GetMethodID(classMethodAccess, "getName",
			"()Ljava/lang/String;");
	jobject resultObj = env->CallObjectMethod(clazz, classNameMethodID);

	return jObject2String(env, resultObj);
}

std::string getParamClassName(JNIEnv *env, jobject obj) {

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/Pointer");

	if (env->ExceptionCheck()) {
		env->ExceptionDescribe();
	}

	jmethodID instanceMethodID = env->GetStaticMethodID(cls, "isInstance",
			"(Ljava/lang/Object;)Z");

	jboolean instanceof = env->CallStaticBooleanMethod(cls, instanceMethodID,
			obj);

	if (boolJ2C(instanceof)) {
		jmethodID classNameMethodID = env->GetMethodID(cls, "getClassName",
				"()Ljava/lang/String;");

		if (env->ExceptionCheck()) {
			env->ExceptionDescribe();
		}
		jobject resultObj = env->CallObjectMethod(obj, classNameMethodID);

		if (resultObj == NULL) {
			return "";
		}

		return jObject2String(env, resultObj);
	}

	return "";
}

// used for debugging only
void printParamType(const uint type, size_t index) {
	using namespace ug::bridge;
	//	iterate through the parameter list and return corresponding int

	switch (type) {
	case ug::Variant::VT_INVALID: {
		UG_LOG("Param " << index << " = VT_INVALID" << std::endl)
	}
		break;
	case ug::Variant::VT_BOOL: {
		UG_LOG("Param " << index << " = VT_BOOL" << std::endl)
	}
		break;
	case ug::Variant::VT_INT: {
		UG_LOG("Param " << index << " = VT_INT" << std::endl)
	}
		break;
	case ug::Variant::VT_SIZE_T: {
		UG_LOG("Param " << index << " = VT_SIZE_T" << std::endl)
	}
		break;
	case ug::Variant::VT_FLOAT:
	case ug::Variant::VT_DOUBLE: {
		UG_LOG("Param " << index << " = VT_NUMBER" << std::endl)
	}
		break;
	case ug::Variant::VT_CSTRING: {
		UG_LOG("Param " << index << " = VT_CSTRING" << std::endl)
	}
		break;
	case ug::Variant::VT_STDSTRING: {
		UG_LOG("Param " << index << " = VT_STDSTRING" << std::endl)
	}
		break;
	case ug::Variant::VT_POINTER: {
		UG_LOG("Param " << index << " = VT_POINTER" << std::endl)
	}
		break;
	case ug::Variant::VT_CONST_POINTER: {
		UG_LOG("Param " << index << " = VT_CONST_POINTER" << std::endl)
	}
		break;
	case ug::Variant::VT_SMART_POINTER: {
		UG_LOG("Param " << index << " = VT_SMART_POINTER" << std::endl)
	}
		break;
	case ug::Variant::VT_CONST_SMART_POINTER: {
		UG_LOG("Param " << index << " = VT_CONST_SMART_POINTER" << std::endl)
	}
		break;
	default:
		break;
	}
}

// used for debugging only
std::string getParamTypeAsString(const uint type) {
	using namespace ug::bridge;
	//	iterate through the parameter list and return corresponding int

	std::string result = "undefined";

	switch (type) {
	case ug::Variant::VT_INVALID: {
		result = "VT_INVALID";
	}
		break;
	case ug::Variant::VT_BOOL: {
		result = "VT_BOOL";
	}
		break;
	case ug::Variant::VT_INT: {
		result = "VT_INT";
	}
		break;
	case ug::Variant::VT_FLOAT:
	case ug::Variant::VT_DOUBLE: {
		result = "VT_NUMBER";
	}
		break;
	case ug::Variant::VT_CSTRING: {
		result = "VT_CSTRING";
	}
		break;
	case ug::Variant::VT_STDSTRING: {
		result = "VT_STD_STRING";
	}
		break;
	case ug::Variant::VT_POINTER: {
		result = "VT_POINTER";
	}
		break;
	case ug::Variant::VT_CONST_POINTER: {
		result = "VT_CONST_POINTER";
	}
		break;
	case ug::Variant::VT_SMART_POINTER: {
		result = "VT_SMART_POINTER";
	}
		break;
	case ug::Variant::VT_CONST_SMART_POINTER: {
		result = "VT_CONST_SMART_POINTER";
	}
		break;
	default:
		break;
	}

	return result;
}

TypeAndArray paramClass2ParamType(JNIEnv *env, jobject obj) {

	TypeAndArray typeAndArray;

	int tmpResultType = ug::Variant::VT_INVALID;

	bool DEBUG = true;

	if (DEBUG) {
		std::cout << "trunk/ugbase/bindings/vrl/type_converter.cpp :"
				<< " paramClass2ParamType() " << std::endl;
		//std::cout << "true = " << true << std::endl;
		//std::cout << "false = " << false << std::endl;
	}

	std::string className = getClassName(env, obj);

	if (DEBUG) {
		std::cout << " paramClass2ParamType : className = " << className
				<< std::endl;

	}

	bool isArrayList = isJObjectAnArrayList(env, obj);

	//UG_LOG("paramClass2ParamType : isArrayList = " << std::boolalpha << isArrayList << std::endl);

	if (isArrayList) {
		UG_THROW(
				"Parameter is a java ArrayList, but we support only arrays." << "MSG from: trunk/ugbase/bindings/vrl/type_converter.cpp :" << " paramClass2ParamType() ");
		typeAndArray.isArray = false;
	}

	bool isArray = isJObjectAnArray(env, obj);
	/*UG_LOG("paramClass2ParamType : isArray = "<< isArray << std::endl);
	 UG_LOG("paramClass2ParamType : isArray = "<< std::boolalpha << isArray << std::endl);
	 */
	typeAndArray.isArray = isArray;

	if (isArray) {
		// if isArray cut the "object array" identifier off: "[L"

		className = className.substr(2, className.length());

		/*std::string lastSign = className.substr( className.length()-1, 1 );

		 // and cut off the tail if there is a semicolon ";"
		 if(lastSign == ";"){
		 className = className.substr(0, className.length()-1 );
		 }*/
		className = cutLastStringPart(className, ";");

		UG_LOG("paramClass2ParamType : isArray cutted className = "<< className << std::endl);
	}

	//bool tmp = className.compare("java.lang.Boolean");
	/*UG_LOG("paramClass2ParamType : className.compare(\"java.lang.Boolean\") = "<< tmp << std::endl);
	 UG_LOG("paramClass2ParamType : className.compare(\"java.lang.Boolean\") = "<< std::boolalpha << tmp << std::endl);
	 */

	if (className.compare("java.lang.Boolean") == 0) {
		tmpResultType = ug::Variant::VT_BOOL;
	} else if (className.compare("java.lang.Integer") == 0) {
		tmpResultType = ug::Variant::VT_INT;
	} else if (className.compare("java.lang.Double") == 0) {
		tmpResultType = ug::Variant::VT_DOUBLE;
	} else if (className.compare("java.lang.String") == 0) {
		tmpResultType = ug::Variant::VT_STDSTRING;
	} else if (className.compare("edu.gcsc.vrl.ug.Pointer") == 0) {
		tmpResultType = ug::Variant::VT_POINTER;
	} else if (className.compare("edu.gcsc.vrl.ug.SmartPointer") == 0) {
		tmpResultType = ug::Variant::VT_SMART_POINTER;
	}

	/*

	 if (className.compare("java.lang.Boolean") == true) {
	 std::cout << "paramClass2ParamType : if (className.compare(\"java.lang.Boolean\") == true)" << std::endl;

	 tmpResultType = ug::Variant::VT_BOOL;
	 } else if (className.compare("java.lang.Integer") == true) {
	 tmpResultType = ug::Variant::VT_INT;
	 } else if (className.compare("java.lang.Double") == true) {
	 tmpResultType = ug::Variant::VT_DOUBLE;
	 } else if (className.compare("java.lang.String") == true) {
	 tmpResultType = ug::Variant::VT_STDSTRING;
	 } else if (className.compare("edu.gcsc.vrl.ug.Pointer") == true) {
	 tmpResultType = ug::Variant::VT_POINTER;
	 } else if (className.compare("edu.gcsc.vrl.ug.SmartPointer") == true) {
	 tmpResultType = ug::Variant::VT_SMART_POINTER;
	 }
	 */

	// What about const pointer?
	// Answer: compare param types allows
	// non-const* to const* conversion
	// That is why we do not check that. We also use the same class for
	// const and non const pointer. Const checking is done via readOnly
	// bit of the pointer instance (Java wrapper).
	typeAndArray.type = tmpResultType;

	if (DEBUG) {
		std::cout << "ug::Variant::VT_INVALID = " << ug::Variant::VT_INVALID
				<< std::endl;
		std::cout << "paramClass2ParamType().tmpResultType = " << tmpResultType
				<< std::endl;
		std::cout << "VT_BOOL = " << ug::Variant::VT_BOOL << std::endl;
		std::cout << "VT_INT = " << ug::Variant::VT_INT << std::endl;
		std::cout << "VT_DOUBLE = " << ug::Variant::VT_DOUBLE << std::endl;
		std::cout << "VT_STDSTRING = " << ug::Variant::VT_STDSTRING
				<< std::endl;
		std::cout << "VT_POINTER = " << ug::Variant::VT_POINTER << std::endl;
		std::cout << "VT_SMART_POINTER = " << ug::Variant::VT_SMART_POINTER
				<< std::endl;
	}

	return typeAndArray;
}

bool compareParamTypes(JNIEnv *env, jobjectArray params,
		ug::bridge::Registry *reg,
		ug::bridge::ParameterInfo const& paramStack) {

	// christian poliwoda ug-log for debug
	UG_LOG(
			"trunk/ugbase/bindings/vrl/type_converter.cpp :: compareParamTypes(...)"<<std::endl);

	//#ifdef UG_DEBUG
	//	UG_LOG("\n -- BEGIN COMPARE --\n")
	//#endif

	// compare array lengths
	jsize len = env->GetArrayLength(params);

	if (len != paramStack.size()) {
		// christian poliwoda ug-log for debug start
		UG_LOG("type_converter.cpp :: if (len != paramStack.size())" << std::endl);
		UG_LOG("type_converter.cpp :: len  = "<< len << std::endl);
		UG_LOG("type_converter.cpp :: paramStack.size()  = "<< paramStack.size() << std::endl);

		// christian poliwoda ug-log for debug end

		return false;
	}

	// iterate over all param stack elements and compare their type with
	// the corresponding elements in the specified Java array
	for (size_t i = 0; i < (size_t) paramStack.size(); i++) {
		jobject param = env->GetObjectArrayElement(params, i);

		// we don't allow null values
		if (param == NULL) {
			std::stringstream ss;
			ss << "Value " << i << " == NULL!";

			jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
			env->ThrowNew(Exception, ss.str().c_str());
		}

		TypeAndArray paramType = paramClass2ParamType(env, param);

		// christian poliwoda ug-log for debug
		UG_LOG("type_converter.cpp :: compareParamTypes(...) : param = "<< param <<std::endl);
		UG_LOG("type_converter.cpp :: compareParamTypes(...) : uint paramType = paramClass2ParamType(env, param) = "<< paramType.type <<std::endl);
		UG_LOG("type_converter.cpp :: compareParamTypes(...) : paramStack.type(i) = "<< paramStack.type(i) <<std::endl);

		// allow non-const * to const *
		if (paramType.type == ug::Variant::VT_POINTER
				&& paramStack.type(i) == ug::Variant::VT_CONST_POINTER) {

			// christian poliwoda ug-log for debug
			UG_LOG(" compareParamTypes(...) if()1 => paramType == ug::Variant::VT_POINTER && paramStack.type(i) == ug::Variant::VT_CONST_POINTER" <<std::endl);

			paramType.type = ug::Variant::VT_CONST_POINTER;
		}

		// allow non-const * to const *
		if (paramType.type == ug::Variant::VT_SMART_POINTER
				&& paramStack.type(i) == ug::Variant::VT_CONST_SMART_POINTER) {

			// christian poliwoda ug-log for debug
			UG_LOG(
					" compareParamTypes(...) if()2 => paramType == ug::Variant::VT_SMART_POINTER && paramStack.type(i) == ug::Variant::VT_CONST_SMART_POINTER" <<std::endl);

			paramType.type = ug::Variant::VT_CONST_SMART_POINTER;
		}

		// allow std::string to c string
		if (paramType.type == ug::Variant::VT_STDSTRING
				&& paramStack.type(i) == ug::Variant::VT_CSTRING) {

			// christian poliwoda ug-log for debug
			UG_LOG(
					" compareParamTypes(...) if()3 => paramType == ug::Variant::VT_STDSTRING && paramStack.type(i) == ug::Variant::VT_CSTRING" <<std::endl);

			paramType.type = ug::Variant::VT_CSTRING;
		}

		// UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!)
		// allow non-const-smart* to non const*
		if (paramType.type == ug::Variant::VT_SMART_POINTER
				&& paramStack.type(i) == ug::Variant::VT_POINTER) {

			// christian poliwoda ug-log for debug
			UG_LOG(
					" compareParamTypes(...) if()4 = paramType == ug::Variant::VT_SMART_POINTER && paramStack.type(i) == ug::Variant::VT_POINTER" <<std::endl);

			paramType.type = ug::Variant::VT_POINTER;
		}

		// allow non-const-smart* to const*
		if (paramType.type == ug::Variant::VT_SMART_POINTER
				&& paramStack.type(i) == ug::Variant::VT_CONST_POINTER) {

			// christian poliwoda ug-log for debug
			UG_LOG(
					" compareParamTypes(...) if()5 => paramType == ug::Variant::VT_SMART_POINTER paramStack.type(i) == ug::Variant::VT_CONST_POINTER" <<std::endl);

			paramType.type = ug::Variant::VT_CONST_POINTER;
		}

		// allow const smart* to const*
		if (paramType.type == ug::Variant::VT_CONST_SMART_POINTER
				&& paramStack.type(i) == ug::Variant::VT_CONST_POINTER) {

			// christian poliwoda ug-log for debug
			UG_LOG(
					" compareParamTypes(...) if()6 => (paramType == ug::Variant::VT_CONST_SMART_POINTER && paramStack.type(i) == ug::Variant::VT_CONST_POINTER" <<std::endl);

			paramType.type = ug::Variant::VT_CONST_POINTER;
		}

		//new by christian poliwoda
		// allow integer to size_t
				if (paramType.type == ug::Variant::VT_INT
						&& paramStack.type(i) == ug::Variant::VT_SIZE_T) {

					// christian poliwoda ug-log for debug
					UG_LOG(" compareParamTypes(...) if()8 => (paramType.type == ug::Variant::VT_INT && paramStack.type(i) == ug::Variant::VT_SIZE_T" <<std::endl);

					paramType.type = ug::Variant::VT_SIZE_T;
				}

		if (paramType.type != (uint) paramStack.type(i)) {

			// christian poliwoda ug-log for debug
			UG_LOG(" compareParamTypes(...) if()7 => paramType != (uint) paramStack.type(i)" <<std::endl);

			// christian poliwoda ug-log for debug
			UG_LOG("compareParamTypes(...) if()7 : paramType = "<< paramType.type <<std::endl);
			UG_LOG("compareParamTypes(...) if()7 : (uint) paramStack.type(i) = "<< (uint) paramStack.type(i) <<std::endl);

			//#ifdef UG_DEBUG
			//			UG_LOG("requested by method:\n")
			//			printParamType(paramStack.type(i), i);
			//			UG_LOG("given as parameter:\n")
			//			printParamType(paramType, i);
			//#endif
			return false;
		}

		UG_LOG("after if()s "<<std::endl);

		//christian poliwoda debug
		//jclass paramClass = env->GetObjectClass(param);
		//jmethodID mid_getName = env->GetMethodID(paramClass, "getName", "()Ljava/lang/String;");
		//BOOM //jstring paramClassName = (jstring) env->CallObjectMethod(paramClass, mid_getName);
		//UG_LOG("paramClassName = "<< stringJ2C(env, paramClassName )<<std::endl);

		UG_LOG(
				"getParamClassName(env, param) = "<< getParamClassName(env, param) <<std::endl);

		// check if param is assignable
		const ug::bridge::ClassNameNode* classNameNode =
				ug::vrl::invocation::getClassNodePtrByName(reg,
						getParamClassName(env, param));

		// christian poliwoda ug-log for debug
		//UG_LOG(	" compareParamTypes(...) classNameNode->name() = "<< classNameNode->name() <<std::endl);

		if (classNameNode != NULL) {
			UG_LOG("if (classNameNode != NULL) "<<std::endl);
			UG_LOG(
					" compareParamTypes(...) classNameNode->name() = "<< classNameNode->name() <<std::endl);

			if (paramStack.class_name(i) != NULL) {
				UG_LOG("if (paramStack.class_name(i) != NULL) "<<std::endl);
				UG_LOG(
						" paramStack.class_name(i) = "<< paramStack.class_name(i) <<std::endl);

				if (!ug::bridge::ClassNameTreeContains(*classNameNode,
						paramStack.class_name(i))) {
					UG_LOG(
							"if (!ug::bridge::ClassNameTreeContains(*classNameNode, paramStack.class_name(i))) "<<std::endl);
					UG_LOG(" return false "<<std::endl);

					return false;
				}
			}
		}		//if (classNameNode != NULL)
				// else for debug by christian poliwoda
		else {
			UG_LOG(" compareParamTypes(...) classNameNode == NULL "<<std::endl);
		}
	}

	//#ifdef UG_DEBUG
	//	UG_LOG(" -- ALL TRUE --\n" << std::endl)
	//#endif

	return true;
}

std::string getParamTypesAsString(JNIEnv *env, jobjectArray const& array) {

	std::stringstream ss;

	size_t length = (size_t) env->GetArrayLength(array);

	//iterate through the parameter array and add the type name to
	// the stream
	for (size_t i = 0; i < length; ++i) {
		jobject value = env->GetObjectArrayElement(array, i);

		TypeAndArray java_value_type = paramClass2ParamType(env, value);

		if (i > 0) {
			ss << ", ";
		}
		bool classValue = java_value_type.type == ug::Variant::VT_CONST_POINTER
				|| java_value_type.type == ug::Variant::VT_POINTER
				|| java_value_type.type == ug::Variant::VT_SMART_POINTER
				|| java_value_type.type == ug::Variant::VT_CONST_SMART_POINTER;

		if (classValue) {
			ss << getParamClassName(env, value);
		} else {
			ss << getParamTypeAsString(java_value_type.type);
		}

	}

	return ss.str();
}

//original Methode
//the modified one can be find below
//christian poliwoda
void jobjectArray2ParamStack(JNIEnv *env, ug::bridge::Registry* reg,
		ug::bridge::ParameterStack& paramsOut,
		const ug::bridge::ParameterInfo& paramsTemplate,
		jobjectArray const& array) {
	using namespace ug::bridge;

	UG_LOG("trunk/ugbase/bindings/vrl/type_converter.cpp ::  jobjectArray2ParamStack( )"<<std::endl);

	//	iterate through the parameter list and copy the value in the
	//  associated stack entry.
	for (size_t i = 0; i < (size_t) paramsTemplate.size(); ++i) {

		uint template_value_Type = paramsTemplate.type(i);

		jobject value = env->GetObjectArrayElement(array, i);

		// only used for
		// UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!)
		//uint java_value_type = paramClass2ParamType(env, value);    //original before TypeAndArray
		uint java_value_type = paramClass2ParamType(env, value).type; //with TypeAndArray

		// we don't allow null values
		if (value == NULL) {
			std::stringstream ss;
			ss << "Value " << i << " == NULL!";

			jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
			env->ThrowNew(Exception, ss.str().c_str());
		}

		//debug
		UG_LOG("template_value_Type = "<< template_value_Type <<std::endl);

		switch (template_value_Type) {
		case ug::Variant::VT_BOOL: {

			paramsOut.push<bool>(jObject2Boolean(env, value));

		}
			break;
		case ug::Variant::VT_INT: {
			paramsOut.push<int>(jObject2Int(env, value));
		}
			break;
		case ug::Variant::VT_SIZE_T: {
			paramsOut.push<size_t>((size_t) jObject2Int(env, value));
		}
			break;
		case ug::Variant::VT_FLOAT:
		case ug::Variant::VT_DOUBLE: {
			paramsOut.push<number>(jObject2Double(env, value));
		}
			break;
		case ug::Variant::VT_CSTRING: {
			paramsOut.push(jObject2String(env, value));
		}
			break;

		case ug::Variant::VT_STDSTRING: {
			paramsOut.push(jObject2String(env, value));
		}
			break;

		case ug::Variant::VT_POINTER: {
			const ug::bridge::ClassNameNode* node =
					ug::vrl::invocation::getClassNodePtrByName(reg,
							jPointerGetName(env, value));
			paramsOut.push(jObject2Pointer(env, value), node);
		}
			break;
		case ug::Variant::VT_CONST_POINTER: {
			const ug::bridge::ClassNameNode* node =
					ug::vrl::invocation::getClassNodePtrByName(reg,
							jPointerGetName(env, value));

			// UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!)
			if (java_value_type == ug::Variant::VT_CONST_SMART_POINTER) {
				paramsOut.push(
						(void*) jObject2ConstSmartPointer(env, value).get(),
						node);
			} else if (java_value_type == ug::Variant::VT_SMART_POINTER) {
				paramsOut.push((void*) jObject2SmartPointer(env, value).get(),
						node);
			} else {
				paramsOut.push(jObject2Pointer(env, value), node);
			}

			//paramsOut.push_const_pointer(
			//		jObject2Pointer(env, value), node);
		}
			break;
		case ug::Variant::VT_SMART_POINTER: {
			const ug::bridge::ClassNameNode* node =
					ug::vrl::invocation::getClassNodePtrByName(reg,
							jPointerGetName(env, value));

			paramsOut.push(jObject2SmartPointer(env, value), node);
			//				}
		}
			break;
		case ug::Variant::VT_CONST_SMART_POINTER: {
			const ug::bridge::ClassNameNode* node =
					ug::vrl::invocation::getClassNodePtrByName(reg,
							jPointerGetName(env, value));

			paramsOut.push(jObject2ConstSmartPointer(env, value), node);
		}
			break;
		}

	} // end for
}

jobject param2JObject(JNIEnv *env, ug::bridge::ParameterStack& params,
		size_t index) {

	// christian poliwoda ug-log for debug
	UG_LOG("trunk/ugbase/bindings/vrl/type_converter.cpp : param2JObject(...)"<<std::endl);

	using namespace ug::bridge;
	//	iterate through the parameter list and copy the value in the
	//	associated stack entry.
	int type = params.type(index);

	switch (type) {
	case ug::Variant::VT_BOOL: {
		//christian poliwoda
		//QUESTION: params is Vector check ??
		UG_LOG("param2JObject(...) : = switch case VT_BOOL"<<std::endl);

		bool isVector = params.is_vector(index);
		UG_LOG(" params.is_vector( "<< index << ") = "<< isVector <<std::endl);
		if (isVector) {
			UG_LOG(
					" param is vector -> calling smartPointer2JObject( env, params.to<SmartPtr<void> >(index) )" << std::endl);
			return smartPointer2JObject(env, params.to<SmartPtr<void> >(index));
		}

		return boolean2JObject(env, params.to<bool>(index));
	}
		break;
	case ug::Variant::VT_INT: {
		UG_LOG("param2JObject(...) : = switch case VT_INT"<<std::endl);
		return int2JObject(env, params.to<int>(index));
	}
		break;
	case ug::Variant::VT_SIZE_T: {
		UG_LOG("param2JObject(...) : = switch case VT_SIZE_T"<<std::endl);
		return int2JObject(env, (int) params.to<size_t>(index));
	}
		break;
	case ug::Variant::VT_FLOAT:
	case ug::Variant::VT_DOUBLE: {
		return double2JObject(env, params.to<number>(index));
	}
		break;
	case ug::Variant::VT_CSTRING: {
		return string2JObject(env, params.to<const char*>(index));
	}
		break;
	case ug::Variant::VT_STDSTRING: {
		return string2JObject(env, params.to<const std::string&>(index).c_str());
	}
		break;
	case ug::Variant::VT_POINTER: {
		return pointer2JObject(env, params.to<void*>(index));
	}
		break;
	case ug::Variant::VT_CONST_POINTER: {
		return pointer2JObject(env, (void*) params.to<const void*>(index));
	}
		break;
	case ug::Variant::VT_SMART_POINTER: {
		return smartPointer2JObject(env, params.to<SmartPtr<void> >(index));
	}
		break;
	case ug::Variant::VT_CONST_SMART_POINTER: {
		return constSmartPointer2JObject(env,
				params.to<ConstSmartPtr<void> >(index));
	}
		break;
	}

	return jobject();
}

//	added by Christian Poliwoda
//	christian.poliwoda@gcsc.uni-frankfurt.de
//	y 13 m 07 d 10
/**
 * This method is designed to compare the first characters (minimum length)
 * of two strings.
 */
bool startsWithSymbolOrderMin(JNIEnv *env, std::string check,
		std::string symbolOrder) {

	bool result = false;

	/*bool DEBUG = true;
	 if (DEBUG) {
	 std::cout << "startsWithSymbolOrderMin() STD::string check = " << check
	 << " , symbolOrder = " << symbolOrder << std::endl;
	 }*/

	size_t checkLenght = check.size();
	size_t symbolOrderLenght = symbolOrder.size();
	size_t minLenght = 0;

	if (symbolOrderLenght <= checkLenght) {
		minLenght = symbolOrderLenght;
	} else {
		minLenght = checkLenght;
	}

	std::string firstSignsCheck = check.substr(0, minLenght);
	std::string firstSignsSymbolOrder = symbolOrder.substr(0, minLenght);

	/*if (DEBUG) {
	 std::cout << "firstSignsCheck: " << firstSignsCheck << std::endl;
	 std::cout << "firstSignsSymbolOrder: " << firstSignsSymbolOrder
	 << std::endl;
	 }*/

	if (firstSignsCheck == firstSignsSymbolOrder) {
		result = true;
	}

	return result;
}

//	added by Christian Poliwoda
//	christian.poliwoda@gcsc.uni-frankfurt.de
//	y 13 m 06 d 06
/**
 * This method is designed to compare the first characters (minimum length)
 * of two strings.
 */
bool startsWithSymbolOrderMin(JNIEnv *env, jstring check, jstring symbolOrder) {
	bool result = false;

	/*bool DEBUG = true;
	 if (DEBUG) {
	 std::cout << "startsWithSymbolOrderMin() Jstring check = " << check
	 << " , symbolOrder = " << symbolOrder << std::endl;
	 }*/

	std::string stdStringCheck = stringJ2C(env, check);
	std::string stdStringSymbolOrder = stringJ2C(env, symbolOrder);

	result = startsWithSymbolOrderMin(env, stdStringCheck,
			stdStringSymbolOrder);

	/*
	 size_t checkLenght = stdStringCheck.size();
	 size_t symbolOrderLenght = stdStringSymbolOrder.size();
	 size_t minLenght = 0;

	 if (symbolOrderLenght <= checkLenght) {
	 minLenght = symbolOrderLenght;
	 } else {
	 minLenght = checkLenght;
	 }

	 std::string firstSignsCheck = stdStringCheck.substr(0, minLenght);
	 std::string firstSignsSymbolOrder = stdStringSymbolOrder.substr(0,
	 minLenght);

	 std::cout << "firstSignsCheck: " << firstSignsCheck << std::endl;
	 std::cout << "firstSignsSymbolOrder: " << firstSignsSymbolOrder
	 << std::endl;

	 if (firstSignsCheck == firstSignsSymbolOrder) {
	 result = true;
	 }
	 */

	return result;
}

//	added by Christian Poliwoda
//	christian.poliwoda@gcsc.uni-frankfurt.de
//	y 13 m 05 d 28
bool isJObjectAnArray(JNIEnv *env, jobject value) {

	/*bool DEBUG = true;
	 if (DEBUG) {
	 std::cout << "trunk/ugbase/bindings/vrl/type_converter.cpp :"
	 << " isJObjectAnArray() " << std::endl;
	 }*/

//	jstring valueClassName = stringC2J(env, getClassName(env, value).c_str());
//	jstring braked = stringC2J(env, "[");
	std::string valueClassName = getClassName(env, value);
	std::string braked = "[";

	/*if (DEBUG) {
	 std::cout << " valueClassName = " << valueClassName << std::endl;
	 std::cout << " braked = " << braked << std::endl;
	 }*/

// if the className starts with [ these is a hint
// that the value-object contains an array
// these is the reason because of jni declaration
// of types and classes
	return startsWithSymbolOrderMin(env, valueClassName, braked);
}

//	added by Christian Poliwoda
//	christian.poliwoda@gcsc.uni-frankfurt.de
//	y 13 m 10 d 30
bool isJObjectAnArrayList(JNIEnv *env, jobject value) {

	/*bool DEBUG = true;

	 if (DEBUG) {
	 std::cout << "trunk/ugbase/bindings/vrl/type_converter.cpp :"
	 << " isJObjectAnArrayList() " << std::endl;
	 }*/

	std::string valueClassName = getClassName(env, value);
	std::string startsWith = "java.util.ArrayList";

	/*if (DEBUG) {
	 std::cout << " valueClassName = " << valueClassName << std::endl;
	 std::cout << " startsWith = " << startsWith << std::endl;
	 }*/

// if the className starts with java.util.ArrayList these is a hint
// that the value-object contains an arraylist
	return startsWithSymbolOrderMin(env, valueClassName, startsWith);
}

//	added by Christian Poliwoda
//	christian.poliwoda@gcsc.uni-frankfurt.de
//	y 13 m 10 d 31
std::string cutLastStringPart(std::string original, std::string toCut) {

	//UG_LOG("trunk/ugbase/bindings/vrl/type_converter.cpp :" << " cutLastStringPart() " << std::endl);

	int cutLength = toCut.length();
	int originalLength = original.length();

	//check if the string that should be cutted is bigger than the removing part
	if (cutLength >= originalLength) {
		return "";
	}

	//check if the to be cutted string contains at the end the part that
	//should be removed
	std::string lastPart = original.substr(originalLength - cutLength,
			cutLength);

	if (toCut == lastPart) {

		return original.substr(0, originalLength - cutLength);
	}

	return original;
}

//	added by Christian Poliwoda
//	christian.poliwoda@gcsc.uni-frankfurt.de
//	y 13 m 05 d 28
jbooleanArray jObject2BooleanArray(JNIEnv *env, jobject object) {

	bool DEBUG = true;
	if (DEBUG) {
		std::cout << "trunk/ugbase/bindings/vrl/type_converter.cpp :"
				<< " jObject2BooleanArray() " << std::endl;
	}

	jbooleanArray newJBoolArray = NULL;

	//jclass argClass = env->GetObjectClass(object);
	//jstring jStringClassName = getClassName(env, argClass);
	//jstring jStringClassName = getClassName(env, object);
	//std::string stdStringClassName = stringJ2C(env, jStringClassName);
	//std::string stdStringClassName = getClassName(env, argClass);

	std::string stdStringClassName = getClassName(env, object);
	std::string stdStringParamClassName = getParamClassName(env, object);

	/*
	 jclass argClass = env->GetObjectClass(object);
	 jstring jStringClassName = getClassName(env, argClass);
	 std::string stdStringClassName = stringJ2C(env, jStringClassName);
	 */

	if (DEBUG) {

		std::cout << " stdStringClassName = " << stdStringClassName
				<< std::endl;
	}

	std::string stdObjectArrayBeginsWith = "[L";

	if (DEBUG) {
		std::cout << " stdObjectArrayBeginsWith = " << stdObjectArrayBeginsWith
				<< std::endl;
	}

	if (startsWithSymbolOrderMin(env, stdStringClassName,
			stdObjectArrayBeginsWith)) {

		if (DEBUG) {
			std::cout
					<< " if (startsWithSymbolOrderMin(env, stdStringClassName, stdObjectArrayBeginsWith)) {"
					<< std::endl;
		}

		size_t stdStringClassNameLenght = stdStringClassName.size();
		size_t objectArrayLenght = stdObjectArrayBeginsWith.size();

		std::string arrayElementClassName = stdStringClassName.substr(
				objectArrayLenght, stdStringClassNameLenght);

		arrayElementClassName = cutLastStringPart(arrayElementClassName, ";");

		if (DEBUG) {
			std::cout << "  object is an array of  (arrayElementClassName =) "
					<< arrayElementClassName << std::endl;
		}

		//cast to object array because we now know that it is an array

		jobjectArray objectArray = (jobjectArray) object;

		jsize arrayLenght = env->GetArrayLength(objectArray);

		if (DEBUG) {
			std::cout << "  arrayLenght = " << arrayLenght << std::endl;
		}

		// create / set a new array for return
		newJBoolArray = env->NewBooleanArray(arrayLenght);
		if (DEBUG) {
			std::cout << "  newJBoolArray = " << newJBoolArray << std::endl;
		}

		// COPY VALUES INTO newJBoolArray
		if (DEBUG) {
			std::cout << "SETting / COPYing VALUES INTO newJBoolArray "
					<< std::endl;
		}

		jboolean *elementsOfNewJBoolArray = env->GetBooleanArrayElements(
				newJBoolArray, false);

		if (DEBUG) {
			std::cout << "  empty booleanArray:" << std::endl;
			for (int i = 0; i < arrayLenght; i++) {
				std::cout << "(bool) elementsOfNewJBoolArray[ " << i << " ] = "
						<< (bool) elementsOfNewJBoolArray[i] << std::endl;
			}

			std::cout << "setting values:" << std::endl;
		}

		/*
		 // howto set a const value
		 for (int i = 0; i < arrayLenght; i++) {
		 elementsOfNewJBoolArray[i] = JNI_TRUE;
		 }
		 */

		//setting/copying values from an other array
		jobject tmpObject = NULL;

		for (int i = 0; i < arrayLenght; i++) {
			//get object from original Array
			tmpObject = env->GetObjectArrayElement(objectArray, i);
			// cast object to jboolean and store the value in the new array
			elementsOfNewJBoolArray[i] = jObject2Boolean(env, tmpObject);
		}

		if (DEBUG) {
			std::cout << "  filled booleanArray:" << std::endl;

			for (int i = 0; i < arrayLenght; i++) {
				std::cout << "(bool) elements[ " << i << " ] = "
						<< (bool) elementsOfNewJBoolArray[i] << std::endl;
			}
		}

		env->ReleaseBooleanArrayElements(newJBoolArray, elementsOfNewJBoolArray,
				0);

	} else {
		std::cout << "trunk/ugbase/bindings/vrl/type_converter.cpp :"
				<< " jObject2BooleanArray()" << std::endl;
		std::cout << "object is NOT an array " << std::endl;
	}

	return newJBoolArray;
}

/*

 //	added by Christian Poliwoda
 //	christian.poliwoda@gcsc.uni-frankfurt.de
 //	y 13 m 05 d 08
 // TODO
 // needed to be merged with or replace original method jobjectArray2ParamStack(...)
 // if finished to reduce code duplication
 // PONDON TO PushLuaStackEntryToParamStack(
 // ParameterStack& ps, lua_State* L, int index, bool bIsVector)
 // in bindings_lua.cpp
 void jobjectArray2ParamStack(JNIEnv *env, ug::bridge::Registry* reg,
 ug::bridge::ParameterStack& paramsOut,
 const ug::bridge::ParameterInfo& paramsTemplate,
 jobjectArray const& array) {

 bool DEBUG = true;

 if (DEBUG) {
 //std::cout << " HALLO MIHI " << std::endl;
 std::cout << "trunk/ugbase/bindings/vrl/type_converter.cpp :"
 << " jobjectArray2ParamStack[NEW]() " << std::endl;
 //std::cout << " DEBUG = true = " << true << std::endl;
 //std::cout << " DEBUG = false = " << false << std::endl;
 }

 // // if no information is available about bIsVector or bIsVector is FALSE
 // // the origin method is called
 // if(!bIsVector){
 //	jobjectArray2ParamStack(env, reg, paramsOut, paramsTemplate, array);
 // } else {
 using namespace ug::bridge;

 //	iterate through the parameter list and copy the value in the
 //  associated stack entry.
 for (size_t i = 0; i < (size_t) paramsTemplate.size(); ++i) {

 uint template_value_Type = paramsTemplate.type(i);

 jobject value = env->GetObjectArrayElement(array, i);

 // only used for
 // UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!)
 TypeAndArray java_value_type = paramClass2ParamType(env, value);

 if (DEBUG) {
 std::cout << "jobjectArray2ParamStack() :" << std::endl;
 std::cout << " for (size_t i = " << i << "; i < "
 << (size_t) paramsTemplate.size() << "; ++i) " << std::endl;
 std::cout
 << " fuer switch-case verantwortlich => uint template_value_Type = paramsTemplate.type(i) = "
 << template_value_Type
 << " = paramsTemplate.class_name(i) = "
 << paramsTemplate.class_name(i) << std::endl;
 std::cout
 << " UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!) => TypeAndArray java_value_type = paramClass2ParamType(env, value) = "
 << java_value_type.type << std::endl;
 }

 // we don't allow null values
 if (value == NULL) {
 std::stringstream ss;
 ss << "Value " << i << " == NULL!";

 jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
 env->ThrowNew(Exception, ss.str().c_str());
 }

 // check if value is an array / a vector
 // if no information is available about isVector or
 // isVector is FALSE
 // the origin method is called

 bool isVector = isJObjectAnArray(env, value);

 if (DEBUG) {
 std::cout << " isVector = " << isVector << std::endl;
 }

 switch (template_value_Type) {
 case ug::Variant::VT_BOOL: {

 if (DEBUG) {
 std::cout << " in switch VT_BOOL" << std::endl;
 }

 // if isVector is FALSE
 // the origin method is called

 if (!isVector) {

 //  original
 paramsOut.push<bool>(jObject2Boolean(env, value));
 } else {

 if (DEBUG) {
 std::cout << "true = " << true << std::endl;
 std::cout << "false = " << false << std::endl;
 std::cout << "!isVector = " << !isVector << std::endl;
 std::cout << " calling jObject2BooleanArray()" << std::endl;
 }

 jbooleanArray jboolArray = jObject2BooleanArray(env, value);

 jsize arrayLength = env->GetArrayLength(jboolArray);

 //create vector with needed length
 std::vector<bool> boolVector(arrayLength);

 std::cout << "std::vector<bool> boolVector(arrayLength);" << std::endl;

 //
 //fill the std::vector<bool> with the values of the jbooleanArray
 //
 jboolean *elements = env->GetBooleanArrayElements(jboolArray,
 false);

 for (int i = 0; i < arrayLength; ++i) {
 boolVector[i] = (bool) elements[i];
 }

 std::cout << "for (int i = 0; i < arrayLength; ++i)" << std::endl;

 //SO_EIN_SCHEISS__WIE_GEHT_DAS_MIT_SMARTptr_???

 SmartPtr<std::vector<bool> > smartPointer =
 ////sollte sein / kompeliert aber NICHT
 //SmartPtr<std::vector<bool> >(boolVector);
 //
 ////kompiliert, hat aber keine daten
 //SmartPtr<std::vector<bool> >(new std::vector<bool>());
 //
 ////VRL-Studio: JavaApplicationStub(36579,0x130883000) malloc: *** error for object 0x12303ec80: pointer being freed was not allocated
 //SmartPtr<std::vector<bool> >(&boolVector);
 //
 //kein speicher fehler -> ? funktioniert ???
 SmartPtr<std::vector<bool> >(new std::vector<bool>());

 //
 //return the smartPointer with the boolean-vector
 //
 UG_LOG( "trunk/ugbase/bindings/vrl/type_converter.cpp :"
 << " jobjectArray2ParamStack[NEW]() "<<
 "BEFORE paramsOut.push(smartPointer);" << std::endl);

 paramsOut.push(smartPointer);

 UG_LOG( "trunk/ugbase/bindings/vrl/type_converter.cpp :"
 << " jobjectArray2ParamStack[NEW]() "<<
 "AFTER paramsOut.push(smartPointer);" << std::endl);
 //HIER_WEITER__GEHT_DAS_SO
 //paramsOut.push<bool[]>(jboolArray);
 //paramsOut.push<SmartPtr<std::vector<bool>>>(jboolArray);

 //paramsOut.push<SmartPtr<std::vector<bool>>>(
 //jObject2BooleanArray(env, value));

 } //END else if (!isVector)
 }
 break;
 case ug::Variant::VT_INT: {

 if (DEBUG) {
 std::cout << " in switch VT_INT" << std::endl;
 }

 paramsOut.push<int>(jObject2Int(env, value));
 }
 break;
 case ug::Variant::VT_SIZE_T: {

 if (DEBUG) {
 std::cout << " in switch VT_SIZE_T" << std::endl;
 }

 paramsOut.push<size_t>((size_t) jObject2Int(env, value));
 }
 break;
 case ug::Variant::VT_FLOAT:
 case ug::Variant::VT_DOUBLE: {
 if (DEBUG) {
 std::cout << " in switch VT_FLOAT / VT_DOUBLE" << std::endl;
 }
 paramsOut.push<number>(jObject2Double(env, value));
 }
 break;
 case ug::Variant::VT_CSTRING: {

 if (DEBUG) {
 std::cout << " in switch VT_CSTRING" << std::endl;
 }
 paramsOut.push(jObject2String(env, value));
 }
 break;

 case ug::Variant::VT_STDSTRING: {
 if (DEBUG) {
 std::cout << " in switch VT_STDSTRING" << std::endl;
 }
 paramsOut.push(jObject2String(env, value));
 }
 break;

 case ug::Variant::VT_POINTER: {

 if (DEBUG) {
 std::cout << " in switch VT_POINTER" << std::endl;
 }

 const ug::bridge::ClassNameNode* node =
 ug::vrl::invocation::getClassNodePtrByName(reg,
 jPointerGetName(env, value));
 paramsOut.push(jObject2Pointer(env, value), node);
 }
 break;
 case ug::Variant::VT_CONST_POINTER: {

 if (DEBUG) {
 std::cout << " in switch VT_CONST_POINTER" << std::endl;
 }

 const ug::bridge::ClassNameNode* node =
 ug::vrl::invocation::getClassNodePtrByName(reg,
 jPointerGetName(env, value));

 if (DEBUG) {
 std::cout << " in switch VT_CONST_POINTER : node = " << node
 << std::endl;
 }

 // UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!)
 if (java_value_type.type == ug::Variant::VT_CONST_SMART_POINTER) {

 if (DEBUG) {
 std::cout
 << " UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!)"
 << std::endl;
 std::cout
 << " in SUB_IF => java_value_type == ug::Variant::VT_CONST_SMART_POINTER"
 << std::endl;
 std::cout << " java_value_type.type = "
 << java_value_type.type << std::endl;
 std::cout << " ug::Variant::VT_CONST_SMART_POINTER = "
 << ug::Variant::VT_CONST_SMART_POINTER << std::endl;
 }

 paramsOut.push(
 (void*) jObject2ConstSmartPointer(env, value).get(),
 node);

 } else if (java_value_type.type == ug::Variant::VT_SMART_POINTER) {

 if (DEBUG) {
 std::cout
 << " UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!)"
 << std::endl;
 std::cout
 << " in SUB_ELSE_IF => java_value_type == ug::Variant::VT_SMART_POINTER"
 << std::endl;
 std::cout << " java_value_type.type = "
 << java_value_type.type << std::endl;
 std::cout << " ug::Variant::VT_SMART_POINTER = "
 << ug::Variant::VT_SMART_POINTER << std::endl;
 }


 //original start
 // paramsOut.push((void*) jObject2SmartPointer(env, value).get(), node);
 //original ende



 //modification01 start by christian poliwoda
 // paramsOut.push(jObject2SmartPointer(env, value).get(), node);
 //modification01 end


 //modification02 start by christian poliwoda
 paramsOut.push(jObject2SmartPointer(env, value).get(), node);
 //modification02 end

 if (DEBUG) {
 std::cout
 << " AFTER paramsOut.push( jObject2SmartPointer() )"
 << std::endl;
 }

 } else {

 if (DEBUG) {
 std::cout << " in SUB_ELSE " << std::endl;
 }

 paramsOut.push(jObject2Pointer(env, value), node);
 }

 //						paramsOut.push_const_pointer(
 //								jObject2Pointer(env, value), node);
 }
 break;
 case ug::Variant::VT_SMART_POINTER: {

 if (DEBUG) {
 std::cout << " in switch VT_SMART_POINTER" << std::endl;
 }

 const ug::bridge::ClassNameNode* node =
 ug::vrl::invocation::getClassNodePtrByName(reg,
 jPointerGetName(env, value));

 paramsOut.push(jObject2SmartPointer(env, value), node);
 //				}
 }
 break;
 case ug::Variant::VT_CONST_SMART_POINTER: {

 if (DEBUG) {
 std::cout << " in switch VT_CONST_SMART_POINTER" << std::endl;
 }

 const ug::bridge::ClassNameNode* node =
 ug::vrl::invocation::getClassNodePtrByName(reg,
 jPointerGetName(env, value));

 paramsOut.push(jObject2ConstSmartPointer(env, value), node);
 }
 break;
 } //switch(template_vale_type)

 } // end for
 //} // end else (!bIsVector)
 }

 */

/*

 //	added by Christian Poliwoda
 //	christian.poliwoda@gcsc.uni-frankfurt.de
 //	y 13 m 05 d xx
 // TODO
 // needed to be merged with original method param2JObject(env, params, index)
 // if finished to reduce code duplication
 jobject param2JObject(JNIEnv *env, ug::bridge::ParameterStack& params,
 size_t index, bool bIsVector) {

 // if no information is available about bIsVector or bIsVector is FALSE
 // the origin method is called
 if (!bIsVector) {
 return param2JObject(env, params, index);
 } else {
 using namespace ug::bridge;
 //	iterate through the parameter list and copy the value in the
 //	associated stack entry.
 int type = params.type(index);

 switch (type) {
 case ug::Variant::VT_BOOL: {
 //return boolean2JObject(env, params.to<bool>(index), bIsVector);

 }
 break;
 case ug::Variant::VT_INT: {
 //return int2JObject(env, params.to<int>(index), bIsVector);
 }
 break;
 case ug::Variant::VT_SIZE_T: {
 //return int2JObject(env, (int) params.to<size_t>(index), bIsVector);
 }
 break;
 case ug::Variant::VT_FLOAT:
 case ug::Variant::VT_DOUBLE: {
 //return double2JObject(env, params.to<number>(index), bIsVector);
 }
 break;
 case ug::Variant::VT_CSTRING: {
 //return string2JObject(env, params.to<const char*>(index), bIsVector);
 }
 break;
 case ug::Variant::VT_STDSTRING: {
 //return string2JObject(env,
 //	params.to<const std::string&>(index).c_str(), bIsVector);
 }
 break;
 case ug::Variant::VT_POINTER: {
 //return pointer2JObject(env, params.to<void*>(index), bIsVector);
 }
 break;
 case ug::Variant::VT_CONST_POINTER: {
 //return pointer2JObject(env, (void*) params.to<const void*>(index),
 //		bIsVector);
 }
 break;
 case ug::Variant::VT_SMART_POINTER: {
 //return smartPointer2JObject(env, params.to<SmartPtr<void> >(index),
 //		bIsVector);
 }
 break;
 case ug::Variant::VT_CONST_SMART_POINTER: {
 //return constSmartPointer2JObject(env,
 //		params.to<ConstSmartPtr<void> >(index), bIsVector);
 }
 break;
 }

 return jobject();
 }

 }
 */

/*
 //	added by Christian Poliwoda
 //	christian.poliwoda@gcsc.uni-frankfurt.de
 //	y13 m05 d07
 // TODO
 // needed to be merged with original method
 // if finished to reduce code duplication
 jobject boolean2JObject(JNIEnv *env, jboolean value, bool bIsVector) {

 // if no information is available about bIsVector or bIsVector is FALSE
 // the origin method is called
 if (!bIsVector) {
 return boolean2JObject(env, value);
 } else {
 if (ggg) {
 SmartPtr<std::vector<bool> > spVec = SmartPtr<std::vector<bool> >(
 new std::vector<bool>());
 //	lua_pushnil(L);
 }
 }
 }
 */

/*
 //	added by Christian Poliwoda
 //	christian.poliwoda@gcsc.uni-frankfurt.de
 //	y 13 m 05 d 28
 SmartPtr<std::vector<bool> > jObject2BooleanVector(JNIEnv *env, jobject obj) {

 SmartPtr<std::vector<bool> > spVec = SmartPtr<std::vector<bool> >(
 new std::vector<bool>());

 int valueIsArraySize = env->GetArrayLength(obj);

 for (int i = 0; i < valueIsArraySize; ++i) {
 spVec->push_back((bool) valueIsArray[i]);
 }

 paramsOut.push(spVec);

 return spVec;
 }*/

int paramType2Int(const ug::bridge::ParameterInfo& params, size_t index) {
	using namespace ug::bridge;

	int type = params.type(index);

	return type;
}

jobjectArray params2NativeParams(JNIEnv *env,
		const ug::bridge::ExportedFunctionBase& func) {

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeParamInfo");

	jobjectArray result = env->NewObjectArray(func.num_parameter(), cls, 0);

	const ug::bridge::ParameterInfo& params = func.params_in();

	for (size_t i = 0; i < func.num_parameter(); i++) {

		// create instance
		jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
		jobject obj = env->NewObject(cls, methodID);

		//assign values

		jmethodID setType = env->GetMethodID(cls, "setType", "(I)V");
		jmethodID setID = env->GetMethodID(cls, "setId", "(I)V");
		jmethodID setClassName = env->GetMethodID(cls, "setClassName",
				"(Ljava/lang/String;)V");
		jmethodID setClassNames = env->GetMethodID(cls, "setClassNames",
				"([Ljava/lang/String;)V");
		//	jmethodID setHelp = env->GetMethodID(cls,
		//		"setHelp", "(Ljava/lang/String;)V");
		//	jmethodID setToolTip = env->GetMethodID(cls,
		//		"setTooltip", "(Ljava/lang/String;)V");
		jmethodID setParamInfo = env->GetMethodID(cls, "setParamInfo",
				"([Ljava/lang/String;)V");

		// christian poliwoda start1
		jmethodID setParamAVector = env->GetMethodID(cls, "setParamAVector",
				"(Ljava/lang/Boolean;)V");
		// christian poliwoda end1

		// check for emptyness
		bool pointerType = params.type(i) == ug::Variant::VT_CONST_POINTER
				|| params.type(i) == ug::Variant::VT_POINTER
				|| params.type(i) == ug::Variant::VT_SMART_POINTER
				|| params.type(i) == ug::Variant::VT_CONST_SMART_POINTER;

		if (pointerType && strlen(params.class_name(i)) == 0) {
			std::cerr << func.name() << ", param(" << i << ")==EMPTY"
					<< std::endl;
			exit(1);
		}

		using namespace ug::bridge;

		// TODO Unfortunately we don't know a better way to convert an enumeration
		// from C++ to Java. Currently we just use integers :(
		int type = paramType2Int(params, i);

		env->CallVoidMethod(obj, setType, type);
		env->CallVoidMethod(obj, setID, i);

		env->CallVoidMethod(obj, setClassName,
				stringC2J(env, params.class_name(i)));

		env->CallVoidMethod(obj, setClassNames,
				stringArrayC2J(env,
						getBaseClassNames(params.class_name_node(i))));

		env->CallVoidMethod(obj, setParamInfo,
				stringArrayC2J(env, func.parameter_info_vec(i)));

		// christian poliwoda start2

		// check if parameter is vector of "basic-parameters"
		// UG_LOG(" trunk/ugbase/bindings/vrl/type_converter.cpp : params2NativeParams() . CallVoidMethod(obj, setParamAVector, isVector)");
		jobject isVector = boolean2JObject(env, params.is_vector(i));

		env->CallVoidMethod(obj, setParamAVector, isVector);

		// christian poliwoda end2

		// set array element
		env->SetObjectArrayElement(result, i, obj);
	} //for (size_t i = 0; i < func.num_parameter(); i++) {

	return result;
}

jobjectArray params2NativeParams(JNIEnv *env,
		const ug::bridge::ExportedConstructor& constructor) {

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeParamInfo");

	jobjectArray result = env->NewObjectArray(constructor.num_parameter(), cls,
			0);

	const ug::bridge::ParameterInfo& params = constructor.params_in();

	for (size_t i = 0; i < constructor.num_parameter(); i++) {

		// create instance
		jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
		jobject obj = env->NewObject(cls, methodID);

		//assign values

		jmethodID setType = env->GetMethodID(cls, "setType", "(I)V");
		jmethodID setID = env->GetMethodID(cls, "setId", "(I)V");
		jmethodID setClassName = env->GetMethodID(cls, "setClassName",
				"(Ljava/lang/String;)V");
		jmethodID setClassNames = env->GetMethodID(cls, "setClassNames",
				"([Ljava/lang/String;)V");
		//	jmethodID setHelp = env->GetMethodID(cls,
		//		"setHelp", "(Ljava/lang/String;)V");
		//	jmethodID setToolTip = env->GetMethodID(cls,
		//		"setTooltip", "(Ljava/lang/String;)V");
		jmethodID setParamInfo = env->GetMethodID(cls, "setParamInfo",
				"([Ljava/lang/String;)V");

		using namespace ug::bridge;

		// TODO Unfortunately we don't know a better way to convert an enumeration
		// from C++ to Java. Currently we just use integers :(
		int type = paramType2Int(params, i);

		env->CallVoidMethod(obj, setType, type);
		env->CallVoidMethod(obj, setID, i);

		env->CallVoidMethod(obj, setClassName,
				stringC2J(env, params.class_name(i)));

		env->CallVoidMethod(obj, setClassNames,
				stringArrayC2J(env,
						getBaseClassNames(params.class_name_node(i))));

		env->CallVoidMethod(obj, setParamInfo,
				stringArrayC2J(env, constructor.parameter_info_vec(i)));

		// set array element
		env->SetObjectArrayElement(result, i, obj);
	}

	return result;
}

jobject retVal2NativeParam(JNIEnv *env,
		const ug::bridge::ExportedFunctionBase& func) {

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeParamInfo");

	unsigned int i = 0; // C/C++/Java only allow one return value

	const ug::bridge::ParameterInfo& params = func.params_out();

	// create instance
	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	jobject obj = env->NewObject(cls, methodID);

	// assign values

	jmethodID setType = env->GetMethodID(cls, "setType", "(I)V");
	jmethodID setID = env->GetMethodID(cls, "setId", "(I)V");
	jmethodID setClassName = env->GetMethodID(cls, "setClassName",
			"(Ljava/lang/String;)V");
	jmethodID setClassNames = env->GetMethodID(cls, "setClassNames",
			"([Ljava/lang/String;)V");
	//	jmethodID setHelp = env->GetMethodID(cls,
	//	"setHelp", "(Ljava/lang/String;)V");
	//	jmethodID setToolTip = env->GetMethodID(cls,
	//	"setTooltip", "(Ljava/lang/String;)V");
	jmethodID setParamInfo = env->GetMethodID(cls, "setParamInfo",
			"([Ljava/lang/String;)V");

	bool returnsVoid = func.params_out().size() == 0;

	// TODO Unfortunately we don't know a better way to convert an enumeration
	// from C++ to Java. Currently we just use integers :(
	int type = -1;

	if (returnsVoid) {
		type = -1; // void
	} else {
		type = paramType2Int(params, i);

		env->CallVoidMethod(obj, setClassName,
				stringC2J(env, params.class_name(i)));

		env->CallVoidMethod(obj, setClassNames,
				stringArrayC2J(env,
						getBaseClassNames(params.class_name_node(i))));
	}

	env->CallVoidMethod(obj, setType, type);
	env->CallVoidMethod(obj, setID, i);

	if (!returnsVoid) {
		env->CallVoidMethod(obj, setParamInfo,
				stringArrayC2J(env, func.return_info_vec()));
	}

	return obj;
}

jobject constructor2NativeConstructor(JNIEnv *env,
		const ug::bridge::ExportedConstructor* constructor) {
	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeConstructorInfo");

	// create instance

	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	jobject obj = env->NewObject(cls, methodID);

	//assign values

	jmethodID setHelp = env->GetMethodID(cls, "setHelp",
			"(Ljava/lang/String;)V");
	jmethodID setToolTip = env->GetMethodID(cls, "setToolTip",
			"(Ljava/lang/String;)V");
	jmethodID setOptions = env->GetMethodID(cls, "setOptions",
			"(Ljava/lang/String;)V");
	jmethodID setParameters = env->GetMethodID(cls, "setParameters",
			"([Ledu/gcsc/vrl/ug/NativeParamInfo;)V");

	using namespace ug::bridge;
	env->CallVoidMethod(obj, setHelp,
			stringC2J(env, constructor->help().c_str()));
	env->CallVoidMethod(obj, setToolTip,
			stringC2J(env, constructor->tooltip().c_str()));
	env->CallVoidMethod(obj, setOptions,
			stringC2J(env, constructor->options().c_str()));

	env->CallVoidMethod(obj, setParameters,
			params2NativeParams(env, *constructor));

	return obj;
}

jobjectArray constructors2NativeConstructors(JNIEnv *env,
		const ug::bridge::IExportedClass& eCls) {
	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeConstructorInfo");

	size_t numConstructors = eCls.num_constructors();

	jobjectArray result = env->NewObjectArray(numConstructors, cls, 0);

	for (size_t i = 0; i < numConstructors; i++) {

		jobject constructor = constructor2NativeConstructor(env,
				&eCls.get_constructor(i));

		env->SetObjectArrayElement(result, i, constructor);
	}

	return result;
}

jobject method2NativeMethod(JNIEnv *env,
		const ug::bridge::ExportedMethod* method) {
	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeMethodInfo");

	// create instance

	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	jobject obj = env->NewObject(cls, methodID);

	//assign values

	jmethodID setName = env->GetMethodID(cls, "setName",
			"(Ljava/lang/String;)V");
	jmethodID setHelp = env->GetMethodID(cls, "setHelp",
			"(Ljava/lang/String;)V");
	jmethodID setToolTip = env->GetMethodID(cls, "setToolTip",
			"(Ljava/lang/String;)V");
	jmethodID setOptions = env->GetMethodID(cls, "setOptions",
			"(Ljava/lang/String;)V");
	jmethodID setRetValue = env->GetMethodID(cls, "setReturnValue",
			"(Ledu/gcsc/vrl/ug/NativeParamInfo;)V");
	jmethodID setParameters = env->GetMethodID(cls, "setParameters",
			"([Ledu/gcsc/vrl/ug/NativeParamInfo;)V");

	using namespace ug::bridge;
	std::string name = method->name(); // TODO pre-rpocessing necessary
	env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
	env->CallVoidMethod(obj, setHelp, stringC2J(env, method->help().c_str()));
	env->CallVoidMethod(obj, setToolTip, stringC2J(env, name.c_str()));
	env->CallVoidMethod(obj, setOptions,
			stringC2J(env, method->options().c_str()));
	env->CallVoidMethod(obj, setRetValue, retVal2NativeParam(env, *method));

	env->CallVoidMethod(obj, setParameters, params2NativeParams(env, *method));

	return obj;
}

jobjectArray methods2NativeGroups(JNIEnv *env,
		const ug::bridge::IExportedClass& eCls, bool constMethods) {
	jclass groupCls = env->FindClass("edu/gcsc/vrl/ug/NativeMethodGroupInfo");
	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeMethodInfo");

	size_t numMethodGroups = 0;

	if (constMethods) {
		numMethodGroups = eCls.num_const_methods();
	} else {
		numMethodGroups = eCls.num_methods();
	}

	jobjectArray result = env->NewObjectArray(numMethodGroups, groupCls, 0);

	unsigned int numberOfMethodsInGroup = 1;

	for (size_t i = 0; i < numMethodGroups; i++) {

		if (constMethods) {
			numberOfMethodsInGroup = eCls.num_const_overloads(i);
		} else {
			numberOfMethodsInGroup = eCls.num_overloads(i);
		}

		jobjectArray methodArray = env->NewObjectArray(numberOfMethodsInGroup,
				cls, 0);

		jmethodID methodID = env->GetMethodID(groupCls, "<init>", "()V");
		jobject groupObj = env->NewObject(groupCls, methodID);
		jmethodID setOverloads = env->GetMethodID(groupCls, "setOverloads",
				"([Ledu/gcsc/vrl/ug/NativeMethodInfo;)V");

		for (size_t j = 0; j < numberOfMethodsInGroup; j++) {

			//--------- METHOD GROUP ---------

			const ug::bridge::ExportedMethod* method;

			if (constMethods) {
				method = &eCls.get_const_overload(i, j);
			} else {
				method = &eCls.get_overload(i, j);
			}

			// set array element, we currently have only one method per group
			env->SetObjectArrayElement(methodArray, j,
					method2NativeMethod(env, method));

		}

		env->CallVoidMethod(groupObj, setOverloads, methodArray);

		env->SetObjectArrayElement(result, i, groupObj);
	}

	return result;
}

jobject function2NativeFunction(JNIEnv *env,
		const ug::bridge::ExportedFunction& func) {
	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeFunctionInfo");

	// create instance

	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	jobject obj = env->NewObject(cls, methodID);

	//assign values

	jmethodID setName = env->GetMethodID(cls, "setName",
			"(Ljava/lang/String;)V");
	jmethodID setConst = env->GetMethodID(cls, "setConst", "(Z)V");
	jmethodID setCategory = env->GetMethodID(cls, "setCategory",
			"(Ljava/lang/String;)V");
	jmethodID setHelp = env->GetMethodID(cls, "setHelp",
			"(Ljava/lang/String;)V");
	jmethodID setToolTip = env->GetMethodID(cls, "setToolTip",
			"(Ljava/lang/String;)V");
	jmethodID setOptions = env->GetMethodID(cls, "setOptions",
			"(Ljava/lang/String;)V");
	jmethodID setRetValue = env->GetMethodID(cls, "setReturnValue",
			"(Ledu/gcsc/vrl/ug/NativeParamInfo;)V");
	jmethodID setParameters = env->GetMethodID(cls, "setParameters",
			"([Ledu/gcsc/vrl/ug/NativeParamInfo;)V");

	using namespace ug::bridge;
	std::string name = func.name(); // TODO pre-rpocessing necessary
	env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
	env->CallVoidMethod(obj, setConst, boolC2J(false));
	env->CallVoidMethod(obj, setCategory, stringC2J(env, func.group().c_str()));
	env->CallVoidMethod(obj, setHelp, stringC2J(env, func.help().c_str()));
	env->CallVoidMethod(obj, setToolTip, stringC2J(env, name.c_str()));
	env->CallVoidMethod(obj, setOptions,
			stringC2J(env, func.options().c_str()));
	env->CallVoidMethod(obj, setRetValue, retVal2NativeParam(env, func));

	env->CallVoidMethod(obj, setParameters, params2NativeParams(env, func));

	return obj;
}

jobjectArray functions2NativeGroups(JNIEnv * env, ug::bridge::Registry * reg) {
	jclass groupArrayCls = env->FindClass(
			"edu/gcsc/vrl/ug/NativeFunctionGroupInfo");

	unsigned int numFunctions = reg->num_functions();

	// create array of functions
	jobjectArray result = env->NewObjectArray(numFunctions, groupArrayCls, 0);

	for (size_t i = 0; i < numFunctions; i++) {

		const ug::bridge::ExportedFunctionGroup& group =
				reg->get_function_group(i);

		size_t numOverloads = group.num_overloads();

		jclass functionsCls = env->FindClass(
				"edu/gcsc/vrl/ug/NativeFunctionInfo");

		jobjectArray functions = env->NewObjectArray(numOverloads, functionsCls,
				0);

		for (size_t j = 0; j < numOverloads; j++) {
			const ug::bridge::ExportedFunction& func = *group.get_overload(j);

			env->SetObjectArrayElement(functions, j,
					function2NativeFunction(env, func));

		} // end for j

		// create function group instance
		jmethodID methodID = env->GetMethodID(groupArrayCls, "<init>", "()V");
		jobject obj = env->NewObject(groupArrayCls, methodID);

		jmethodID setOverloads = env->GetMethodID(groupArrayCls, "setOverloads",
				"([Ledu/gcsc/vrl/ug/NativeFunctionInfo;)V");

		env->CallVoidMethod(obj, setOverloads, functions);

		env->SetObjectArrayElement(result, i, obj);
	} // end for i

	return result;
}

jobjectArray classes2NativeClasses(JNIEnv *env,
		const ug::bridge::Registry* reg) {

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeClassInfo");

	jobjectArray result = env->NewObjectArray(reg->num_classes(), cls, 0);

	for (size_t i = 0; i < reg->num_classes(); i++) {

		const ug::bridge::IExportedClass& eCls = reg->get_class(i);

		// create instance

		jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
		jobject obj = env->NewObject(cls, methodID);

		// assign values

		jmethodID setName = env->GetMethodID(cls, "setName",
				"(Ljava/lang/String;)V");
		jmethodID setCategory = env->GetMethodID(cls, "setCategory",
				"(Ljava/lang/String;)V");
		jmethodID setClassNames = env->GetMethodID(cls, "setClassNames",
				"([Ljava/lang/String;)V");
		jmethodID setInstantiable = env->GetMethodID(cls, "setInstantiable",
				"(Z)V");
		jmethodID setConstructors = env->GetMethodID(cls, "setConstructors",
				"([Ledu/gcsc/vrl/ug/NativeConstructorInfo;)V");
		jmethodID setMethods = env->GetMethodID(cls, "setMethods",
				"([Ledu/gcsc/vrl/ug/NativeMethodGroupInfo;)V");
		jmethodID setConstMethods = env->GetMethodID(cls, "setConstMethods",
				"([Ledu/gcsc/vrl/ug/NativeMethodGroupInfo;)V");

		std::string name = eCls.name(); // TODO pre-rpocessing necessary

		// these lines check for empty class names. we really want this exit()
		// command as empty names will mess up everything.
		std::vector<std::string> baseClasses;

		for (size_t j = 0; j < eCls.class_names()->size(); j++) {

			if (eCls.class_names()->at(j) == NULL) {
				std::cerr << name << ", baseCls(" << j << ")==NULL"
						<< std::endl;
				exit(1);
			}

			if (strlen(eCls.class_names()->at(j)) == 0) {
				std::cerr << name << ", baseCls(" << j << ")==empty"
						<< std::endl;
				exit(1);
			}

			baseClasses.push_back(std::string(eCls.class_names()->at(j)));
		}

		env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
		env->CallVoidMethod(obj, setCategory,
				stringC2J(env, eCls.group().c_str()));
		env->CallVoidMethod(obj, setClassNames,
				stringArrayC2J(env, eCls.class_names()));
		//		env->CallVoidMethod(obj, setClassNames,
		//				stringArrayC2J(env, baseClasses));
		env->CallVoidMethod(obj, setInstantiable,
				boolC2J(eCls.is_instantiable()));
		env->CallVoidMethod(obj, setConstructors,
				constructors2NativeConstructors(env, eCls));
		env->CallVoidMethod(obj, setMethods,
				methods2NativeGroups(env, eCls, false));
		env->CallVoidMethod(obj, setConstMethods,
				methods2NativeGroups(env, eCls, true));

		// set array element
		env->SetObjectArrayElement(result, i, obj);
	}

	//	env->ExceptionCheck();
	//	env->ExceptionDescribe();

	return result;
}

jobjectArray classGroups2NativeClassGroups(JNIEnv *env,
		const ug::bridge::Registry* reg) {

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeClassGroupInfo");

	jobjectArray result = env->NewObjectArray(reg->num_class_groups(), cls, 0);

	for (size_t i = 0; i < reg->num_class_groups(); i++) {

		const ug::bridge::ClassGroupDesc* clsGrp = reg->get_class_group(i);

		// we ignore empty groups. actually, empty groups shouldn't exist.
		if (clsGrp->empty()) {
			continue;
		}

		// create string list containing class names of group i
		std::vector<std::string> class_names;

		for (size_t j = 0; j < clsGrp->num_classes(); j++) {
			class_names.push_back(std::string(clsGrp->get_class(j)->name()));
		}

		// initialize default class name with empty string. when trying to
		// instanciate from java we have to check that!
		std::string defaultClassName = "";

		if (clsGrp->get_default_class() != NULL) {
			defaultClassName = clsGrp->get_default_class()->name();
		}

		// create instance
		jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
		jobject obj = env->NewObject(cls, methodID);

		// assign values
		jmethodID setName = env->GetMethodID(cls, "setName",
				"(Ljava/lang/String;)V");
		jmethodID setClasses = env->GetMethodID(cls, "setClasses",
				"([Ljava/lang/String;)V");
		jmethodID setDefaultClass = env->GetMethodID(cls, "setDefaultClass",
				"(Ljava/lang/String;)V");

		// calls the java methods
		std::string name = clsGrp->name();
		env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
		env->CallVoidMethod(obj, setClasses, stringArrayC2J(env, class_names));
		env->CallVoidMethod(obj, setDefaultClass,
				stringC2J(env, defaultClassName.c_str()));

		// set array element
		env->SetObjectArrayElement(result, i, obj);
	}

	//		env->ExceptionCheck();
	//		env->ExceptionDescribe();

	return result;
}

jobject registry2NativeAPI(JNIEnv * env, ug::bridge::Registry * reg) {

	jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeAPIInfo");

	// create instance

	jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
	jobject obj = env->NewObject(cls, methodID);

	//assign values

	jmethodID setClassGroups = env->GetMethodID(cls, "setClassGroups",
			"([Ledu/gcsc/vrl/ug/NativeClassGroupInfo;)V");
	jmethodID setClasses = env->GetMethodID(cls, "setClasses",
			"([Ledu/gcsc/vrl/ug/NativeClassInfo;)V");
	jmethodID setFunctions = env->GetMethodID(cls, "setFunctions",
			"([Ledu/gcsc/vrl/ug/NativeFunctionGroupInfo;)V");

	env->CallVoidMethod(obj, setClassGroups,
			classGroups2NativeClassGroups(env, reg));
	env->CallVoidMethod(obj, setClasses, classes2NativeClasses(env, reg));
	env->CallVoidMethod(obj, setFunctions, functions2NativeGroups(env, reg));

	return obj;
}

} // end vrl::
} // end ug::
