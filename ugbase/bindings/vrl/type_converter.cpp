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
	std::string result(tmpStr);
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

bool isConstJPtr(JNIEnv* env, jobject obj)
{
	jclass classMethodAccess = env->FindClass("edu/gcsc/vrl/ug/Pointer");

	jmethodID classNameMethodID = env->GetMethodID(classMethodAccess,
			"isConst", "()Z");

	return (bool) env->CallBooleanMethod(obj, classNameMethodID);
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

	jclass argClass = env->GetObjectClass(obj);

	jmethodID getSmartPointer = env->GetMethodID(argClass, "getSmartPointer", "()[B");
	jbyteArray mem = (jbyteArray) env->CallObjectMethod(obj, getSmartPointer);
	jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

	// temporarily it is safe to use the memory provided by the Java byte array
	SmartPtr<void>* smartPtr = reinterpret_cast<SmartPtr<void>*>((void*) memPtr);

	// but we have to use a new instance for use outside of this function
	SmartPtr<void> result(*smartPtr);

	env->ReleaseByteArrayElements(mem, memPtr, 0);

	return result;
}

ConstSmartPtr<void> jObject2ConstSmartPointer(JNIEnv *env, jobject obj) {

	ConstSmartPtr<void> result;

	// We allow conversion from SmartPtr to ConstSmartPtr. As these classes
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

jobject constPointer2JObject(JNIEnv *env, const void* value) {

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

void throwUgErrorAsJavaException(JNIEnv *env, std::string error) {
	jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");

	std::stringstream ss;

	ss << "<!--NumMsg:" << 1 << ":-->";

	ss << "<div><pre>\n";

	ss << "<!--StartMsg:" << 0 << ":-->";
	ss << "\n" << error << "\n";
	ss << "<!--EndMsg:" << 0 << ":-->";

	if (error.size() != 0 && error[error.size() - 1] != '\n')
		ss << "\n";

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

	ug::Variant::Type tmpResultType = ug::Variant::VT_INVALID;

	std::string className = getClassName(env, obj);

	bool isArray = isJObjectAnArray(env, obj);
	typeAndArray.isArray = isArray;

	jobject sub_obj = NULL;

	if (isArray) {
		// if array is not empty, get first elem to find out type
		// this is necessary as the elements of this array could all be of
		// a inherited type, notably SmartPointer instead of Pointer!
		jobjectArray arr = (jobjectArray) obj;

		if (env->GetArrayLength(arr))
		{
			jobject sub_obj = env->GetObjectArrayElement(arr, 0);

			className = getClassName(env, sub_obj);
		}
		// if array empty, take the type the array is built on as type
		// FIXME: This will definitely cause problems if the expected UG type
		// is std::vector<SmartPtr<T> > and the passed VRL array is empty;
		// in this case, we will identify our non-existing object as typeof
		// "VT_POINTER" and this cannot be cast to VT_SMART_POINTER.
		else
		{// cut off the "object array" identifier: "[L"
			className = className.substr(2);

			// and cut off the tail if there is a semicolon ";"
			if (className.at(className.size()-1) == ';')
				className.erase(className.end()-1);
		}
	}


	if (className.compare("java.lang.Boolean") == 0) {
		tmpResultType = ug::Variant::VT_BOOL;
	} else if (className.compare("java.lang.Integer") == 0) {
		tmpResultType = ug::Variant::VT_INT;
	} else if (className.compare("java.lang.Double") == 0) {
		tmpResultType = ug::Variant::VT_DOUBLE;
	} else if (className.compare("java.lang.String") == 0) {
		tmpResultType = ug::Variant::VT_STDSTRING;
	} else if (className.compare("edu.gcsc.vrl.ug.Pointer") == 0)
	{
		if (isArray)
		{
			// if first object in array exists, determine its constness
			if (sub_obj)
			{
				if (isConstJPtr(env, sub_obj))
					tmpResultType = ug::Variant::VT_CONST_POINTER;
				else
					tmpResultType = ug::Variant::VT_POINTER;
			}
			// else (empty array) always take non-const,
			// since it can be cast to const
			else
			{
				tmpResultType = ug::Variant::VT_POINTER;
			}
		}
		else
		{
			if (isConstJPtr(env, obj))
				tmpResultType = ug::Variant::VT_CONST_POINTER;
			else
				tmpResultType = ug::Variant::VT_POINTER;
		}
	}
	else if (className.compare("edu.gcsc.vrl.ug.SmartPointer") == 0)
	{
		if (isArray)
		{
			// if first object in array exists, determine its constness
			if (sub_obj)
			{
				if (isConstJPtr(env, sub_obj))
					tmpResultType = ug::Variant::VT_CONST_SMART_POINTER;
				else
					tmpResultType = ug::Variant::VT_SMART_POINTER;
			}
			// else (empty array) always take non-const,
			// since it can be cast to const
			else
			{
				tmpResultType = ug::Variant::VT_SMART_POINTER;
			}
		}
		else
		{
			if (isConstJPtr(env, obj))
				tmpResultType = ug::Variant::VT_CONST_SMART_POINTER;
			else
				tmpResultType = ug::Variant::VT_SMART_POINTER;
		}
	}
	else
	{
		UG_LOG("\nparamClass2ParamType: Encountered invalid type " << className << ".\n");
	}
	// What about const pointers?
	// Answer: compareParamTypes() allows non-const* to const* conversion
	// That is why we do not check that. We also use the same class for
	// const and non const pointer. Const checking is done via readOnly
	// bit of the pointer instance (Java wrapper).

	// No, as compareParamTypes() only uses the ug::Variant::Type, but not the
	// wrappers' isConst() method, we better distinguish here!

	typeAndArray.type = tmpResultType;

//	if (DEBUG) {
//		std::cout << "ug::Variant::VT_INVALID = " << ug::Variant::VT_INVALID
//				<< std::endl;
//		std::cout << "paramClass2ParamType().tmpResultType = " << tmpResultType
//				<< std::endl;
//		std::cout << "VT_BOOL = " << ug::Variant::VT_BOOL << std::endl;
//		std::cout << "VT_INT = " << ug::Variant::VT_INT << std::endl;
//		std::cout << "VT_DOUBLE = " << ug::Variant::VT_DOUBLE << std::endl;
//		std::cout << "VT_STDSTRING = " << ug::Variant::VT_STDSTRING
//				<< std::endl;
//		std::cout << "VT_POINTER = " << ug::Variant::VT_POINTER << std::endl;
//		std::cout << "VT_SMART_POINTER = " << ug::Variant::VT_SMART_POINTER
//				<< std::endl;
//	}

	return typeAndArray;
}

bool compareParamTypes(JNIEnv *env, jobjectArray params,
		ug::bridge::Registry *reg,
		ug::bridge::ParameterInfo const& paramStack,
        bool allowSmartToRawPtrConversion) {

	// compare array lengths
	jsize len = env->GetArrayLength(params);

	if (len != paramStack.size()) {
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

			ug::vrl::throwUgErrorAsJavaException(env, ss.str());
		}

		TypeAndArray paramType = paramClass2ParamType(env, param);

		// check that is_vector == isArray
		if (paramType.isArray != paramStack.is_vector(i))
			return false;

		// allow non-const * to const *
		if (paramType.type == ug::Variant::VT_POINTER
				&& paramStack.type(i) == ug::Variant::VT_CONST_POINTER) {
			paramType.type = ug::Variant::VT_CONST_POINTER;
		}

		// allow non-const * to const *
		if (paramType.type == ug::Variant::VT_SMART_POINTER
				&& paramStack.type(i) == ug::Variant::VT_CONST_SMART_POINTER) {
			paramType.type = ug::Variant::VT_CONST_SMART_POINTER;
		}

		// allow std::string to c string
		if (paramType.type == ug::Variant::VT_STDSTRING
				&& paramStack.type(i) == ug::Variant::VT_CSTRING) {
			paramType.type = ug::Variant::VT_CSTRING;
		}

        if (allowSmartToRawPtrConversion) {

			// allow non-const-smart* to non-const*
			if (paramType.type == ug::Variant::VT_SMART_POINTER
					&& paramStack.type(i) == ug::Variant::VT_POINTER) {
				paramType.type = ug::Variant::VT_POINTER;
			}

			// allow non-const-smart* to const*
			if (paramType.type == ug::Variant::VT_SMART_POINTER
					&& paramStack.type(i) == ug::Variant::VT_CONST_POINTER) {
				paramType.type = ug::Variant::VT_CONST_POINTER;
			}

			// allow const smart* to const*
			if (paramType.type == ug::Variant::VT_CONST_SMART_POINTER
					&& paramStack.type(i) == ug::Variant::VT_CONST_POINTER) {

				paramType.type = ug::Variant::VT_CONST_POINTER;
			}
        }

		// allow integer to size_t
		if (paramType.type == ug::Variant::VT_INT
				&& paramStack.type(i) == ug::Variant::VT_SIZE_T) {
			paramType.type = ug::Variant::VT_SIZE_T;
		}

		if (paramType.type != (uint) paramStack.type(i)) {

			//#ifdef UG_DEBUG
			//			UG_LOG("requested by method:\n")
			//			printParamType(paramStack.type(i), i);
			//			UG_LOG("given as parameter:\n")
			//			printParamType(paramType, i);
			//#endif
			return false;
		}

		// check if param is assignable
		const ug::bridge::ClassNameNode* classNameNode =
				ug::vrl::invocation::getClassNodePtrByName(reg,
						getParamClassName(env, param));

		if (classNameNode != NULL) {

			if (paramStack.class_name(i) != NULL) {

				if (!ug::bridge::ClassNameTreeContains(*classNameNode,
						paramStack.class_name(i))) {

					return false;
				}
			}
		}
	}

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


void jobjectArray2ParamStack
(
	JNIEnv *env,
	ug::bridge::Registry* reg,
	ug::bridge::ParameterStack& paramsOut,
	const ug::bridge::ParameterInfo& paramsTemplate,
	jobjectArray const& array
)
{
	using namespace ug::bridge;

	//	iterate through the parameter list and copy the value in the
	//  associated stack entry.
	for (size_t i = 0; i < (size_t) paramsTemplate.size(); ++i) {

		ug::Variant::Type template_value_type = paramsTemplate.type(i);

		jobject value = env->GetObjectArrayElement(array, i);

		// we do not allow null values
		if (value == NULL)
		{
			std::stringstream ss;
			ss << "Value " << i << " == NULL!";

			jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
			env->ThrowNew(Exception, ss.str().c_str());
		}

		TypeAndArray jv_taa = paramClass2ParamType(env, value);
		ug::Variant::Type java_value_type = jv_taa.type;

		// if requested UG type is not a vector
		if (!paramsTemplate.is_vector(i))
		{
			switch (template_value_type)
			{
				case ug::Variant::VT_BOOL:
				{
					paramsOut.push<bool>(jObject2Boolean(env, value));
					break;
				}
				case ug::Variant::VT_INT:
				{
					paramsOut.push<int>(jObject2Int(env, value));
					break;
				}
				case ug::Variant::VT_SIZE_T:
				{
					paramsOut.push<size_t>((size_t) jObject2Int(env, value));
					break;
				}
				case ug::Variant::VT_FLOAT:
				case ug::Variant::VT_DOUBLE:
				{
					paramsOut.push<number>(jObject2Double(env, value));
					break;
				}
				case ug::Variant::VT_CSTRING:	// std::string can be converted to cstring by ug registry
				case ug::Variant::VT_STDSTRING:
				{
					paramsOut.push<std::string>(jObject2String(env, value));
					break;
				}
				case ug::Variant::VT_POINTER:
				{
					const ug::bridge::ClassNameNode* node =
						ug::vrl::invocation::getClassNodePtrByName(reg,
							jPointerGetName(env, value));
					paramsOut.push(jObject2Pointer(env, value), node);
					break;
				}
				case ug::Variant::VT_CONST_POINTER:
				{
					const ug::bridge::ClassNameNode* node =
						ug::vrl::invocation::getClassNodePtrByName(reg,
							jPointerGetName(env, value));

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
					break;
				}
				case ug::Variant::VT_SMART_POINTER:
				{
					const ug::bridge::ClassNameNode* node =
						ug::vrl::invocation::getClassNodePtrByName(reg,
							jPointerGetName(env, value));

					paramsOut.push(jObject2SmartPointer(env, value), node);
					break;
				}
				case ug::Variant::VT_CONST_SMART_POINTER:
				{
					const ug::bridge::ClassNameNode* node =
						ug::vrl::invocation::getClassNodePtrByName(reg,
							jPointerGetName(env, value));

					paramsOut.push(jObject2ConstSmartPointer(env, value), node);
					break;
				}
				case ug::Variant::VT_INVALID:
				{
					UG_THROW("Parameter in method template is invalid. This must not happen.");
				}
			} // end switch
		}

		// if requested UG type IS a vector
		else
		{
			switch (template_value_type)
			{
				case ug::Variant::VT_BOOL:
				{
					// cast to bool vector
					SmartPtr<std::vector<bool> > sp_bv(new std::vector<bool>());
					jObject2BoolVector(env, value, *(sp_bv.get()));

					// return a smart pointer
					paramsOut.push(sp_bv);

					break;
				}
				case ug::Variant::VT_INT:
				{
					// cast to number vector
					SmartPtr<std::vector<int> > sp_iv(new std::vector<int>());
					jObject2IntVector(env, value, *(sp_iv.get()));

					// return a smart pointer
					paramsOut.push(sp_iv);

					break;
				}
				case ug::Variant::VT_SIZE_T:
				{
					// cast to number vector
					SmartPtr<std::vector<size_t> > sp_stv(new std::vector<size_t>());
					jObject2SizetVector(env, value, *(sp_stv.get()));

					// return a smart pointer
					paramsOut.push(sp_stv);

					break;
				}
				case ug::Variant::VT_FLOAT:
				case ug::Variant::VT_DOUBLE:
				{
					// cast to number vector
					SmartPtr<std::vector<number> > sp_dv(new std::vector<number>());
					jObject2NumberVector(env, value, *(sp_dv.get()));

					// return a smart pointer
					paramsOut.push(sp_dv);

					break;
				}
				case ug::Variant::VT_CSTRING:	// std::string can be converted to cstring by ug registry
				case ug::Variant::VT_STDSTRING:
				{
					// cast to number vector
					SmartPtr<std::vector<std::string> > sp_sv(new std::vector<std::string>());
					jObject2stdStringVector(env, value, *(sp_sv.get()));

					// return a smart pointer
					paramsOut.push(sp_sv);

					break;
				}
				case ug::Variant::VT_POINTER:
				case ug::Variant::VT_CONST_POINTER:
				{
					// cast to void* vector
					SmartPtr<std::vector<std::pair<void*, const ClassNameNode*> > >
						sp_pv(new std::vector<std::pair<void*, const ClassNameNode*> >());

					jObject2PtrVector(env, value, java_value_type, reg, *(sp_pv.get()));

					// return a smart pointer
					paramsOut.push(sp_pv);

					break;
				}
				case ug::Variant::VT_SMART_POINTER:
				{
					// cast to void* vector
					SmartPtr<std::vector<std::pair<SmartPtr<void>, const ClassNameNode*> > >
						sp_pv(new std::vector<std::pair<SmartPtr<void>, const ClassNameNode*> >());

					jObject2SmartPtrVector(env, value, reg, *(sp_pv.get()));

					// return a smart pointer
					paramsOut.push(sp_pv);

					break;
				}
				case ug::Variant::VT_CONST_SMART_POINTER:
				{
					// cast to void* vector
					SmartPtr<std::vector<std::pair<ConstSmartPtr<void>, const ClassNameNode*> > >
						sp_pv(new std::vector<std::pair<ConstSmartPtr<void>, const ClassNameNode*> >());

					jObject2ConstSmartPtrVector(env, value, java_value_type, reg, *(sp_pv.get()));

					// return a smart pointer
					paramsOut.push(sp_pv);

					break;
				}
				case ug::Variant::VT_INVALID:
				{
					UG_THROW("Parameter in method template is invalid. This must not happen.");
				}
			} // end switch
		}
	} // end for
}

jobject param2JObject(JNIEnv *env, ug::bridge::ParameterStack& params, size_t index)
{
	using namespace ug::bridge;

	ug::Variant::Type type = params.type(index);
	bool isVector = params.is_vector(index);

	// distinguish between vector and "scalar" types
	switch (type)
	{
		case ug::Variant::VT_BOOL:
			if (!isVector) return boolean2JObject(env, params.to<bool>(index));
			return boolVector2JObject(env, params.to<std::vector<bool>& >(index));

		case ug::Variant::VT_INT:
			if (!isVector) return int2JObject(env, params.to<int>(index));
			return intVector2JObject(env, params.to<std::vector<int>& >(index));

		case ug::Variant::VT_SIZE_T:
			if (!isVector) return int2JObject(env, (int) params.to<size_t>(index));
			return sizetVector2JObject(env, params.to<std::vector<size_t>& >(index));

		case ug::Variant::VT_FLOAT:
		case ug::Variant::VT_DOUBLE:
			if (!isVector) return double2JObject(env, params.to<number>(index));
			return numberVector2JObject(env, params.to<std::vector<number>& >(index));

		case ug::Variant::VT_CSTRING:
			if (!isVector) return string2JObject(env, params.to<const char*>(index));
			return cStringVector2JObject(env, params.to<std::vector<const char*>& >(index));

		case ug::Variant::VT_STDSTRING:
			if (!isVector) return string2JObject(env, params.to<const std::string&>(index).c_str());
			return stdStringVector2JObject(env, params.to<std::vector<std::string>& >(index));

		case ug::Variant::VT_POINTER:
			if (!isVector) return pointer2JObject(env, params.to<void*>(index));
			return ptrVector2JObject(env, params.to<std::vector<void*>& >(index));

		case ug::Variant::VT_CONST_POINTER:
			if (!isVector) return constPointer2JObject(env, params.to<const void*>(index));
			return constPtrVector2JObject(env, params.to<std::vector<const void*>& >(index));

		case ug::Variant::VT_SMART_POINTER:
			if (!isVector) return smartPointer2JObject(env, params.to<SmartPtr<void> >(index));
			return smartPtrVector2JObject(env, params.to<std::vector<SmartPtr<void> >& >(index));

		case ug::Variant::VT_CONST_SMART_POINTER:
			if (!isVector) return constSmartPointer2JObject(env, params.to<ConstSmartPtr<void> >(index));
			return constSmartPtrVector2JObject(env, params.to<std::vector<ConstSmartPtr<void> >& >(index));

		case ug::Variant::VT_INVALID:
			UG_THROW("Parameter in method template is invalid. This must not happen.");
	}

	return jobject();
}


//	added by Christian Poliwoda
//	christian.poliwoda@gcsc.uni-frankfurt.de
//	y 13 m 05 d 28
bool isJObjectAnArray(JNIEnv *env, jobject value) {

	std::string valueClassName = getClassName(env, value);

// if the className starts with [ this is a hint
// that the value-object contains an array
	if (valueClassName.size() < 1) return false;
	if (valueClassName.at(0) == '[')
		return true;
	return false;
}


void jObject2BoolVector(JNIEnv *env, jobject object, std::vector<bool>& bv)
{
	// convert to array
	jobjectArray objectArray = (jobjectArray) object;
	jsize arrayLength = env->GetArrayLength(objectArray);

	// copy values to vector
	bv.resize((size_t) arrayLength);
	for (jsize i = 0; i < arrayLength; ++i)
		bv[i] = (bool) jObject2Boolean(env, env->GetObjectArrayElement(objectArray, i));

	return;
}

void jObject2IntVector(JNIEnv *env, jobject object, std::vector<int>& iv)
{
	// convert to array
	jobjectArray objectArray = (jobjectArray) object;
	jsize arrayLength = env->GetArrayLength(objectArray);

	// copy values to vector
	iv.resize((size_t) arrayLength);
	for (jsize i = 0; i < arrayLength; ++i)
		iv[i] = (int) jObject2Int(env, env->GetObjectArrayElement(objectArray, i));

	return;
}

void jObject2SizetVector(JNIEnv *env, jobject object, std::vector<size_t>& stv)
{
	// convert to array
	jobjectArray objectArray = (jobjectArray) object;
	jsize arrayLength = env->GetArrayLength(objectArray);

	// copy values to vector
	stv.resize((size_t) arrayLength);
	for (jsize i = 0; i < arrayLength; ++i)
		stv[i] = (size_t) jObject2Int(env, env->GetObjectArrayElement(objectArray, i));

	return;
}

void jObject2NumberVector(JNIEnv *env, jobject object, std::vector<number>& nv)
{
	// convert to array
	jobjectArray objectArray = (jobjectArray) object;
	jsize arrayLength = env->GetArrayLength(objectArray);

	// copy values to vector
	nv.resize((size_t) arrayLength);
	for (jsize i = 0; i < arrayLength; ++i)
		nv[i] = (number) jObject2Double(env, env->GetObjectArrayElement(objectArray, i));

	return;
}

void jObject2stdStringVector(JNIEnv *env, jobject object, std::vector<std::string>& sv)
{
	// convert to array
	jobjectArray objectArray = (jobjectArray) object;
	jsize arrayLength = env->GetArrayLength(objectArray);

	// copy values to vector
	sv.resize((size_t) arrayLength);
	for (jsize i = 0; i < arrayLength; ++i)
		sv[i] = jObject2String(env, env->GetObjectArrayElement(objectArray, i));

	return;
}

void jObject2PtrVector
(
	JNIEnv *env,
	jobject object,
	ug::Variant::Type jo_type,
	ug::bridge::Registry* reg,
	std::vector<std::pair<void*, const ug::bridge::ClassNameNode*> >& pv
)
{
	// convert to array
	jobjectArray objectArray = (jobjectArray) object;
	jsize arrayLength = env->GetArrayLength(objectArray);

	// copy values to vector
	pv.resize((size_t) arrayLength);
	for (jsize i = 0; i < arrayLength; ++i)
	{
		const ug::bridge::ClassNameNode* node =
			ug::vrl::invocation::getClassNodePtrByName(reg,
				jPointerGetName(env, env->GetObjectArrayElement(objectArray, i)));

		if (jo_type == ug::Variant::VT_SMART_POINTER)
		{
			pv[i] = std::make_pair<>
				((void*) jObject2SmartPointer(env, env->GetObjectArrayElement(objectArray, i)).get(), node);
		}
		else if (jo_type == ug::Variant::VT_CONST_SMART_POINTER)
		{
			pv[i] = std::make_pair<>
				((void*) jObject2ConstSmartPointer(env, env->GetObjectArrayElement(objectArray, i)).get(), node);
		}
		else	// VT_POINTER or VT_CONST_POINTER
		{
			pv[i] = std::make_pair<>
				(jObject2Pointer(env, env->GetObjectArrayElement(objectArray, i)), node);
		}
	}

	return;
}

void jObject2SmartPtrVector
(
	JNIEnv *env,
	jobject object,
	ug::bridge::Registry* reg,
	std::vector<std::pair<SmartPtr<void>, const ug::bridge::ClassNameNode*> >& pv
)
{
	// convert to array
	jobjectArray objectArray = (jobjectArray) object;
	jsize arrayLength = env->GetArrayLength(objectArray);

	// copy values to vector
	pv.resize((size_t) arrayLength);
	for (jsize i = 0; i < arrayLength; ++i)
	{
		const ug::bridge::ClassNameNode* node =
			ug::vrl::invocation::getClassNodePtrByName(reg,
				jPointerGetName(env, env->GetObjectArrayElement(objectArray, i)));

		pv[i] = std::make_pair<>
				(jObject2SmartPointer(env, env->GetObjectArrayElement(objectArray, i)), node);
	}

	return;
}

void jObject2ConstSmartPtrVector
(
	JNIEnv *env,
	jobject object,
	ug::Variant::Type jo_type,
	ug::bridge::Registry* reg,
	std::vector<std::pair<ConstSmartPtr<void>, const ug::bridge::ClassNameNode*> >& pv
)
{
	// convert to array
	jobjectArray objectArray = (jobjectArray) object;
	jsize arrayLength = env->GetArrayLength(objectArray);

	// copy values to vector
	pv.resize((size_t) arrayLength);
	for (jsize i = 0; i < arrayLength; ++i)
	{
		const ug::bridge::ClassNameNode* node =
			ug::vrl::invocation::getClassNodePtrByName(reg,
				jPointerGetName(env, env->GetObjectArrayElement(objectArray, i)));

		if (jo_type == ug::Variant::VT_SMART_POINTER)
		{
			pv[i] = std::make_pair<>
				(jObject2SmartPointer(env, env->GetObjectArrayElement(objectArray, i)), node);
		}
		else	// VT_CONST_SMART_POINTER
		{
			pv[i] = std::make_pair<>
				(jObject2ConstSmartPointer(env, env->GetObjectArrayElement(objectArray, i)), node);
		}
	}

	return;
}


jobject boolVector2JObject(JNIEnv* env, const std::vector<bool>& bv)
{
	size_t size = bv.size();
	jclass cls = env->FindClass("Ljava/lang/Boolean;");

	// allocate array of Booleans (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, boolean2JObject(env, (jboolean) bv[i]));

	return array;
}


jobject intVector2JObject(JNIEnv* env, const std::vector<int>& iv)
{
	size_t size = iv.size();
	jclass cls = env->FindClass("Ljava/lang/Integer;");

	// allocate array of Integers (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, int2JObject(env, (jint) iv[i]));

	return array;
}

jobject sizetVector2JObject(JNIEnv* env, const std::vector<size_t>& iv)
{
	size_t size = iv.size();
	jclass cls = env->FindClass("Ljava/lang/Integer;");

	// allocate array of Integers (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, int2JObject(env, (jint) iv[i]));

	return array;
}

jobject numberVector2JObject(JNIEnv* env, const std::vector<number>& nv)
{
	size_t size = nv.size();
	jclass cls = env->FindClass("Ljava/lang/Double;");

	// allocate array of Doubles (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, double2JObject(env, (jdouble) nv[i]));

	return array;
}

jobject cStringVector2JObject(JNIEnv* env, const std::vector<const char*>& sv)
{
	size_t size = sv.size();
	jclass cls = env->FindClass("Ljava/lang/String;");

	// allocate array of Strings (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, string2JObject(env, sv[i]));

	return array;
}

jobject stdStringVector2JObject(JNIEnv* env, const std::vector<std::string>& sv)
{
	size_t size = sv.size();
	jclass cls = env->FindClass("Ljava/lang/String;");

	// allocate array of Strings (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, string2JObject(env, sv[i].c_str()));

	return array;
}

jobject ptrVector2JObject(JNIEnv* env, const std::vector<void*>& pv)
{
	size_t size = pv.size();
	jclass cls = env->FindClass("Ledu/gcsc/vrl/ug/Pointer;");

	// allocate array of Pointers (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, pointer2JObject(env, pv[i]));

	return array;
}

jobject constPtrVector2JObject(JNIEnv* env, const std::vector<const void*>& pv)
{
	size_t size = pv.size();
	jclass cls = env->FindClass("Ledu/gcsc/vrl/ug/Pointer;");

	// allocate array of Pointers (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, constPointer2JObject(env, pv[i]));

	return array;
}

jobject smartPtrVector2JObject(JNIEnv* env, const std::vector<SmartPtr<void> >& pv)
{
	size_t size = pv.size();
	jclass cls = env->FindClass("Ledu/gcsc/vrl/ug/SmartPointer;");

	// allocate array of SmartPointers (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, smartPointer2JObject(env, pv[i]));

	return array;
}

jobject constSmartPtrVector2JObject(JNIEnv* env, const std::vector<ConstSmartPtr<void> >& pv)
{
	size_t size = pv.size();
	jclass cls = env->FindClass("Ledu/gcsc/vrl/ug/SmartPointer;");

	// allocate array of SmartPointers (constructing entries using NULL)
	jobjectArray array = env->NewObjectArray((jsize) size, cls, (jobject) NULL);

	for (size_t i = 0; i < size; ++i)
		env->SetObjectArrayElement(array, (jsize) i, constSmartPointer2JObject(env, pv[i]));

	return array;
}


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

		jmethodID setVectorFlagMethod = env->GetMethodID(cls, "setVectorFlag", "(Z)V");

		// check for emptiness
		bool pointerType = params.type(i) == ug::Variant::VT_CONST_POINTER
				|| params.type(i) == ug::Variant::VT_POINTER
				|| params.type(i) == ug::Variant::VT_SMART_POINTER
				|| params.type(i) == ug::Variant::VT_CONST_SMART_POINTER;

		if (pointerType && strlen(params.class_name(i)) == 0) {
			std::cerr << "ERROR: empty class name in "
                                  << func.name() << ", param(" << i << ")"
				  << std::endl;
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

		env->CallVoidMethod(obj, setVectorFlagMethod, (jboolean) params.is_vector(i));

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

		jmethodID setVectorFlagMethod = env->GetMethodID(cls, "setVectorFlag", "(Z)V");

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

		env->CallVoidMethod(obj, setVectorFlagMethod, (jboolean) params.is_vector(i));

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

	jmethodID setVectorFlagMethod = env->GetMethodID(cls, "setVectorFlag", "(Z)V");

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
		env->CallVoidMethod(obj, setVectorFlagMethod, (jboolean) params.is_vector(i));
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
	std::string name = method->name(); // TODO pre-processing necessary
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
	std::string name = func.name(); // TODO pre-processing necessary
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

		std::string name = eCls.name(); // TODO pre-processing necessary

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
		// instantiate from java we have to check that!
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

	if (!reg->check_consistency()) {
	ug::GetLogAssistant().flush_error_log();
	UG_LOG("trunk/ugbase/bindings/vrl/type_converter.cpp : registry2NativeAPI() calling check_consistency() \n");
	//return 1;
	}

	return obj;
}

} // end vrl::
} // end ug::
