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
			std::string result = (std::string)tmpStr;
			env->ReleaseStringUTFChars(s, tmpStr);
			return result;
		}

		jobjectArray stringArrayC2J(
				JNIEnv *env,
				const std::string* strings,
				const size_t array_length) {
			jclass stringClass = env->FindClass("java/lang/String");

			jobjectArray result =
					env->NewObjectArray(array_length, stringClass, 0);

			// convert array elements
			for (size_t i = 0; i < array_length; i++) {
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
			for (size_t i = 0; i < array_length; i++) {
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

		jobjectArray stringArrayC2J(
				JNIEnv *env, const std::vector<const char*>* strings) {

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

		jobjectArray stringArrayC2J(
				JNIEnv *env, const std::vector<const char*> strings) {

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

		std::vector<std::string> stringArrayJ2C(
				JNIEnv *env, jobjectArray const& array) {

			std::vector<std::string> result;

			size_t length = env->GetArrayLength(array);

			// convert each element of the java object array to a std string
			// and add it to the result vector
			for (size_t i = 0; i < length; i++) {
				result.push_back(stringJ2C(env,
						(jstring) env->GetObjectArrayElement(array, i)));
			}

			return result;
		}

		const std::vector<const ug::bridge::IExportedClass*> getParentClasses(
				ug::bridge::Registry* reg,
				const ug::bridge::IExportedClass* clazz) {
			std::vector<const ug::bridge::IExportedClass*> result;
			// search registered classes by name as specified in the
			// class_names vector and add them to the result vector
			for (size_t i = 0; i < clazz->class_names()->size(); i++) {
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
			jmethodID methodID = env->GetMethodID(argClass, "booleanValue", "()Z");
			return env->CallBooleanMethod(obj, methodID);
		}

		std::vector<const char*> getBaseClassNames(const ug::bridge::ClassNameNode* node) {
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
			jmethodID methodID = env->GetMethodID(argClass, "getClassName", "()Ljava/lang/String;");

			std::string name =
					stringJ2C(env, (jstring) env->CallObjectMethod(obj, methodID));

			return name;
		}

		SmartPtr<void> jObject2SmartPointer(JNIEnv *env, jobject obj) {

			jclass argClass = env->GetObjectClass(obj);

			jmethodID getSmartPointer =
					env->GetMethodID(argClass, "getSmartPointer", "()[B");
			jbyteArray mem = (jbyteArray) env->CallObjectMethod(obj, getSmartPointer);
			jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

			// temporarily it is save to use the memory provided by the Java
			// byte array
			SmartPtr<void>* smartPtr =
					reinterpret_cast<SmartPtr<void>*> ((void*) memPtr);

			// but we have to use a new instance for use outside of this function
			SmartPtr<void> result(*smartPtr);

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
				jmethodID getSmartPointer =
						env->GetMethodID(argClass, "getSmartPointer", "()[B");
				jbyteArray mem =
						(jbyteArray) env->CallObjectMethod(obj, getSmartPointer);
				jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

				ConstSmartPtr<void>* smartPtr =
						reinterpret_cast<ConstSmartPtr<void>*> ((void*) memPtr);
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

			jmethodID getSmartPointer =
					env->GetMethodID(argClass, "getSmartPointer", "()[B");
			jbyteArray mem = (jbyteArray) env->CallObjectMethod(obj, getSmartPointer);
			jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

			// It is important to call the invalidate method directly on the
			// smart-pointer instance that has been retrieved from the java object.
			// It is not valid to use copy-constructor etc. and to work on a stack
			// instance because invalidation of the stack instance won't allow the
			// destructor to decrease the reference count. This is equivalent to
			// not invoking the invalidate() method.
			(reinterpret_cast<SmartPtr<void>*> ((void*) memPtr))->invalidate();

			env->ReleaseByteArrayElements(mem, memPtr, 0);
		}

		void invalidateJConstSmartPointer(JNIEnv *env, jobject obj) {

			jclass argClass = env->GetObjectClass(obj);

			jmethodID getSmartPointer =
					env->GetMethodID(argClass, "getSmartPointer", "()[B");
			jbyteArray mem = (jbyteArray) env->CallObjectMethod(obj, getSmartPointer);
			jbyte* memPtr = env->GetByteArrayElements(mem, NULL);

			// It is important to call the invalidate method directly on the
			// smart-pointer instance that has been retrieved from the java object.
			// It is not valid to use copy-constructor etc. and to work on a stack
			// instance because invalidation of the stack instance won't allow the
			// destructor to decrease the reference count. This is equivalent to
			// not invoking the invalidate() method.
			(reinterpret_cast<ConstSmartPtr<void>*> ((void*) memPtr))->invalidate();

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

			jsize size = sizeof (SmartPtr<void>);
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
			jobject result = env->NewObject(cls, methodID, (jlong) value.get(),
					mem, boolC2J(false));

			return result;
		}

		jobject constSmartPointer2JObject(JNIEnv *env, ConstSmartPtr<void> value) {

			jsize size = sizeof (ConstSmartPtr<void>);
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
			jobject result = env->NewObject(cls, methodID, (jlong) value.get(),
					mem, boolC2J(true));

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

			ss << "<!--NumMsg:"<<error.num_msg()<<":-->";
			ss << "<div><pre>\n";
			for(size_t i = 0; i < error.num_msg(); i++) {

				ss << "<!--StartMsg:"<<i<<":-->";
				std::string msg = error.get_msg(i);
				ss << "\n" << msg << "\n";
				ss << "<!--EndMsg:"<<i<<":-->";

				if (msg.size() !=0 && msg[msg.size()-1] != '\n') {
					ss << "\n";
				}

				ss	<< ">> File:\t";
				ss << "<!--StartFile:"<<i<<":-->";
				ss	<< error.get_file(i);
				ss << "<!--EndFile:"<<i<<":-->";
				ss << "\n";

				ss	<< ">> Line:\t";
				ss << "<!--StartLine:"<<i<<":-->";
				ss  << error.get_line(i);
				ss << "<!--EndLine:"<<i<<":-->";
				ss  << "\n\n";

				if (i < error.num_msg()-1) {
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

		std::string getParamClassName(JNIEnv *env, jobject obj) {

			jclass cls = env->FindClass("edu/gcsc/vrl/ug/Pointer");

			if (env->ExceptionCheck()) {
				env->ExceptionDescribe();
			}

			jmethodID instanceMethodID =
					env->GetStaticMethodID(cls, "isInstance", "(Ljava/lang/Object;)Z");

			jboolean instanceof =
					env->CallStaticBooleanMethod(cls, instanceMethodID, obj);

			if (boolJ2C(instanceof)) {
				jmethodID classNameMethodID = env->GetMethodID(cls,
						"getClassName", "()Ljava/lang/String;");

				if (env->ExceptionCheck()) {
					env->ExceptionDescribe();
				}
				jobject resultObj = env->CallObjectMethod(obj,
						classNameMethodID);

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
				case PT_UNKNOWN:
				{
					UG_LOG("Param " << index << " = PT_UNKNOWN" << std::endl)
				}
					break;
				case PT_BOOL:
				{
					UG_LOG("Param " << index << " = PT_BOOL" << std::endl)
				}
					break;
				case PT_INTEGER:
				{
					UG_LOG("Param " << index << " = PT_INTEGER" << std::endl)
				}
					break;
				case PT_NUMBER:
				{
					UG_LOG("Param " << index << " = PT_NUMBER" << std::endl)
				}
					break;
				case PT_CSTRING:
				{
					UG_LOG("Param " << index << " = PT_CSTRING" << std::endl)
				}
					break;
				case PT_STD_STRING:
				{
					UG_LOG("Param " << index << " = PT_STD_STRING" << std::endl)
				}
					break;
				case PT_POINTER:
				{
					UG_LOG("Param " << index << " = PT_POINTER" << std::endl)
				}
					break;
				case PT_CONST_POINTER:
				{
					UG_LOG("Param " << index << " = PT_CONST_POINTER" << std::endl)
				}
					break;
				case PT_SMART_POINTER:
				{
					UG_LOG("Param " << index << " = PT_SMART_POINTER" << std::endl)
				}
					break;
				case PT_CONST_SMART_POINTER:
				{
					UG_LOG("Param " << index <<
							" = PT_CONST_SMART_POINTER" << std::endl)
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
			case PT_UNKNOWN:
			{
				result = "PT_UNKNOWN";
			}
				break;
			case PT_BOOL:
			{
				result = "PT_BOOL";
			}
				break;
			case PT_INTEGER:
			{
				result = "PT_INTEGER";
			}
				break;
			case PT_NUMBER:
			{
				result = "PT_NUMBER";
			}
				break;
			case PT_CSTRING:
			{
				result ="PT_CSTRING";
			}
				break;
			case PT_STD_STRING:
			{
				result = "PT_STD_STRING";
			}
				break;
			case PT_POINTER:
			{
				result = "PT_POINTER";
			}
				break;
			case PT_CONST_POINTER:
			{
				result ="PT_CONST_POINTER";
			}
				break;
			case PT_SMART_POINTER:
			{
				result = "PT_SMART_POINTER";
			}
				break;
			case PT_CONST_SMART_POINTER:
			{
				result = "PT_CONST_SMART_POINTER";
			}
				break;
			default:
				break;
			}

			return result;
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
				result = ug::bridge::PT_STD_STRING;
			} else if (className.compare("edu.gcsc.vrl.ug.Pointer") == 0) {
				result = ug::bridge::PT_POINTER;
			} else if (className.compare("edu.gcsc.vrl.ug.SmartPointer") == 0) {
				result = ug::bridge::PT_SMART_POINTER;
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
				ug::bridge::Registry *reg,
				ug::bridge::ParameterStack const& paramStack) {

			//#ifdef UG_DEBUG
			//	UG_LOG("\n -- BEGIN COMPARE --\n")
			//#endif

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

					jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
					env->ThrowNew(Exception, ss.str().c_str());
				}

				uint paramType = paramClass2ParamType(env, param);

				// allow non-const * to const *
				if (paramType == ug::bridge::PT_POINTER &&
						paramStack.get_type(i) == ug::bridge::PT_CONST_POINTER) {
					paramType = ug::bridge::PT_CONST_POINTER;
				}

				// allow non-const * to const *
				if (paramType == ug::bridge::PT_SMART_POINTER &&
						paramStack.get_type(i) == ug::bridge::PT_CONST_SMART_POINTER) {
					paramType = ug::bridge::PT_CONST_SMART_POINTER;
				}


				// allow std::string to c string
				if (paramType == ug::bridge::PT_STD_STRING &&
						paramStack.get_type(i) == ug::bridge::PT_CSTRING) {
					paramType = ug::bridge::PT_CSTRING;
				}


				// UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!) 
				// allow non-const-smart* to non const*
				if (paramType == ug::bridge::PT_SMART_POINTER &&
						paramStack.get_type(i) == ug::bridge::PT_POINTER) {
					paramType = ug::bridge::PT_POINTER;
				}

				// allow non-const-smart* to const*
				if (paramType == ug::bridge::PT_SMART_POINTER &&
						paramStack.get_type(i) == ug::bridge::PT_CONST_POINTER) {
					paramType = ug::bridge::PT_CONST_POINTER;
				}

				// allow const smart* to const*
				if (paramType == ug::bridge::PT_CONST_SMART_POINTER &&
						paramStack.get_type(i) == ug::bridge::PT_CONST_POINTER) {
					paramType = ug::bridge::PT_CONST_POINTER;
				}

				if (paramType != paramStack.get_type(i)) {
					//#ifdef UG_DEBUG
					//			UG_LOG("requested by method:\n")
					//			printParamType(paramStack.get_type(i), i);
					//			UG_LOG("given as parameter:\n")
					//			printParamType(paramType, i);
					//#endif
					return false;
				}

				// check if param is assignable
				const ug::bridge::ClassNameNode* classNameNode =
						ug::vrl::invocation::getClassNodePtrByName(
						reg, getParamClassName(env, param));

				if (classNameNode != NULL) {
					if (paramStack.class_name(i) != NULL) {
						if (!ug::bridge::ClassNameTreeContains(
								*classNameNode, paramStack.class_name(i))) {
							return false;
						}
					}
				}
			}

			//#ifdef UG_DEBUG
			//	UG_LOG(" -- ALL TRUE --\n" << std::endl)
			//#endif

			return true;
		}

		std::string getParamTypesAsString(JNIEnv *env, jobjectArray const& array) {

			std::stringstream ss;

			size_t length = (size_t)env->GetArrayLength(array);

			//iterate through the parameter array and add the type name to
			// the stream
			for (size_t i = 0; i < length; ++i) {
				jobject value = env->GetObjectArrayElement(array, i);

				uint java_value_type = paramClass2ParamType(env, value);

				if (i > 0) {
					ss << ", ";
				}
				bool classValue =
						java_value_type == ug::bridge::PT_CONST_POINTER ||
						java_value_type == ug::bridge::PT_POINTER ||
						java_value_type == ug::bridge::PT_SMART_POINTER ||
						java_value_type == ug::bridge::PT_CONST_SMART_POINTER;

				if (classValue) {
					ss << getParamClassName(env,value);
				} else {
					ss << getParamTypeAsString(java_value_type);
				}

			}

			return ss.str();
		}

		void jobjectArray2ParamStack(
				JNIEnv *env, ug::bridge::Registry* reg,
				ug::bridge::ParameterStack& paramsOut,
				const ug::bridge::ParameterStack& paramsTemplate,
				jobjectArray const& array) {
			using namespace ug::bridge;

			//	iterate through the parameter list and copy the value in the
			//  associated stack entry.
			for (size_t i = 0; i < (size_t) paramsTemplate.size(); ++i) {

				uint template_value_Type = paramsTemplate.get_type(i);

				jobject value = env->GetObjectArrayElement(array, i);

				// only used for
				// UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!) 
				uint java_value_type = paramClass2ParamType(env, value);

				// we don't allow null values
				if (value == NULL) {
					std::stringstream ss;
					ss << "Value " << i << " == NULL!";

					jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
					env->ThrowNew(Exception, ss.str().c_str());
				}

				switch (template_value_Type) {
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
					case PT_CSTRING:
					{
						paramsOut.push_std_string(jObject2String(env, value));
					}
						break;

					case PT_STD_STRING:
					{
						paramsOut.push_std_string(jObject2String(env, value));
					}
						break;

					case PT_POINTER:
					{
						const ug::bridge::ClassNameNode* node =
								ug::vrl::invocation::getClassNodePtrByName(reg,
								jPointerGetName(env, value));
						paramsOut.push_pointer(jObject2Pointer(env, value), node);
					}
						break;
					case PT_CONST_POINTER:
					{
						const ug::bridge::ClassNameNode* node =
								ug::vrl::invocation::getClassNodePtrByName(reg,
								jPointerGetName(env, value));

						// UGLY SMART-PTR to RAW-PTR CONVERSION (don't use this!) 
						if (java_value_type == PT_CONST_SMART_POINTER) {
							paramsOut.push_const_pointer(
									(void*) jObject2ConstSmartPointer(
									env, value).get(), node);
						} else if (java_value_type == PT_SMART_POINTER) {
							paramsOut.push_const_pointer(
									(void*) jObject2SmartPointer(
									env, value).get(), node);
						} else {
							paramsOut.push_const_pointer(
									jObject2Pointer(env, value), node);
						}

//						paramsOut.push_const_pointer(
//								jObject2Pointer(env, value), node);
					}
						break;
					case PT_SMART_POINTER:
					{
						const ug::bridge::ClassNameNode* node =
								ug::vrl::invocation::getClassNodePtrByName(reg,
								jPointerGetName(env, value));

						paramsOut.push_smart_pointer(
								jObject2SmartPointer(env, value), node);
						//				}
					}
						break;
					case PT_CONST_SMART_POINTER:
					{
						const ug::bridge::ClassNameNode* node =
								ug::vrl::invocation::getClassNodePtrByName(reg,
								jPointerGetName(env, value));

						paramsOut.push_const_smart_pointer(
								jObject2ConstSmartPointer(env, value), node);
					}
						break;
				}

			} // end for
		}

		jobject param2JObject(
				JNIEnv *env, ug::bridge::ParameterStack& params, size_t index) {
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
				case PT_CSTRING:
				{
					return string2JObject(env, params.to_cstring(index));
				}
					break;
				case PT_STD_STRING:
				{
					return string2JObject(env, params.to_std_string(index).c_str());
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
				case PT_SMART_POINTER:
				{
					return smartPointer2JObject(env, params.to_smart_pointer(index));
				}
					break;
				case PT_CONST_SMART_POINTER:
				{
					return constSmartPointer2JObject(
							env, params.to_const_smart_pointer(index));
				}
					break;
			}

			return jobject();
		}

		int paramType2Int(const ug::bridge::ParameterStack& params, size_t index) {
			using namespace ug::bridge;

			int type = params.get_type(index);

			return type;
		}

		jobjectArray params2NativeParams(JNIEnv *env,
				const ug::bridge::ExportedFunctionBase& func) {

			jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeParamInfo");

			jobjectArray result =
					env->NewObjectArray(func.num_parameter(), cls, 0);

			const ug::bridge::ParameterStack& params = func.params_in();

			for (size_t i = 0; i < func.num_parameter(); i++) {

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


				// check for emptyness
				bool pointerType = params.get_type(i) == ug::bridge::PT_CONST_POINTER ||
						params.get_type(i) == ug::bridge::PT_POINTER ||
						params.get_type(i) == ug::bridge::PT_SMART_POINTER ||
						params.get_type(i) == ug::bridge::PT_CONST_SMART_POINTER;

				if (pointerType && strlen(params.class_name(i)) == 0) {
					std::cerr << func.name() << ", param(" << i << ")==EMPTY" << std::endl;
					exit(1);
				}


				using namespace ug::bridge;

				// TODO Unfortunately we don't know a better way to convert an enumeration
				// from C++ to Java. Currently we just use integers :(
				int type = paramType2Int(params, i);

				env->CallVoidMethod(obj, setType, type);
				env->CallVoidMethod(obj, setID, i);

				env->CallVoidMethod(obj, setClassName, stringC2J(env, params.class_name(i)));

				env->CallVoidMethod(obj, setClassNames, stringArrayC2J(env,
						getBaseClassNames(params.class_name_node(i))));

				env->CallVoidMethod(obj, setParamInfo,
						stringArrayC2J(env, func.parameter_info_vec(i)));

				// set array element
				env->SetObjectArrayElement(result, i, obj);
			}

			return result;
		}

		jobjectArray params2NativeParams(JNIEnv *env,
				const ug::bridge::ExportedConstructor& constructor) {

			jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeParamInfo");

			jobjectArray result =
					env->NewObjectArray(constructor.num_parameter(), cls, 0);

			const ug::bridge::ParameterStack& params = constructor.params_in();

			for (size_t i = 0; i < constructor.num_parameter(); i++) {

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

				env->CallVoidMethod(obj, setType, type);
				env->CallVoidMethod(obj, setID, i);

				env->CallVoidMethod(obj, setClassName, stringC2J(env, params.class_name(i)));

				env->CallVoidMethod(obj, setClassNames, stringArrayC2J(env,
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

			const ug::bridge::ParameterStack& params = func.params_out();

			// create instance
			jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
			jobject obj = env->NewObject(cls, methodID);

			// assign values

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

			jmethodID setHelp = env->GetMethodID(cls,
					"setHelp", "(Ljava/lang/String;)V");
			jmethodID setToolTip = env->GetMethodID(cls,
					"setToolTip", "(Ljava/lang/String;)V");
			jmethodID setOptions = env->GetMethodID(cls,
					"setOptions", "(Ljava/lang/String;)V");
			jmethodID setParameters = env->GetMethodID(cls,
					"setParameters", "([Ledu/gcsc/vrl/ug/NativeParamInfo;)V");

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

			jobjectArray result =
					env->NewObjectArray(numConstructors, cls, 0);

			for (size_t i = 0; i < numConstructors; i++) {

				jobject constructor =
						constructor2NativeConstructor(env, &eCls.get_constructor(i));

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

			jmethodID setName = env->GetMethodID(cls,
					"setName", "(Ljava/lang/String;)V");
			jmethodID setHelp = env->GetMethodID(cls,
					"setHelp", "(Ljava/lang/String;)V");
			jmethodID setToolTip = env->GetMethodID(cls,
					"setToolTip", "(Ljava/lang/String;)V");
			jmethodID setOptions = env->GetMethodID(cls,
					"setOptions", "(Ljava/lang/String;)V");
			jmethodID setRetValue = env->GetMethodID(cls,
					"setReturnValue", "(Ledu/gcsc/vrl/ug/NativeParamInfo;)V");
			jmethodID setParameters = env->GetMethodID(cls,
					"setParameters", "([Ledu/gcsc/vrl/ug/NativeParamInfo;)V");

			using namespace ug::bridge;
			std::string name = method->name(); // TODO pre-rpocessing necessary
			env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
			env->CallVoidMethod(obj, setHelp, stringC2J(env, method->help().c_str()));
			env->CallVoidMethod(obj, setToolTip, stringC2J(env, name.c_str()));
			env->CallVoidMethod(obj, setOptions,
					stringC2J(env, method->options().c_str()));
			env->CallVoidMethod(obj, setRetValue,
					retVal2NativeParam(env, *method));

			env->CallVoidMethod(obj, setParameters,
					params2NativeParams(env, *method));

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

			jobjectArray result =
					env->NewObjectArray(numMethodGroups, groupCls, 0);

			unsigned int numberOfMethodsInGroup = 1;

			for (size_t i = 0; i < numMethodGroups; i++) {

				if (constMethods) {
					numberOfMethodsInGroup = eCls.num_const_overloads(i);
				} else {
					numberOfMethodsInGroup = eCls.num_overloads(i);
				}

				jobjectArray methodArray =
						env->NewObjectArray(numberOfMethodsInGroup, cls, 0);

				jmethodID methodID = env->GetMethodID(groupCls, "<init>", "()V");
				jobject groupObj = env->NewObject(groupCls, methodID);
				jmethodID setOverloads = env->GetMethodID(groupCls,
						"setOverloads", "([Ledu/gcsc/vrl/ug/NativeMethodInfo;)V");

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
					"setReturnValue", "(Ledu/gcsc/vrl/ug/NativeParamInfo;)V");
			jmethodID setParameters = env->GetMethodID(cls,
					"setParameters", "([Ledu/gcsc/vrl/ug/NativeParamInfo;)V");

			using namespace ug::bridge;
			std::string name = func.name(); // TODO pre-rpocessing necessary
			env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
			env->CallVoidMethod(obj, setConst, boolC2J(false));
			env->CallVoidMethod(obj, setCategory,
					stringC2J(env, func.group().c_str()));
			env->CallVoidMethod(obj, setHelp, stringC2J(env, func.help().c_str()));
			env->CallVoidMethod(obj, setToolTip, stringC2J(env, name.c_str()));
			env->CallVoidMethod(obj, setOptions,
					stringC2J(env, func.options().c_str()));
			env->CallVoidMethod(obj, setRetValue,
					retVal2NativeParam(env, func));

			env->CallVoidMethod(obj, setParameters,
					params2NativeParams(env, func));

			return obj;
		}

		jobjectArray functions2NativeGroups(JNIEnv *env, ug::bridge::Registry* reg) {
			jclass groupArrayCls =
					env->FindClass("edu/gcsc/vrl/ug/NativeFunctionGroupInfo");

			unsigned int numFunctions = reg->num_functions();

			// create array of functions
			jobjectArray result =
					env->NewObjectArray(numFunctions, groupArrayCls, 0);

			for (size_t i = 0; i < numFunctions; i++) {

				const ug::bridge::ExportedFunctionGroup& group =
						reg->get_function_group(i);

				size_t numOverloads = group.num_overloads();

				jclass functionsCls =
						env->FindClass("edu/gcsc/vrl/ug/NativeFunctionInfo");

				jobjectArray functions =
						env->NewObjectArray(numOverloads, functionsCls, 0);

				for (size_t j = 0; j < numOverloads; j++) {
					const ug::bridge::ExportedFunction& func = *group.get_overload(j);

					env->SetObjectArrayElement(functions, j,
							function2NativeFunction(env, func));

				} // end for j

				// create function group instance
				jmethodID methodID = env->GetMethodID(groupArrayCls, "<init>", "()V");
				jobject obj = env->NewObject(groupArrayCls, methodID);

				jmethodID setOverloads = env->GetMethodID(groupArrayCls,
						"setOverloads", "([Ledu/gcsc/vrl/ug/NativeFunctionInfo;)V");

				env->CallVoidMethod(obj, setOverloads, functions);

				env->SetObjectArrayElement(result, i, obj);
			} // end for i

			return result;
		}

		jobjectArray classes2NativeClasses(JNIEnv *env,
				const ug::bridge::Registry* reg) {

			jclass cls = env->FindClass("edu/gcsc/vrl/ug/NativeClassInfo");

			jobjectArray result =
					env->NewObjectArray(reg->num_classes(), cls, 0);

			for (size_t i = 0; i < reg->num_classes(); i++) {

				const ug::bridge::IExportedClass& eCls = reg->get_class(i);

				// create instance

				jmethodID methodID = env->GetMethodID(cls, "<init>", "()V");
				jobject obj = env->NewObject(cls, methodID);

				// assign values

				jmethodID setName = env->GetMethodID(cls,
						"setName", "(Ljava/lang/String;)V");
				jmethodID setCategory = env->GetMethodID(cls,
						"setCategory", "(Ljava/lang/String;)V");
				jmethodID setClassNames = env->GetMethodID(
						cls, "setClassNames", "([Ljava/lang/String;)V");
				jmethodID setInstantiable = env->GetMethodID(cls,
						"setInstantiable", "(Z)V");
				jmethodID setConstructors = env->GetMethodID(cls,
						"setConstructors", "([Ledu/gcsc/vrl/ug/NativeConstructorInfo;)V");
				jmethodID setMethods = env->GetMethodID(cls,
						"setMethods", "([Ledu/gcsc/vrl/ug/NativeMethodGroupInfo;)V");
				jmethodID setConstMethods = env->GetMethodID(cls,
						"setConstMethods",
						"([Ledu/gcsc/vrl/ug/NativeMethodGroupInfo;)V");

				std::string name = eCls.name(); // TODO pre-rpocessing necessary

				// these lines check for empty class names. we really want this exit()
				// command as empty names will mess up everything.
				std::vector<std::string> baseClasses;

				for (size_t j = 0; j < eCls.class_names()->size(); j++) {

					if (eCls.class_names()->at(j) == NULL) {
						std::cerr << name << ", baseCls(" << j << ")==NULL" << std::endl;
						exit(1);
					}

					if (strlen(eCls.class_names()->at(j)) == 0) {
						std::cerr << name << ", baseCls(" << j << ")==empty" << std::endl;
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

			jobjectArray result =
					env->NewObjectArray(reg->num_class_groups(), cls, 0);

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
				jmethodID setName = env->GetMethodID(cls,
						"setName", "(Ljava/lang/String;)V");
				jmethodID setClasses = env->GetMethodID(
						cls, "setClasses", "([Ljava/lang/String;)V");
				jmethodID setDefaultClass = env->GetMethodID(
						cls, "setDefaultClass", "(Ljava/lang/String;)V");

				// calls the java methods
				std::string name = clsGrp->name();
				env->CallVoidMethod(obj, setName, stringC2J(env, name.c_str()));
				env->CallVoidMethod(obj, setClasses,
						stringArrayC2J(env, class_names));
				env->CallVoidMethod(obj, setDefaultClass,
						stringC2J(env, defaultClassName.c_str()));

				// set array element
				env->SetObjectArrayElement(result, i, obj);
			}

			//		env->ExceptionCheck();
			//		env->ExceptionDescribe();

			return result;
		}

		jobject registry2NativeAPI(JNIEnv *env, ug::bridge::Registry* reg) {

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

			env->CallVoidMethod(obj, setClassGroups, classGroups2NativeClassGroups(env, reg));
			env->CallVoidMethod(obj, setClasses, classes2NativeClasses(env, reg));
			env->CallVoidMethod(obj, setFunctions, functions2NativeGroups(env, reg));

			return obj;
		}


	} // end vrl::
}// end ug::
