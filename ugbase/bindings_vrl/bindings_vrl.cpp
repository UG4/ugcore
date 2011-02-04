#include <string>
#include <vector>
#include <string>

#include "ug.h"
#include "ug_bridge/registry.h"
#include "ug_bridge/class.h"


#include "type_converter.h"
#include "messaging.h"
#include "canvas.h"
#include "bindings_vrl.h"
#include "bindings_vrl_native.h"

#include "lib_grid/lib_grid.h"
#include "compiledate.h"
#include "vrl_user_number.h"
#include "invocation.h"

namespace ug {
namespace vrl {
static ug::bridge::Registry* vrlRegistry = NULL;
static JNIEnv* jniEnv = NULL;
static JavaVM* javaVM = NULL;

void SetVRLRegistry(ug::bridge::Registry* pReg) {
	vrlRegistry = pReg;
}

void SetJNIEnv(JNIEnv* env) {
	jniEnv = env;
	env->GetJavaVM(&javaVM);
}

JNIEnv* getJNIEnv() {
	return jniEnv;
}

JavaVM* getJavaVM() {
	return javaVM;
}

} // end vrl::
}// end ug::

class TestClass {
public:

	TestClass() {
		//
	}

	std::string getRev() {
		return ug::vrl::svnRevision();
	}

	int add(int a, int b) {
		return a + b;
	}

	std::string getString() {
		UG_LOG("Test123" << std::endl);
		return "Test123";
	}

	~TestClass() {
		UG_LOG("Destructor called:" << (long) this << std::endl);
	}
};

//*********************************************************
//* JNI METHODS
//*********************************************************

JNIEXPORT jint JNICALL Java_edu_gcsc_vrl_ug4_UG4_ugInit
(JNIEnv *env, jobject obj, jobjectArray args) {

	ug::vrl::SetJNIEnv(env);

	std::vector<std::string> arguments = ug::vrl::stringArrayJ2C(env, args);

	std::vector<char*> argv(arguments.size());
	for (unsigned int i = 0; i < arguments.size(); i++) {
		argv[i] = (char*) arguments[i].c_str();
	}

	//	static ug::bridge::Registry testReg;

	//	Choose registry used.
	ug::bridge::Registry& reg = ug::bridge::GetUGRegistry();
	//	ug::bridge::Registry& reg = testReg;

	using namespace ug;

	int argc = arguments.size();
	char** pargv = &argv[0];
	int retVal = ug::UGInit(&argc, &pargv);


//	reg.add_class_<TestClass > ("TestClass", "testing")
//			.add_constructor()
//			.add_method("svnRevision", &TestClass::getRev)
//			.add_method("add", &TestClass::add, "result",
//			"a|default|min=-3;max=5;value=-12#b|default|min=-1;max=1;value=23")
//			.add_method("getString", &TestClass::getString);

	//	Register Standard Interfaces (excluding algebra)
	//		ug::bridge::RegisterStandardInterfaces(reg);

	//	Register algebra
	CPUAlgebraChooser chooser;
	ug::bridge::RegisterDynamicLibAlgebraInterface(reg, chooser.get_algebra_type());
	ug::bridge::RegisterDynamicLibDiscretizationInterface(reg, chooser.get_algebra_type());

//	ug::vrl::RegisterVRLUserNumber(reg, "testing");
	//				ug::bridge::RegisterTestInterface(reg);

	//	ug::bridge::RegisterLibGridInterface(testReg);

	//	ug::vrl::SetVRLRegistry(&ug::GetUGRegistry());
	ug::vrl::SetVRLRegistry(&reg);

	return (jint) retVal;
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug4_UG4_invokeMethod
(JNIEnv *env, jobject obj,
		jstring exportedClassName, jlong objPtr, jboolean readOnly,
		jstring methodName, jobjectArray params) {

	const ug::bridge::IExportedClass* clazz =
			ug::vrl::invocation::getExportedClassPtrByName(
			ug::vrl::vrlRegistry, ug::vrl::stringJ2C(env, exportedClassName));

	ug::bridge::ParameterStack paramsIn;
	ug::bridge::ParameterStack paramsOut;

	std::string name = ug::vrl::stringJ2C(env, methodName);

	jobject result = NULL;

	try {

		const ug::bridge::ExportedMethod* method =
				ug::vrl::invocation::getMethodBySignature(
				env, ug::vrl::vrlRegistry,
				clazz, ug::vrl::boolJ2C(readOnly), name, params);

		if (method == NULL && readOnly == false) {
			method = ug::vrl::invocation::getMethodBySignature(
					env, ug::vrl::vrlRegistry,
					clazz, ug::vrl::boolJ2C(true), name, params);
		}

		ug::vrl::jobjectArray2ParamStack(
				env, paramsIn, method->params_in(), params);

		method->execute((void*) objPtr, paramsIn, paramsOut);


		if (paramsOut.size() > 0) {
			result = ug::vrl::param2JObject(env, paramsOut, 0);
		}

	} catch (ug::bridge::ERROR_IncompatibleClasses ex) {
		UG_LOG("Incopatible Conversion from " <<
				ex.m_from << " : " << ex.m_to << std::endl);
	} catch (ug::bridge::ERROR_BadConversion ex) {
		UG_LOG("Incopatible Conversion from " <<
				ex.m_from << " : " << ex.m_to << std::endl);
	} catch (...) {
		UG_LOG("Unknown exception thrown while"
				<< " trying to invoke method!" << std::endl);
	}

	return result;
}

JNIEXPORT jlong JNICALL Java_edu_gcsc_vrl_ug4_UG4_newInstance
(JNIEnv *env, jobject obj, jlong objPtr) {

	long result = NULL;

	ug::bridge::IExportedClass* clazz = NULL;

	try {
		ug::bridge::IExportedClass* clazz =
				(ug::bridge::IExportedClass*) objPtr;
		result = (long) clazz->create();
	} catch (...) {
		std::string className = "Unknown class";
		if (clazz != NULL) {
			className = clazz->name();
		}

		UG_LOG("Unknown exception thrown while"
				<< " trying to instanciate class \""
				<< className << "\"!" << std::endl);
	}

	return result;
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug4_UG4_invokeFunction
(JNIEnv *env, jobject obj, jstring fName, jboolean readOnly, jobjectArray params) {



	const ug::bridge::ExportedFunction* func =
			ug::vrl::invocation::getFunctionBySignature(
			env, ug::vrl::vrlRegistry, ug::vrl::stringJ2C(env, fName), params);

	ug::bridge::ParameterStack paramsIn;
	ug::bridge::ParameterStack paramsOut;

	jobject result = NULL;

	try {
		ug::vrl::jobjectArray2ParamStack(
				env, paramsIn, func->params_in(), params);

		func->execute(paramsIn, paramsOut);

		if (paramsOut.size() > 0) {
			result = ug::vrl::param2JObject(env, paramsOut, 0);
		}

	} catch (ug::bridge::ERROR_IncompatibleClasses ex) {
		UG_LOG("Incopatible Conversion from " <<
				ex.m_from << " : " << ex.m_to << std::endl);
	} catch (ug::bridge::ERROR_BadConversion ex) {
		UG_LOG("Incopatible Conversion from " <<
				ex.m_from << " : " << ex.m_to << std::endl);
	} catch (...) {
		UG_LOG("Unknown exception thrown while"
				<< " trying to invoke function!" << std::endl);
	}

	return result;
}

JNIEXPORT jobjectArray JNICALL Java_edu_gcsc_vrl_ug4_UG4_createJavaBindings
(JNIEnv *env, jobject obj) {

	std::vector<std::string> cResult;
	jobjectArray jResult =
			ug::vrl::createEmptyJavaArray(env, "java/lang/String");

	try {

		for (unsigned int i = 0; i < ug::vrl::vrlRegistry->num_classes(); i++) {

			const ug::bridge::IExportedClass& clazz =
					ug::vrl::vrlRegistry->get_class(i);

			if (clazz.is_instantiable()) {
				cResult.push_back(
						ug::vrl::exportedClass2Groovy(
						ug::vrl::vrlRegistry, clazz));
			}

		}

		for (unsigned int i = 0; i < ug::vrl::vrlRegistry->num_functions(); i++) {
			ug::bridge::ExportedFunction& func =
					ug::vrl::vrlRegistry->get_function(i);
			cResult.push_back(ug::vrl::exportedFunction2Groovy(func));
		}

		jResult = ug::vrl::stringArrayC2J(env, cResult);

	} catch (...) {
		UG_LOG("Unknown exception thrown while"
				<< " trying to convert registered classes to Groovy code!"
				<< std::endl);
	}

	// print groovy code
	//	for(unsigned int i = 0; i < cResult.size();i++) {
	//		std::cout << "\n\n***********************************************\n\n";
	//		std::cout << cResult[i];
	//	}

	return jResult;
}

JNIEXPORT jlong JNICALL Java_edu_gcsc_vrl_ug4_UG4_getExportedClassPtrByName
(JNIEnv *env, jobject obj, jstring name) {
	return (long) ug::vrl::invocation::getExportedClassPtrByName(
			ug::vrl::vrlRegistry, ug::vrl::stringJ2C(env, name));
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug4_UG4_getMessages
(JNIEnv *env, jobject obj) {
	return ug::vrl::stringC2J(env, ug::vrl::MessageBuffer::getMessages().c_str());
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug4_UG4_getSvnRevision
(JNIEnv *env, jobject obj) {
	std::string revision = ug::vrl::svnRevision();
	return ug::vrl::stringC2J(env, revision.c_str());
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug4_UG4_getCompileDate
(JNIEnv *env, jobject obj) {
	return ug::vrl::stringC2J(env, COMPILE_DATE);
}

JNIEXPORT void JNICALL Java_edu_gcsc_vrl_ug4_MemoryManager_delete
(JNIEnv * env, jclass cls, jlong objPtr, jlong exportedClsPtr) {

//	if (((void*)objPtr) != NULL && ((void*)exportedClsPtr) != NULL) {
//
//		ug::bridge::IExportedClass* clazz =
//				(ug::bridge::IExportedClass*) exportedClsPtr;
//		clazz->destroy((void*) objPtr);
//	}
}

//JNIEXPORT void JNICALL Java_edu_gcsc_vrl_ug4_UG4_attachCanvas
//(JNIEnv *env, jobject obj, jobject canvas) {
//	ug::vrl::Canvas::getInstance()->setJObject(env, canvas);
//
//	//	ug::vrl::Canvas::getInstance()->addObject(ug::vrl::string2JObject(env,"Test_String"));
//}
