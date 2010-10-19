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

namespace ug {
	namespace vrl {
		static ug::bridge::Registry* vrlRegistry = NULL;

		void SetVRLRegistry(ug::bridge::Registry* pReg) {
			vrlRegistry = pReg;
		}
	} // end vrl::
}// end ug::


typedef int (*FunctionPtr)(int, int);

std::string TestHello() {
	return "Hello, World!";
}

int TestAdd(int a, int b) {
	return a + b;
}

int TestMult(int a, int b) {
	return a*b;
}

void* GetTestMult() {
	return (void*) TestMult;
}

void* GetTestAdd() {
	return (void*) TestAdd;
}

bool TestBoolAnd(bool a, bool b) {
	return a && b;
}

double TestDoubleAdd(double a, double b) {
	return a + b;
}

std::string TestStringAdd1(std::string a, std::string b) {
	return a;
}

std::string TestStringAdd2(std::string a, std::string b) {
	return b;
}

std::string TestStringAdd3(std::string a, std::string b, std::string c) {
	return a + b + c;
}

int TestFuncPointer(void* func, int a, int b) {
	FunctionPtr function = (FunctionPtr) func;
	return function(a, b);
}

class TestClass1 {
public:

	TestClass1() {
		//
	}

	std::string hello(int a, int b) {
		std::stringstream result;
		result << "Hello, World! I'm a class!" << a << ", " << b;
		return result.str();
	}
};



//*********************************************************
//* JNI METHODS
//*********************************************************

JNIEXPORT jint JNICALL Java_edu_gcsc_vrl_ug4_UG4_ugInit
(JNIEnv *env, jobject obj, jobjectArray args) {
	std::vector<std::string> arguments = ug::vrl::stringArrayJ2C(env, args);

	char* argv[arguments.size()];
	for (unsigned int i = 0; i < arguments.size(); i++) {
		argv[i] = (char*) arguments[i].c_str();
	}


	static ug::bridge::Registry testReg;
	testReg.add_class_<TestClass1 > ("TestClass")
			.add_constructor()
			.add_method("hello", &TestClass1::hello);

	//	testReg.add_function("UGAddInt", &TestAdd);
	//	testReg.add_function("UGMultInt", &TestMult);
	//	testReg.add_function("TestFuncPointer", &TestFuncPointer);
	//	testReg.add_function("UGHello", &TestHello);
	//	testReg.add_function("GetUgAddInt", &GetTestAdd);
	//	testReg.add_function("GetUgMultInt", &GetTestMult);
	//	testReg.add_function("UGAndBool", &TestBoolAnd);
	//	testReg.add_function("UGDouble", &TestDoubleAdd);
	//testReg.add_function("UGAddString1", &TestStringAdd1);
	//testReg.add_function("UGAddString2", &TestStringAdd2);
	//	testReg.add_function("UGAddString3", &TestStringAdd3);

	int retVal = ug::UGInit(arguments.size(), argv);

	ug::bridge::RegisterTestInterface(testReg);

	//	ug::vrl::SetVRLRegistry(&ug::GetUGRegistry());
	ug::vrl::SetVRLRegistry(&testReg);

	return (jint) retVal;
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug4_UG4_invokeMethod
(JNIEnv *env, jobject obj,
		jstring exportedClassName, jlong objPtr, jboolean readOnly,
		jstring methodName, jobjectArray params) {
	//	ug::bridge::IExportedClass* clazz = (ug::bridge::IExportedClass*) objPtr;

	ug::bridge::IExportedClass* clazz = (ug::bridge::IExportedClass*)
			ug::vrl::getExportedClassPtrByName(
			env, ug::vrl::vrlRegistry, ug::vrl::stringJ2C(env, exportedClassName));

	ug::bridge::ParameterStack paramsIn;
	ug::bridge::ParameterStack paramsOut;

	VRL_DBG("BEFORE_GET_METHOD 1", 1);

	std::string name = ug::vrl::stringJ2C(env, methodName);

	VRL_DBG("BEFORE_GET_METHOD 2", 1);

	const ug::bridge::ExportedMethod* method =
			ug::vrl::getMethodBySignature(
			env, clazz, ug::vrl::boolJ2C(readOnly), name, params);

	VRL_DBG("AFTER_GET_METHOD", 1);

	VRL_DBG(method->name(), 1);

	VRL_DBG("BEFORE_OBJECT_ARRAY_TO_STACK", 1);

	ug::vrl::jobjectArray2ParamStack(env, paramsIn, method->params_in(), params);

	VRL_DBG("AFTER_OBJECT_ARRAY_TO_STACK", 1);

	//	ug::vrl::displayMessage("Test-Message",">> Hello from UG4",ug::vrl::INFO);


	try {

		method->execute((void*) objPtr, paramsIn, paramsOut);

	} catch (ug::bridge::ERROR_IncompatibleClasses ex) {
		std::cout << "Incopatible Conversion from " << ex.m_from << " : " << ex.m_to << std::endl;
	}

	jobject result = NULL;

	if (paramsOut.size() > 0) {
		result = ug::vrl::param2JObject(env, paramsOut, 0);
	}
	return result;
}

JNIEXPORT jlong JNICALL Java_edu_gcsc_vrl_ug4_UG4_newInstance
(JNIEnv *env, jobject obj, jlong objPtr) {
	ug::bridge::IExportedClass* clazz = (ug::bridge::IExportedClass*) objPtr;
	return (long) clazz->create();
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug4_UG4_invokeFunction
(JNIEnv *env, jobject obj, jlong fPtr, jboolean readOnly, jobjectArray params) {
	ug::bridge::ExportedFunction* func = (ug::bridge::ExportedFunction*) fPtr;

	ug::bridge::ParameterStack paramsIn;
	ug::bridge::ParameterStack paramsOut;

	ug::vrl::jobjectArray2ParamStack(env, paramsIn, func->params_in(), params);

	func->execute(paramsIn, paramsOut);

	jobject result = NULL;

	if (paramsOut.size() > 0) {
		result = ug::vrl::param2JObject(env, paramsOut, 0);
	}
	return result;
}

JNIEXPORT jobjectArray JNICALL Java_edu_gcsc_vrl_ug4_UG4_createJavaBindings
(JNIEnv *env, jobject obj) {

	std::vector<std::string> result;

	VRL_DBG("BEFORE_FUNCTIONS",1);

	for (unsigned int i = 0; i < ug::vrl::vrlRegistry->num_functions(); i++) {
		ug::bridge::ExportedFunction& func = ug::vrl::vrlRegistry->get_function(i);
		result.push_back(ug::vrl::exportedFunction2Groovy(func));
	}

	VRL_DBG("FUNCTIONS_DONE",1);

	for (unsigned int i = 0; i < ug::vrl::vrlRegistry->num_classes(); i++) {
		const ug::bridge::IExportedClass& clazz = ug::vrl::vrlRegistry->get_class(i);
		result.push_back(ug::vrl::exportedClass2Groovy(clazz));
	}

	VRL_DBG("CLASSES_DONE",1);

	return ug::vrl::stringArrayC2J(env, result);
}

JNIEXPORT jlong JNICALL Java_edu_gcsc_vrl_ug4_UG4_getExportedClassPtrByName
(JNIEnv *env, jobject obj, jstring name) {
	return (long) ug::vrl::getExportedClassPtrByName(
			env, ug::vrl::vrlRegistry, ug::vrl::stringJ2C(env, name));
}

JNIEXPORT void JNICALL Java_edu_gcsc_vrl_ug4_UG4_attachCanvas
(JNIEnv *env, jobject obj, jobject canvas) {
	ug::vrl::Canvas::getInstance()->setJObject(env, canvas);

	//	ug::vrl::Canvas::getInstance()->addObject(ug::vrl::string2JObject(env,"Test_String"));
}