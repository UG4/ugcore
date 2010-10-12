#include <string>
#include <vector>
#include <string>

#include "ug.h"
#include "ug_bridge/registry.h"
#include "ug_bridge/class.h"

#include "type_converter.h"
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
	return a+b;
}

int TestMult(int a, int b) {
	return a*b;
}

void* GetTestMult() {
	return (void*)TestMult;
}

void* GetTestAdd() {
	return (void*)TestAdd;
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
	return a+b+c;
}

int TestFuncPointer(void* func, int a, int b) {
	FunctionPtr function = (FunctionPtr)func;
	return function(a,b);
}



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
	testReg.add_class_<FunctionPtr>("FunctionPtr");

	testReg.add_function("UGAddInt", &TestAdd);
	testReg.add_function("UGMultInt", &TestMult);
	testReg.add_function("TestFuncPointer", &TestFuncPointer);
	testReg.add_function("UGHello", &TestHello);
	testReg.add_function("GetUgAddInt", &GetTestAdd);
	testReg.add_function("GetUgMultInt", &GetTestMult);
	testReg.add_function("UGAndBool", &TestBoolAnd);
	testReg.add_function("UGDouble", &TestDoubleAdd);
	//testReg.add_function("UGAddString1", &TestStringAdd1);
	//testReg.add_function("UGAddString2", &TestStringAdd2);
	testReg.add_function("UGAddString3", &TestStringAdd3);

	int retVal = ug::UGInit(arguments.size(), argv);
	//ug::vrl::SetVRLRegistry(&ug::GetUGRegistry());
	ug::vrl::SetVRLRegistry(&testReg);

	return (jint) retVal;
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug4_UG4_invokeFunction
(JNIEnv *env, jobject obj, jlong fPtr, jobjectArray params) {
	ug::bridge::ExportedFunction* func = (ug::bridge::ExportedFunction*) fPtr;

	ug::bridge::ParameterStack paramsIn;
	ug::bridge::ParameterStack paramsOut;

	std::cout << "PARAMS_BEFORE" << std::endl;

	ug::vrl::jobjectArray2ParamStack(env, paramsIn, func->params_in(), params);

	std::cout << "PARAMS_AFTER" << std::endl;

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

	for (unsigned int i = 0; i < ug::vrl::vrlRegistry->num_functions(); i++) {
		ug::bridge::ExportedFunction& func = ug::vrl::vrlRegistry->get_function(i);
		result.push_back(ug::vrl::exportedFunction2Groovy(func));
	}

	return ug::vrl::stringArrayC2J(env, result);
}
