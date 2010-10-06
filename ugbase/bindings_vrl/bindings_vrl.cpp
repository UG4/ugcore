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

int TestAdd(int a, int b) {
	return a + b;
}

//*********************************************************
//* JNI METHODS
//*********************************************************

JNIEXPORT jobjectArray JNICALL Java_edu_gcsc_vrl_ug4_UG4_helloUG
(JNIEnv *env, jobject obj) {
	int argc = 0;
	char* argv[argc];

	ug::UGInit(argc, argv);
	ug::vrl::SetVRLRegistry(&ug::GetUGRegistry());

	std::vector<std::string> classes;
	classes.push_back("class Hello_From_UG4{String hello(){return \"Hello from UG4!\"}}");
	classes.push_back("***TESTING***");

	const std::vector< const char*> * classNames = ug::vrl::vrlRegistry->get_class(0).class_names();

	for (unsigned int i = 0; i < classNames->size(); i++) {
		classes.push_back(std::string((*classNames)[i]));
	}

	return ug::vrl::stringArrayC2J(env, classes);
}

JNIEXPORT jint JNICALL Java_edu_gcsc_vrl_ug4_UG4_ugInit
(JNIEnv *env, jobject obj, jobjectArray args) {
	std::vector<std::string> arguments = ug::vrl::stringArrayJ2C(env, args);

	char* argv[arguments.size()];
	for (unsigned int i = 0; i < arguments.size(); i++) {
		argv[i] = (char*) arguments[i].c_str();
	}


	static ug::bridge::Registry testReg;
	testReg.add_function("UGAdd", &TestAdd);

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

	for (unsigned int i = 0; i < ug::vrl::vrlRegistry->num_functions(); i++) {
		ug::bridge::ExportedFunction& func = ug::vrl::vrlRegistry->get_function(i);
		result.push_back(ug::vrl::exportedFunction2Groovy(func));
	}

	return ug::vrl::stringArrayC2J(env, result);
}
