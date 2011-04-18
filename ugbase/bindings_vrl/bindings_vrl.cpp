#include <string>
#include <vector>
#include <string>

#include "ug.h"
#include "ugbase.h"
#include "ug_bridge/registry.h"
#include "ug_bridge/class.h"


#include "type_converter.h"
#include "messaging.h"
#include "canvas.h"
#include "bindings_vrl.h"
#include "bindings_vrl_native.h"

#include "lib_grid/lib_grid.h"
#include "compiledate.h"
#include "user_data.h"
#include "invocation.h"
#include "playground.h"

namespace ug {
namespace vrl {
static ug::bridge::Registry* vrlRegistry = NULL;
static JavaVM* javaVM = NULL;

void SetVRLRegistry(ug::bridge::Registry* pReg) {
	vrlRegistry = pReg;
}

void initJavaVM(JNIEnv* env) {
	if (javaVM == NULL) {
		env->GetJavaVM(&javaVM);
	} else {
		UG_LOG("UG-VRL: JavaVM already initialized!"
				" JavaVM can be initialized only once!");
	}
}

JavaVM* getJavaVM() {
	return javaVM;
}

}// end vrl::
}// end ug::


//*********************************************************
//* JNI METHODS
//*********************************************************

JNIEXPORT jint JNICALL Java_edu_gcsc_vrl_ug4_UG4_ugInit
(JNIEnv *env, jobject obj, jobjectArray args) {

	ug::vrl::initJavaVM(env);

	std::vector<std::string> arguments = ug::vrl::stringArrayJ2C(env, args);

	std::vector<char*> argv(arguments.size());
	for (unsigned int i = 0; i < arguments.size(); i++) {
		argv[i] = (char*) arguments[i].c_str();
	}

	// Choose registry used.
	ug::bridge::Registry& reg = ug::bridge::GetUGRegistry();

	using namespace ug;

	int argc = arguments.size();
	char** pargv = &argv[0];
	int retVal = ug::UGInit(&argc, &pargv);

	// Register Playground if we are in debug mode

	//#ifdef UG_DEBUG
	//	registerPlayground(reg);
	//#endif

	// Register algebra
	CPUAlgebraSelector selector;
	ug::bridge::RegisterDynamicLibAlgebraInterface(
			reg, selector.get_algebra_type());
	ug::bridge::RegisterDynamicLibDiscretizationInterface(
			reg, selector.get_algebra_type());

	ug::vrl::RegisterUserData(reg,"UG4/VRL");


	if (!reg.check_consistency()) {
		UG_LOG("UG-VRL: cannot compile code due to registration error.");
		return 1;
	}

	ug::vrl::SetVRLRegistry(&reg);

	return (jint) retVal;
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug4_UG4_invokeMethod
(JNIEnv *env, jobject obj,
		jstring exportedClassName, jlong objPtr, jboolean readOnly,
		jstring methodName, jobjectArray params) {

	std::string className = ug::vrl::stringJ2C(env, exportedClassName);

	const ug::bridge::IExportedClass* clazz =
			ug::vrl::invocation::getExportedClassPtrByName(
			ug::vrl::vrlRegistry, className);

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

		if (method == NULL) {
			UG_LOG("Method not found: " << name <<
					"()" << " : " << std::endl << VRL_CRITICAL_ERROR);
			return NULL;
		}

		ug::vrl::jobjectArray2ParamStack(
				env, paramsIn, method->params_in(), params);

		const ug::bridge::ClassNameNode* clsNode =
				ug::vrl::invocation::getClassNodePtrByName(
				ug::vrl::vrlRegistry, className);

		void* finalObjPtr = ug::bridge::ClassCastProvider::cast_to_base_class(
				(void*) objPtr,
				clsNode, method->class_name());

		method->execute(finalObjPtr, paramsIn, paramsOut);

		if (paramsOut.size() > 0) {
			result = ug::vrl::param2JObject(env, paramsOut, 0);
		}

	} catch (ug::bridge::ERROR_IncompatibleClasses ex) {
		UG_LOG("Incompatible Conversion from " <<
				ex.m_from << " : " << ex.m_to << std::endl
				<< VRL_CRITICAL_ERROR);
	} catch (ug::bridge::ERROR_BadConversion ex) {
		UG_LOG("Incompatible Conversion from " <<
				ex.m_from << " : " << ex.m_to << std::endl
				<< VRL_CRITICAL_ERROR);
	} catch (...) {
		UG_LOG("Unknown exception thrown while"
				<< " trying to invoke method!" << std::endl
				<< VRL_CRITICAL_ERROR);
	}

	return result;
}

JNIEXPORT jlong JNICALL Java_edu_gcsc_vrl_ug4_UG4_newInstance
(JNIEnv *env, jobject obj, jlong objPtr) {

	long result = 0;
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

		if (func == NULL) {
			UG_LOG("Function not found: " << ug::vrl::stringJ2C(env, fName) <<
					"()" << " : " << std::endl << VRL_CRITICAL_ERROR);
			return NULL;
		}

		ug::vrl::jobjectArray2ParamStack(
				env, paramsIn, func->params_in(), params);

		func->execute(paramsIn, paramsOut);

		if (paramsOut.size() > 0) {
			result = ug::vrl::param2JObject(env, paramsOut, 0);
		}

	} catch (ug::bridge::ERROR_IncompatibleClasses ex) {
		UG_LOG("Incopatible Conversion from " <<
				ex.m_from << " : " << ex.m_to << std::endl
				<< VRL_CRITICAL_ERROR);
	} catch (ug::bridge::ERROR_BadConversion ex) {
		UG_LOG("Incopatible Conversion from " <<
				ex.m_from << " : " << ex.m_to << std::endl
				<< VRL_CRITICAL_ERROR);
	} catch (...) {
		UG_LOG("Unknown exception thrown while"
				<< " trying to invoke function!" << std::endl);
	}

	return result;
}

JNIEXPORT jlong JNICALL Java_edu_gcsc_vrl_ug4_UG4_getExportedClassPtrByName
(JNIEnv *env, jobject obj, jstring name) {
	return (long) ug::vrl::invocation::getExportedClassPtrByName(
			ug::vrl::vrlRegistry, ug::vrl::stringJ2C(env, name));
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

	if (((void*) objPtr) != NULL && ((void*) exportedClsPtr) != NULL) {
		ug::bridge::IExportedClass* clazz =
				(ug::bridge::IExportedClass*) exportedClsPtr;
		clazz->destroy((void*) objPtr);
	}
}

JNIEXPORT void JNICALL Java_edu_gcsc_vrl_ug4_MemoryManager_invalidate
(JNIEnv * env, jclass cls, jobject smartPtr) {

	if (ug::vrl::isJSmartPointerConst(env, smartPtr)) {
		ug::vrl::invalidateJConstSmartPointer(env, smartPtr);
	} else {
		ug::vrl::invalidateJSmartPointer(env, smartPtr);
	}
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug4_UG4_convertRegistryInfo
(JNIEnv * env, jobject obj) {
	ug::vrl::invocation::initClasses(*ug::vrl::vrlRegistry);
	return ug::vrl::registry2NativeAPI(env, ug::vrl::vrlRegistry);
}



//JNIEXPORT void JNICALL Java_edu_gcsc_vrl_ug4_UG4_attachCanvas
//(JNIEnv *env, jobject obj, jobject canvas) {
//	ug::vrl::Canvas::getInstance()->setJObject(env, canvas);
//
//	//	ug::vrl::Canvas::getInstance()->addObject(ug::vrl::string2JObject(env,"Test_String"));
//}
