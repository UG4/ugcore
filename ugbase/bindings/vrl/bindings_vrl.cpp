#include <string>
#include <vector>

#include "ug.h"
#include "ugbase.h"
#include "registry/registry.h"
#include "registry/class.h"
#include "common/util/path_provider.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

#include "common/common.h"
#include "lib_algebra/operator/convergence_check.h"
#include "common/authors.h"
#include "common/util/string_util.h"

#include "type_converter.h"
#include "messaging.h"
#include "canvas.h"
#include "bindings_vrl.h"
#include "bindings_vrl_native.h"

#include "user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

#include "compile_info/compile_info.h"

#include "invocation.h"
#include "playground.h"
#include "threading.h"

//#include "bindings/lua/externals/lua/lstate.h"
#include "basicTest.h"

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

void Log(std::string s) {
	UG_LOG(s);
}

void Logln(std::string s) {
	UG_LOG(s << std::endl);
}

void ThrowIf(bool b, std::string s) {
	if (!b) {
		throw(ug::UGError(s.c_str()));
	}
}

void ThrowIfNot(bool b, std::string s) {
	if (!b) {
		throw(ug::UGError(s.c_str()));
	}
}

void registerMessaging(ug::bridge::Registry & reg) {
	reg.add_function("print", &Log, "UG4/Util/Messaging");
	reg.add_function("println", &Logln, "UG4/Util/Messaging");
}

void registerThrowUtil(ug::bridge::Registry & reg) {
	reg.add_function("throwIf", &ThrowIf, "UG4/Util");
	reg.add_function("throwIfNot", &ThrowIfNot, "UG4/Util");
}

class NumberArray {
private:
	std::vector<number> _vec;
public:

	NumberArray() {
	}

	NumberArray(std::vector<number> vec) {
		_vec = vec;
	}

	void setArray(std::vector<number> vec) {
		_vec = vec;
	}

	int size() const {
		return _vec.size();
	}

	number get(int i) const {

		if (i < 0 || (size_t) i >= _vec.size()) {
			throw UGError("NumberArray: index out of Bounds!");
		}

		return _vec[i];
	}
};

template<typename TVector>
SmartPtr<NumberArray> getDefects(const ug::StdConvCheck<TVector>* convCheck) {

	return SmartPtr<NumberArray>(new NumberArray(convCheck->get_defects()));
}

void registerNumberArray(ug::bridge::Registry & reg) {
	reg.add_class_<NumberArray>("NumberArray", "UG4/Util").add_constructor().add_method(
			"get", &NumberArray::get).add_method("size", &NumberArray::size);
}

void registerUGFinalize(ug::bridge::Registry & reg) {
	reg.add_function("UGFinalize", &ug::UGFinalize, "UG4/Util");
}

class VTest {
public:

	VTest() {
		UG_LOG("VTest::VTest() constructor used.\n")
	}

	VTest(const char* msg) {
		UG_LOG("VTest::VTest(const char*) constructor used.\n")
		UG_LOG("Message is: '" << msg << "'.\n");
	}

	std::string hello() {
		return "i am instantiated!";
	}
};

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality {

	/**
	 * Function called for the registration of Algebra dependent parts.
	 * All Functions and Classes depending on Algebra
	 * are to be placed here when registering. The method is called for all
	 * available Algebra types, based on the current build options.
	 *
	 * @param reg				registry
	 * @param parentGroup		group for sorting of functionality
	 */
	template<typename TAlgebra>
	static void Algebra(ug::bridge::Registry& reg, std::string parentGroup) {
//	typedefs for Vector and Matrix
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename TAlgebra::matrix_type matrix_type;

//	suffix and tag
		std::string suffix = ug::bridge::GetAlgebraSuffix<TAlgebra>();
		std::string tag = ug::bridge::GetAlgebraTag<TAlgebra>();

		reg.add_function("GetDefects", &getDefects<vector_type>, "UG4/Util",
				"Defects");
	}

};
// end Functionality

void RegisterVRLFunctionality(ug::bridge::Registry& reg, std::string grp) {
	typedef ug::vrl::Functionality Functionality;

	ug::bridge::RegisterAlgebraDependent<Functionality>(reg, grp);
}

} // end vrl::
} // end ug::


//*********************************************************
//* JNI METHODS
//*********************************************************
JNIEXPORT jint JNICALL Java_edu_gcsc_vrl_ug_UG__1ugInit(JNIEnv *env, jclass cls,
		jobjectArray args) {

	ug::vrl::initJavaVM(env);

	std::vector<std::string> arguments = ug::vrl::stringArrayJ2C(env, args);

	std::vector<char*> argv(arguments.size());
	for (unsigned int i = 0; i < arguments.size(); i++) {
		argv[i] = (char*) arguments[i].c_str();
	}

	// Choose registry used.
	ug::bridge::Registry& reg = ug::bridge::GetUGRegistry();

	//	reg.add_callback(&ug::vrl::registryChanged);

	using namespace ug;

	// define paths
	ug::PathProvider::set_path(PLUGIN_PATH, arguments[0]);

	int argc = arguments.size();
	char** pargv = &argv[0];
	//\todo: generalize outputproc rank
	// isn't this possible already via SetOuputRank() ?
	int retVal = ug::UGInit(&argc, &pargv, 0);

	// Register Playground if we are in debug mode

#ifdef UG_DEBUG
	ug::vrl::registerPlayground(reg);
#endif

	ug::vrl::registerBasicTest(reg);

	ug::vrl::RegisterUserData(reg, "UG4/VRL");
	ug::vrl::registerMessaging(reg);
	ug::vrl::registerThrowUtil(reg);
	ug::vrl::registerNumberArray(reg);
	ug::vrl::RegisterVRLFunctionality(reg, "UG4/VRL");
	ug::vrl::registerUGFinalize(reg);

	reg.add_class_<ug::vrl::VTest>("VTest", "UG4/VRL/Testing").add_constructor().add_constructor<
			void (*)(const char*)>().add_method("hello",
			&ug::vrl::VTest::hello);

	if (!reg.check_consistency()) {
		ug::GetLogAssistant().flush_error_log();
		UG_LOG("UG-VRL: cannot compile code due to registration error.\n");
		return 1;
	}

	ug::vrl::SetVRLRegistry(&reg);

	ug::vrl::invocation::initClasses(*ug::vrl::vrlRegistry);

	if (!ug::vrl::vrlRegistry->check_consistency()) {
		ug::GetLogAssistant().flush_error_log();
		UG_LOG("UG-VRL: cannot compile code due to registration error after initClasses.\n");
		return 1;
	}

	return (jint) retVal;
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug_UG__1invokeMethod(JNIEnv *env,
		jobject obj, jstring exportedClassName, jlong objPtr, jboolean readOnly,
		jstring methodName, jobjectArray params) {

	bool DEBUG = true;
	if (DEBUG) {
		std::cout << "trunk/ugbase/bindings/vrl/bindings_vrl.cpp :"
				<< " Java_edu_gcsc_vrl_ug_UG__1invokeMethod() " << std::endl;
	}

	std::string className = ug::vrl::stringJ2C(env, exportedClassName);

	const ug::bridge::IExportedClass* clazz =
			ug::vrl::invocation::getExportedClassPtrByName(ug::vrl::vrlRegistry,
					className);

	ug::bridge::ParameterStack paramsIn;
	ug::bridge::ParameterStack paramsOut;

	std::string name = ug::vrl::stringJ2C(env, methodName);

	jobject result = NULL;

	try {
		const ug::bridge::ExportedMethod* exMethod =
				ug::vrl::invocation::getMethodBySignature(env,
						ug::vrl::vrlRegistry, clazz, ug::vrl::boolJ2C(readOnly),
						name, params);

		if (exMethod == NULL && readOnly == false) {
			exMethod = ug::vrl::invocation::getMethodBySignature(env,
					ug::vrl::vrlRegistry, clazz, ug::vrl::boolJ2C(true), name,
					params);
		}

		if (exMethod == NULL) {

			std::stringstream ss;

			ss << "No method found that matches the given signature: "
					<< EMPHASIZE_BEGIN << clazz->name() << "." << name
					<< "(" + ug::vrl::getParamTypesAsString(env, params) + ")"
					<< EMPHASIZE_END << ".";

			jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
			env->ThrowNew(Exception, ss.str().c_str());
			return NULL;
		}

		ug::vrl::jobjectArray2ParamStack(env, ug::vrl::vrlRegistry, paramsIn,
				exMethod->params_in(), params);

		const ug::bridge::ClassNameNode* clsNode =
				ug::vrl::invocation::getClassNodePtrByName(ug::vrl::vrlRegistry,
						className);

		void* finalObjPtr = ug::bridge::ClassCastProvider::cast_to_base_class(
				(void*) objPtr, clsNode, exMethod->class_name());

		exMethod->execute(finalObjPtr, paramsIn, paramsOut);

		if (paramsOut.size() > 0) {
			result = ug::vrl::param2JObject(env, paramsOut, 0);
		}

	} catch (ug::bridge::UGError_ClassCastFailed& ex) {

		std::stringstream ss;

		ss << "Incompatible conversion in method " << className << "."
				<< methodName << "(): from " << ex.m_from << " to " << ex.m_to;

		jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
		env->ThrowNew(Exception, ss.str().c_str());
	} catch (ug::UGError& ex) {

		UG_LOG("bindings_vrl.cpp : _1invokeMethod() : catch (ug::UGError& ex)" << std::endl);

		ug::vrl::throwUgErrorAsJavaException(env, ex);
	} catch (...) {

		std::stringstream ss;

		ss << "Unknown exception thrown while" << " trying to invoke method: "
				<< clazz->name() << "." << name << "().";

		jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
		env->ThrowNew(Exception, ss.str().c_str());
	}

	return result;
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug_UG__1newInstance(JNIEnv *env,
		jobject obj, jlong exportedClassPointer, jobjectArray params) {

	bool DEBUG = true;
	if (DEBUG) {
		std::cout << "trunk/ugbase/bindings/vrl/bindings_vrl.cpp :"
				<< " Java_edu_gcsc_vrl_ug_UG__1newInstance() " << std::endl;
	}

	ug::bridge::IExportedClass* clazz =
			(ug::bridge::IExportedClass*) exportedClassPointer;

	ug::bridge::ParameterStack paramsIn;

	std::string name = "constructor";

	try {
		const ug::bridge::ExportedConstructor* constructor =
				ug::vrl::invocation::getConstructorBySignature(env,
						ug::vrl::vrlRegistry, clazz, params);

		if (constructor == NULL) {

			std::stringstream ss;

			ss << "No constructor found that matches the given signature: "
					<< EMPHASIZE_BEGIN << clazz->name() << "." << name
					<< "(" + ug::vrl::getParamTypesAsString(env, params) + ")"
					<< EMPHASIZE_END << ".";

			jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
			env->ThrowNew(Exception, ss.str().c_str());
			return (jlong) NULL;
		}

		ug::vrl::jobjectArray2ParamStack(env, ug::vrl::vrlRegistry, paramsIn,
				constructor->params_in(), params);

		if (clazz->construct_as_smart_pointer()) {
			SmartPtr<void> instance = SmartPtr<void>(
					constructor->create(paramsIn),
					clazz->get_delete_function());

			return ug::vrl::smartPointer2JObject(env, instance);
		} else {
			void *ptr = constructor->create(paramsIn);
			return ug::vrl::pointer2JObject(env, ptr);
		}

	} catch (ug::bridge::UGError_ClassCastFailed& ex) {

		std::stringstream ss;

		ss << "Incompatible conversion in constructor of " << clazz->name()
				<< ": from " << ex.m_from << " to " << ex.m_to;

		jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
		env->ThrowNew(Exception, ss.str().c_str());
	} catch (ug::UGError& ex) {

		ug::vrl::throwUgErrorAsJavaException(env, ex);
	} catch (...) {

		std::stringstream ss;

		ss << "Unknown exception thrown while" << " trying to invoke method: "
				<< name << "().";

		jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
		env->ThrowNew(Exception, ss.str().c_str());
	}

	return (jlong) NULL;
}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug_UG__1invokeFunction(JNIEnv *env,
		jobject obj, jstring fName, jboolean readOnly, jobjectArray params) {

	bool DEBUG = true;
	if (DEBUG) {
		std::cout << "trunk/ugbase/bindings/vrl/bindings_vrl.cpp :"
				<< " Java_edu_gcsc_vrl_ug_UG__1invokeFunction() " << std::endl;
		std::cout << "Java_edu_gcsc_vrl_ug_UG__1invokeFunction() : fName = "
				<< fName << std::endl;
	}

	std::string name = ug::vrl::stringJ2C(env, fName);

	if (DEBUG) {
		std::cout << "trunk/ugbase/bindings/vrl/bindings_vrl.cpp :"
				<< " Java_edu_gcsc_vrl_ug_UG__1invokeFunction() " << std::endl;

		std::cout << "UG__1invokeFunction(): name = " << name << std::endl;

		std::cout << "function-parameters are: " << std::endl;
		jsize size = env->GetArrayLength(params);

		std::cout << "paramCount = " << size << std::endl;

		/*for (int i = 0; i < size; ++i) {
			jobject param = env->GetObjectArrayElement(params,i);

			//std::cout << "paramName = " << ug::vrl::jPointerGetName(env, param) << std::endl;

			    jclass argClass = env->GetObjectClass(param);

			    std::cout << "argClass = " << argClass << std::endl;

				jmethodID methodID = env->GetMethodID(argClass, "getName",
						"()Ljava/lang/String;");

				std::cout << "methodID = " << methodID << std::endl;

				std::string paramName = ug::vrl::stringJ2C(env,
						(jstring) env->CallObjectMethod(param, methodID));

				std::cout << "paramName = " << paramName << std::endl;
		}*/
	}

	const ug::bridge::ExportedFunction* func =
			ug::vrl::invocation::getFunctionBySignature(env,
					ug::vrl::vrlRegistry, name, params);

	if(func == NULL){
		std::cout << "func is null "<< std::endl;
	}else if (DEBUG) {
		std::cout << "UG__1invokeFunction(): func->m_name = " << func->name()
				<< std::endl;
	}

	ug::bridge::ParameterStack paramsIn;
	ug::bridge::ParameterStack paramsOut;

	jobject result = NULL;

	try {

		if (func == NULL) {
			std::stringstream ss;

			ss << "No function found that matches the given signature: "
					<< EMPHASIZE_BEGIN << name
					<< "(" + ug::vrl::getParamTypesAsString(env, params) + ")"
					<< EMPHASIZE_END << ".";

			jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
			env->ThrowNew(Exception, ss.str().c_str());

			return NULL;
		}

		//zum debuggen christian poliwoda ?? HIER passiert der Fehler
		 //UG_LOG("UG__1invokeFunction(): ?? HIER passiert der Fehler"<<std::endl);

		 //UG_LOG("UG__1invokeFunction(): paramsIn = "<< paramsIn<<std::endl);

		 for (int i = 0; i < paramsIn.size(); ++i) {
		 UG_LOG("UG__1invokeFunction(): paramsIn.class_name( "<<i<<" ) = "<< paramsIn.class_name(i)<<std::endl);
		 UG_LOG("UG__1invokeFunction(): paramsIn.class_name_node( "<<i<<" ) = "<< paramsIn.class_name_node(i)<<std::endl);
		 }

		 //UG_LOG(" func->params_in() = "<< func->params_in() <<std::endl);
		 for (int i = 0; i < func->params_in().size(); ++i) {
		 UG_LOG("UG__1invokeFunction(): func->params_in().class_name( "<<i<<" ) = "<< func->params_in().class_name(i)<<std::endl);
		 UG_LOG("UG__1invokeFunction(): func->params_in().class_name_node( "<<i<<" ) = "<< func->params_in().class_name_node(i)<<std::endl);
		 }

		ug::vrl::jobjectArray2ParamStack(env, ug::vrl::vrlRegistry, paramsIn,
				func->params_in(), params);

		func->execute(paramsIn, paramsOut);

		if (paramsOut.size() > 0) {
			//christian poliwoda
			//question: need here to distinct between Jobject and JobjectARRAY ?!??
			UG_LOG("UG__1invokeFunction(): result = ug::vrl::param2JObject("<<std::endl);

			result = ug::vrl::param2JObject(env, paramsOut, 0);
		}

	} catch (ug::bridge::UGError_ClassCastFailed& ex) {
		std::stringstream ss;
		ss << "Incompatible conversion in function " << func->name()
				<< "(): from " << ex.m_from << " to " << ex.m_to;

		jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
		env->ThrowNew(Exception, ss.str().c_str());

	} catch (ug::UGError& ex) {

		ug::vrl::throwUgErrorAsJavaException(env, ex);
	} catch (...) {
		std::stringstream ss;

		ss << "Unknown exception thrown while" << " trying to invoke function: "
				<< ug::vrl::stringJ2C(env, fName) << "().";

		jclass Exception = env->FindClass("edu/gcsc/vrl/ug/UGException");
		env->ThrowNew(Exception, ss.str().c_str());
	}

	return result;
}

JNIEXPORT jlong JNICALL Java_edu_gcsc_vrl_ug_UG__1getExportedClassPtrByName(
		JNIEnv *env, jobject obj, jstring name, jboolean classGrp) {

	if (ug::vrl::boolJ2C(classGrp)) {

		const ug::bridge::ClassGroupDesc* grpDesc =
				ug::vrl::vrlRegistry->get_class_group(
						ug::vrl::stringJ2C(env, name).c_str());

		if (grpDesc == NULL || grpDesc->get_default_class() == NULL) {
			return (jlong) NULL;
		}

		return (jlong) grpDesc->get_default_class();

	} else {
		return (jlong) ug::vrl::invocation::getExportedClassPtrByName(
				ug::vrl::vrlRegistry, ug::vrl::stringJ2C(env, name));
	}

	return (jlong) NULL;
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug_UG__1getDefaultClassNameFromGroup(
		JNIEnv *env, jobject obj, jstring grpName) {
	const ug::bridge::ClassGroupDesc* grpDesc =
			ug::vrl::vrlRegistry->get_class_group(
					ug::vrl::stringJ2C(env, grpName).c_str());

	if (grpDesc == NULL) {
		return ug::vrl::stringC2J(env, "");
	}

	if (grpDesc->get_default_class() == NULL) {
		return ug::vrl::stringC2J(env, "");
	}

	return ug::vrl::stringC2J(env, grpDesc->get_default_class()->name().c_str());
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug_UG__1getSvnRevision(JNIEnv *env,
		jobject obj) {
	std::string revision = ug::vrl::svnRevision();
	return ug::vrl::stringC2J(env, revision.c_str());
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug_UG__1getCompileDate(JNIEnv *env,
		jobject obj) {
	return ug::vrl::stringC2J(env, ug::UGCompileDate());
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug_UG__1getUGVersion(JNIEnv *env,
		jobject obj) {
	return ug::vrl::stringC2J(env, ug::UGGetVersionString().c_str());
}

JNIEXPORT void JNICALL Java_edu_gcsc_vrl_ug_UG__1delete(JNIEnv * env,
		jclass cls, jlong objPtr, jlong exportedClsPtr) {

	if (((void*) objPtr) != NULL && ((void*) exportedClsPtr) != NULL) {
		ug::bridge::IExportedClass* clazz =
				(ug::bridge::IExportedClass*) exportedClsPtr;
		clazz->destroy((void*) objPtr);
	}
}

JNIEXPORT void JNICALL Java_edu_gcsc_vrl_ug_UG__1invalidate(JNIEnv * env,
		jclass cls, jobject smartPtr) {

	if (ug::vrl::isJSmartPointerConst(env, smartPtr)) {
		ug::vrl::invalidateJConstSmartPointer(env, smartPtr);
	} else {
		ug::vrl::invalidateJSmartPointer(env, smartPtr);
	}

}

JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug_UG__1convertRegistryInfo(
		JNIEnv * env, jobject obj) {
	return ug::vrl::registry2NativeAPI(env, ug::vrl::vrlRegistry);
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug_UG__1getDescription(JNIEnv *env,
		jobject obj) {
	std::string desc = "UG is a general platform for the numerical solution<br>"
			" of partial differential equations.";

	return ug::vrl::stringC2J(env, desc.c_str());
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug_UG__1getAuthors(JNIEnv *env,
		jobject obj) {

	// START origin method DO NOT delete
	return ug::vrl::stringC2J(env, ug::UG_AUTHORS.c_str());
	// END origin method DO NOT delete

	/*
	 std::string tmpinfo = "project: VRL-UG \n  package: edu.gcsc.vrl.ug \n class: UG.java \n method: getAuthors() \n";
	 //std::string tmpresult = strncat( tmpinfo , " \n tmpresult :-)", tmpinfo.length());

	 std::string tmpresult = " ### ### ### tmpresult :-)";

	 std::cout << tmpresult << std::endl;

	 return ug::vrl::stringC2J(env, tmpresult.c_str());
	 */
}

JNIEXPORT jstring JNICALL Java_edu_gcsc_vrl_ug_UG__1getBinaryLicense(
		JNIEnv *env, jobject obj) {
	return ug::vrl::stringC2J(env, ug::UG_BINARY_LICENSE.c_str());
}

//JNIEXPORT void JNICALL Java_edu_gcsc_vrl_ug_UG_attachCanvas
//(JNIEnv *env, jobject obj, jobject canvas) {
//	ug::vrl::Canvas::getInstance()->setJObject(env, canvas);
//
//	//	ug::vrl::Canvas::getInstance()->addObject(ug::vrl::string2JObject(env,"Test_String"));
//}

//method needed only for debug
//after debugging is finished this method should be commented out
jobject bool_array(JNIEnv *env, jobject obj, jstring jName,
		jobjectArray params) {

	std::cout << "trunk/ugbase/bindings/vrl/bindings_vrl.cpp : bool_array() "
			<< std::endl;

	std::cout << " jName = " << ug::vrl::stringJ2C(env, jName) << std::endl;

	/*
	 jobject firstElement = env->GetObjectArrayElement(params, 0);

	 jobjectArray nestedObjectArray = (jobjectArray) firstElement;

	 jbooleanArray boolArray = (jbooleanArray) nestedObjectArray;

	 jsize size = env->GetArrayLength(boolArray);
	 jboolean *elements = env->GetBooleanArrayElements(boolArray, false);

	 std::cout << " true = " << true << std::endl;
	 std::cout << " false = " << false << std::endl;

	 //do something with the entries
	 for (int i = 0; i < size; i++) {

	 std::cout << " elements[ " << i << " ] = " << elements[i] << std::endl;
	 }

	 env->ReleaseBooleanArrayElements(boolArray, elements, 0);

	 return firstElement;
	 */

	//
	//
	//
	//
	jbooleanArray jBoolArray = ug::vrl::jObject2BooleanArray(env, obj);

	std::cout << "CPP: jBoolArray: " << jBoolArray << std::endl;

	jsize arrayLenght = env->GetArrayLength(jBoolArray);
	std::cout << "arrayLenght = " << arrayLenght << std::endl;

	jboolean *arrayBoolElements = env->GetBooleanArrayElements(jBoolArray,
			NULL);

	// getVALUEs from jBooleanArray

	for (int arrayindex = 0; arrayindex < arrayLenght; ++arrayindex) {

		////
		//// ATTENTION !!!
		////
		////we need to cast manuel between jboolean and bool
		std::cout << arrayindex << " (bool) arrayBoolElements = "
				<< (bool) arrayBoolElements[arrayindex] << std::endl;
	}

	env->ReleaseBooleanArrayElements(jBoolArray, arrayBoolElements, JNI_ABORT);

	return jBoolArray;

}

//method needed only for debug
//after debugging is finished this method should be commented out
jobject array_of_bool_arrays(JNIEnv *env, jobject obj, jstring jName,
		jobjectArray params) {
 jobject result = NULL;

	std::cout
			<< "trunk/ugbase/bindings/vrl/bindings_vrl.cpp : array_of_bool_arrays() "
			<< std::endl;

	std::cout << " array_of_bool_arrays().jName = "
			<< ug::vrl::stringJ2C(env, jName) << std::endl;

	std::cout
			<< " array_of_bool_arrays() calling  ug::vrl::jObject2BooleanArray(env, PARAMS)"
			<< std::endl;

	//
	// here i need to know too that there are arrays in "params"
	//

	//RICHTIG ?? START

	// for (array.size)
	jsize arraySize = env->GetArrayLength(params);

	std::cout << "for (int i = 0; i < " << arraySize << "; ++i)" << std::endl;
	for (int i = 0; i < arraySize; ++i) {
		//std::cout << "Array i = " << i << ": "<< std::endl;

		std::cout << "jobject value = env->GetObjectArrayElement(array, " << i
				<< " );" << std::endl;
		jobject value = env->GetObjectArrayElement(params, i);

		//check if the parameter is an array

		bool isArray = ug::vrl::isJObjectAnArray(env, value);

		std::cout << "value " << i << " is an array: "
				<< (isArray ? "TRUE" : "FALSE") << std::endl;

		//here we KNOW that the jobject "value" contains also an array!!! therefor calling
		std::cout
				<< "here we KNOW that the jobject value contains also an array!!! therefor calling"
				<< std::endl;

		//std::cout <<"Java_de_tutorials_NativeExample_nativePrintArray(env, clazz, value);"<< std::endl;
		////Java_de_tutorials_NativeExample_nativePrintArray(env, clazz, value);

		//ersetzte den aufruf von bsp _nativePrintArray() durch den code/inhalt der methode

		jbooleanArray jBoolArray = ug::vrl::jObject2BooleanArray(env, value);

		// SETTING RESULT
		result = jBoolArray;

		jsize jBoolArrayLenght = env->GetArrayLength(jBoolArray);
		//cout << "jBoolArrayLenght = " << jBoolArrayLenght << endl;

		jboolean *jBoolArrayElements = env->GetBooleanArrayElements(jBoolArray,
				NULL);

		std::cout << "trunk/ugbase/bindings/vrl/bindings_vrl.cpp : array_of_bool_arrays() ."
				<< " getVALUEs from jBooleanArray" << std::endl;

		for (int arrayindex = 0; arrayindex < jBoolArrayLenght; ++arrayindex) {

			////
			//// ATTENTION !!!
			////
			////we need to cast manuel between jboolean and bool
			std::cout << " (bool) jBoolArrayElements[ " << arrayindex << " ] = "
					<< (bool) jBoolArrayElements[arrayindex] << std::endl;
		}

		env->ReleaseBooleanArrayElements(jBoolArray, jBoolArrayElements,
				JNI_ABORT);

	}					//for (int i = 0; i < arraySize; ++i)
	//RICHTIG ?? END

	return result;
}

/*
 * this method is only for debug issues and should NOT be used later
 */JNIEXPORT jobject JNICALL Java_edu_gcsc_vrl_ug_UG__1test_1debug(JNIEnv *env,
		jobject obj, jstring jName, jobjectArray params) {

	std::cout << "CPP: trunk/ugbase/bindings/vrl/bindings_vrl.cpp :"
			<< " Java_edu_gcsc_vrl_ug_UG__1test_1debug() " << std::endl;
	std::cout << " Java_edu_gcsc_vrl_ug_UG__1test_1debug() : jName = " << jName
			<< std::endl;

	/* By selecting the debug-method you choose the test that should
	 * by trigert.
	 */

	/* std::cout << " bool_array()" << std::endl;
	 jobject result = bool_array(env, obj, jName, params);
	 */

	jobject result;

	std::cout << " UG__1test_1debug() bottle neck."
			<< " redirecting method call to array_of_bool_arrays()"
			<< std::endl;
	result = array_of_bool_arrays(env, obj, jName, params); //DID WORK

	std::cout << " result = " << result << std::endl; //what is result array_of_bool_arrays does not return anything
	//std::cout << " jObject2String(result) = " << ug::vrl::jObject2String(env, result) << std::endl;

	//
	//NOW TRY CALLING "in a real scenario"
	// result = Java_edu_gcsc_vrl_ug_UG__1invokeFunction(env, obj, jName, true, params);

	return result;
}
