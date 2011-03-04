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

#include "ug_script/user_data/user_data.h"

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


class TestClass {
public:

	TestClass() {
		//
	}

	int performTest() {

		std::cout << "***0\n";

		const ug::bridge::IExportedClass* cls =
				ug::vrl::invocation::getExportedClassPtrByName(
				ug::vrl::vrlRegistry, "Domain2d");

		std::cout << "***1\n";

		void* obj = cls->create();

		UG_LOG("Domain2d: " << (long) obj << ", " << obj << "\n");

		std::cout << "***2\n";

		int methodID = -1;

		for (unsigned int i = 0; i < cls->num_methods(); i++) {
			std::cout << "***3:" << i << std::endl;
			if (cls->get_method(i).name() == "get_grid") {
				UG_LOG("method found!\n");
				std::cout << "***3:found, ID=" << i << "\n";
				methodID = i;
				break;
			}
		}

		std::cout << "***4\n";

		const ug::bridge::ExportedMethod& method = cls->get_method(methodID);

		UG_LOG("Call Method:" << method.name() << "\n");

		ug::bridge::ParameterStack paramsIn;
		ug::bridge::ParameterStack paramsOut;

		std::cout << "***5\n";

		method.execute(obj, paramsIn, paramsOut);

		if (paramsOut.size() > 0
				&& paramsOut.get_type(0) == ug::bridge::PT_POINTER) {
			UG_LOG("Output: " << (long) paramsOut.to_pointer(0) << ", " <<
					paramsOut.to_pointer(0) << "\n");
		}

		std::cout << "***6\n";

		return 23; // :D
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

	SmartPtr<TestClass> smartTestImpl() {
		return SmartPtr<TestClass > (new TestClass());
	}

	ConstSmartPtr<TestClass> constSmartTestImpl() {
		return ConstSmartPtr<TestClass > (new TestClass());
	}

	int print_name() {
		UG_LOG("Name is Test\n");
		return 1;
	}

	int print() {
		UG_LOG("Test::print()\n");
		return 0;
	}

	int print2() {
		UG_LOG("Test::print2()\n");
		return 1;
	}

	int print2() const {
		UG_LOG("Test::print2() const\n");
		return 1;
	}

	~TestClass() {
		UG_LOG("~TestClass" << std::endl);
	}
};

int SmartTestFunction(SmartPtr<TestClass> test) {
	UG_LOG("SmartTestFunc: ");

	test->print2();

	return test.get_refcount();
}

int ConstSmartTestFunction(ConstSmartPtr<TestClass> test) {
	UG_LOG("ConstSmartTestFunc: ");

	test->print2();

	return test.get_refcount();
}

template <typename TVector>
TVector* CreateVector(int dim) {
	return new TVector(dim);
}

template <typename TVector>
SmartPtr<TVector> CreateVectorSmart(int dim) {
	return SmartPtr<TVector > (new TVector(dim));
}




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

	reg.add_class_<TestClass > ("TestClass", "testing")
			.add_constructor()
			.add_method("svnRevision|hide=true,interactive=false", &TestClass::getRev)
			.add_method("add", &TestClass::add, "result",
			"a|default|min=-3;max=5;value=-12#b|default|min=-1;max=1;value=23")
			.add_method("getString", &TestClass::getString)
			.add_method("performTest", &TestClass::performTest)
			.add_method("print", &TestClass::print)
			.add_method("smartTestImpl", &TestClass::smartTestImpl)
			.add_method("constSmartTestImpl", &TestClass::constSmartTestImpl);

	reg.add_function("SmartTestFunction", &SmartTestFunction, "testing");
	reg.add_function("ConstSmartTestFunction", &ConstSmartTestFunction, "testing");

	ug::vrl::RegisterVRLUserNumber(reg, "testing");

	/************************************/

	reg.add_class_<Vector<double> >("DVector", "testing")
			.add_constructor()
			.add_method("set|hide=true", (bool (Vector<double> ::*)(number)) & Vector<double> ::set,
			"Success", "Number")
			.add_method("size|hide=true", (size_t(Vector<double> ::*)()) & Vector<double> ::size,
			"Size", "")
			.add_method("set_random|hide=true", (bool (Vector<double> ::*)(number)) & Vector<double> ::set_random,
			"Success", "Number")
			.add_method("print|hide=true", &Vector<double> ::p);

	reg.add_function("VecScaleAdd2", &VecScaleAdd2<Vector<double> >, "",
			"dest, alpha1, vec1, alpha2, vec2", "dest = alpha1*vec1 + alpha2*vec2");
	reg.add_function("VecScaleAdd3", &VecScaleAdd3<Vector<double> >, "",
			"dest, alpha1, vec1, alpha2, vec2, alpha3, vec3", "dest = alpha1*vec1 + alpha2*vec2 + alpha3*vec3");

	reg.add_function("CreateVector",
			&CreateVector<Vector<double> >, "testing", "DVector", "Dimension");

	/************************************/

	// Register algebra
	CPUAlgebraChooser chooser;
	ug::bridge::RegisterDynamicLibAlgebraInterface(
			reg, chooser.get_algebra_type());
	ug::bridge::RegisterDynamicLibDiscretizationInterface(
			reg, chooser.get_algebra_type());

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
	return ug::vrl::registry2NativeAPI(env, ug::vrl::vrlRegistry);
}



//JNIEXPORT void JNICALL Java_edu_gcsc_vrl_ug4_UG4_attachCanvas
//(JNIEnv *env, jobject obj, jobject canvas) {
//	ug::vrl::Canvas::getInstance()->setJObject(env, canvas);
//
//	//	ug::vrl::Canvas::getInstance()->addObject(ug::vrl::string2JObject(env,"Test_String"));
//}
