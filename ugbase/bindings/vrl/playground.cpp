///*
// * File:   Playground.cpp
// * Author: Michael Hoffer <info@michaelhoffer.de>
// *
// * Created on 12. April 2011, 12:49
// */
//
//#include "playground.h"
//
//#include "compiledate.h"
//#include "bindings_vrl.h"
//#include "vrl_user_number.h"
//
//
//
//namespace ug {
//namespace vrl {
//
//class A {
//public:
//	std::string name;
//
//	A() {
//		name = "A";
//	}
//
//	void say() {
//		UG_LOG("Name:" << A::name << std::endl);
//	}
//};
//
//class B : public A {
//public:
//	static std::string name;
//
//	B()
//	{
//		name = "B";
//	}
//	void say() {
//		UG_LOG("Name:" << B::name << std::endl);
//	}
//};
//
//class C {
//public:
//	static std::string name;
//
//	C() {
//		name = "C";
//	}
//
//	void say() {
//		UG_LOG("Name:" << C::name << std::endl);
//	}
//};
//
//
//class D : public B, public C {
//public:
//	std::string name;
//
//	D() {
//		name = "D";
//	}
//
//	void say() {
//		UG_LOG("Name:" << D::name << std::endl);
//	}
//};
//
//
//
////**************** PLAYGROUND CLASSES
//
//class TestClass {
//public:
//
//	TestClass() {
//		//
//	}
//
//	//	int performTest() {
//	//
//	//		std::cout << "***0\n";
//	//
//	//		const ug::bridge::IExportedClass* cls =
//	//				ug::vrl::invocation::getExportedClassPtrByName(
//	//				ug::vrl::vrlRegistry, "Domain2d");
//	//
//	//		std::cout << "***1\n";
//	//
//	//		void* obj = cls->create();
//	//
//	//		UG_LOG("Domain2d: " << (long) obj << ", " << obj << "\n");
//	//
//	//		std::cout << "***2\n";
//	//
//	//		int methodID = -1;
//	//
//	//		for (unsigned int i = 0; i < cls->num_metho    } // end vrl::
//	//} // end ug::ds(); i++) {
//	//			std::cout << "***3:" << i << std::endl;
//	//			if (cls->get_method(i).name() == "get_grid") {
//	//				UG_LOG("method found!\n");
//	//				std::cout << "***3:found, ID=" << i << "\n";
//	//				methodID = i;
//	//				break;
//	//			}
//	//		}
//	//
//	//		std::cout << "***4\n";
//	//
//	//		const ug::bridge::ExportedMethod& method = cls->get_method(methodID);
//	//
//	//		UG_LOG("Call Method:" << method.name() << "\n");
//	//
//	//		ug::bridge::ParameterStack paramsIn;
//	//		ug::bridge::ParameterStack paramsOut;
//	//
//	//		std::cout << "***5\n";
//	//
//	//		method.execute(obj, paramsIn, paramsOut);
//	//
//	//		if (paramsOut.size() > 0
//	//				&& paramsOut.get_type(0) == ug::bridge::PT_POINTER) {
//	//			UG_LOG("Output: " << (long) paramsOut.to_pointer(0) << ", " <<
//	//					paramsOut.to_pointer(0) << "\n");
//	//		}
//	//
//	//		std::cout << "***6\n";
//	//
//	//		return 23; // :D
//	//	}
//
//	std::string getRev() {
//		return ug::vrl::svnRevision();
//	}
//
//	int add(int a, int b) {
//		return a + b;
//	}
//
//	std::string getString() {
//		UG_LOG("Test123" << std::endl);
//		return "Test123";
//	}
//
//	SmartPtr<TestClass> smartTestImpl() {
//		return SmartPtr<TestClass > (new TestClass());
//	}
//
//	ConstSmartPtr<TestClass> constSmartTestImpl() {
//		return ConstSmartPtr<TestClass > (new TestClass());
//	}
//
//	int print_name() {
//		UG_LOG("Name is Test\n");
//		return 1;
//	}
//
//	int print() {
//		UG_LOG("Test::print()\n");
//		return 0;
//	}
//
//	int print2() {
//		UG_LOG("Test::print2()\n");
//		return 1;
//	}
//
//	int print2() const {
//		UG_LOG("Test::print2() const\n");
//		return 1;
//	}
//
//	~TestClass() {
//		UG_LOG("~TestClass" << std::endl);
//	}
//};
//
//int SmartTestFunction(SmartPtr<TestClass> test) {
//	UG_LOG("SmartTestFunc: ");
//
//	test->print2();
//
//	return test.get_refcount();
//}
//
//int ConstSmartTestFunction(ConstSmartPtr<TestClass> test) {
//	UG_LOG("ConstSmartTestFunc: ");
//
//	test->print2();
//
//	return test.get_refcount();
//}
//
//template <typename TVector>
//TVector* CreateVector(int dim) {
//	return new TVector(dim);
//}
//
//template <typename TVector>
//SmartPtr<TVector> CreateVectorSmart(int dim) {
//	return SmartPtr<TVector > (new TVector(dim));
//}
//
//void registerPlayground(ug::bridge::Registry& reg) {
//	reg.add_class_<TestClass > ("TestClass", "ug4/testing")
//			.add_constructor()
//			.add_method("svnRevision|hide=true,interactive=true", &TestClass::getRev)
//			.add_method("add", &TestClass::add, "result",
//			"a|default|min=-3;max=5;value=-12#b|default|min=-1;max=1;value=23")
//			.add_method("getString", &TestClass::getString)
//			//			.add_method("performTest", &TestClass::performTest)
//			.add_method("print", &TestClass::print)
//			.add_method("smartTestImpl", &TestClass::smartTestImpl)
//			.add_method("constSmartTestImpl", &TestClass::constSmartTestImpl);
//
//	reg.add_function("SmartTestFunction", &SmartTestFunction, "ug4/testing");
//	reg.add_function("ConstSmartTestFunction", &ConstSmartTestFunction, "ug4/testing");
//
//	reg.add_class_<A > ("A", "ug4/testing").
//			add_constructor().
//			add_method("say", &A::say);
//
//	reg.add_class_<B, A > ("B", "ug4/testing").
//			add_constructor().
//			add_method("say", &B::say);
//
//	reg.add_class_<C > ("C", "ug4/testing").
//			add_constructor().
//			add_method("say", &C::say);
//
//	reg.add_class_<D, B, C > ("D", "ug4/testing").
//			add_constructor().
//			add_method("say", &D::say);
//	
//}
//
//} // end vrl::
//}// end ug::
//
//
