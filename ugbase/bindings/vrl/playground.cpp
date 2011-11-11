/*
 * File:   Playground.cpp
 * Author: Michael Hoffer <info@michaelhoffer.de>
 *
 * Created on 12. April 2011, 12:49
 */

#include "playground.h"

#include "compiledate.h"
#include "build_hostname.h"
#include "bindings_vrl.h"
#include "vrl_user_number.h"



namespace ug {
namespace vrl {


TestClass::TestClass() {
	UG_LOG("Constructor TestClass() called." << std::endl);
}

TestClass::TestClass(std::string name) {
	UG_LOG("Constructor TestClass(std::string name) called." << std::endl);
}

std::string TestClass::getRev() {
	return ug::vrl::svnRevision();
}

int TestClass::add(int a, int b) {
	return a + b;
}

std::string TestClass::getString() {
	UG_LOG("Test123" << std::endl);
	return "Test123";
}

SmartPtr<TestClass> TestClass::smartTestImpl() {
	return SmartPtr<TestClass > (new TestClass());
}

ConstSmartPtr<TestClass> TestClass::constSmartTestImpl() {
	return ConstSmartPtr<TestClass > (new TestClass());
}

int TestClass::print_name() {
	UG_LOG("Name is Test\n");
	return 1;
}

int TestClass::print() {
	UG_LOG("Test::print()\n");
	return 0;
}

int TestClass::print2() {
	UG_LOG("Test::print2()\n");
	return 1;
}

int TestClass::print2() const {
	UG_LOG("Test::print2() const\n");
	return 1;
}

TestClass::~TestClass() {
	UG_LOG("~TestClass" << std::endl);
}

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

void registerPlayground(ug::bridge::Registry& reg) {
	reg.add_class_<TestClass > ("TestClass", "ug4/testing")
			.add_constructor()
			.add_constructor<void(*)(std::string)>()
			.add_method("svnRevision|hide=true,interactive=true", &TestClass::getRev)
			.add_method("add", &TestClass::add, "result",
			"a|default|min=-3;max=5;value=-12#b|default|min=-1;max=1;value=23")
			.add_method("getString", &TestClass::getString)
			//			.add_method("performTest", &TestClass::performTest)
			.add_method("print", &TestClass::print)
			.add_method("smartTestImpl", &TestClass::smartTestImpl)
			.add_method("constSmartTestImpl", &TestClass::constSmartTestImpl);

	reg.add_function("SmartTestFunction", &SmartTestFunction, "ug4/testing");
	reg.add_function("ConstSmartTestFunction", &ConstSmartTestFunction, "ug4/testing");
}

} // end vrl::
}// end ug::


