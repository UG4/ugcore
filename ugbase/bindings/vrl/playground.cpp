/*
 * File:   Playground.cpp
 * Author: Michael Hoffer <info@michaelhoffer.de>
 *
 * Created on 12. April 2011, 12:49
 */

#include "playground.h"
#include "bindings_vrl.h"



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

	return test.refcount();
}

int ConstSmartTestFunction(ConstSmartPtr<TestClass> test) {
	UG_LOG("ConstSmartTestFunc: ");

	test->print2();

	return test.refcount();
}

void TestSmartPtr2ConstPtr(const TestClass* t) {
	t->print2();
}

SmartPtrCls::SmartPtrCls() {
	UG_LOG("Constructor SmartPtrCls() called." << std::endl);
	_name = "noname";
	_data = NULL;
}

SmartPtrCls::SmartPtrCls(std::string name) {
	_name = name;
	UG_LOG("Constructor SmartPtrCls(std::string name) called." << std::endl);
	_data = NULL;
}

void SmartPtrCls::print_name() {
	UG_LOG("Name is " << _name << std::endl);
}

void SmartPtrCls::create_data(int size) {
	if (_data ==NULL) {
		_data = new byte[size];
	} else {
		UG_LOG(_name << ": Data already created! " << std::endl);
	}
}

SmartPtrCls::~SmartPtrCls() {
	UG_LOG("~SmartPtrCls:" << _name << std::endl);

	if (_data !=NULL ){
		delete[] _data;
		_data = NULL;
	}
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
	reg.add_function("TestSmartPtr2ConstPtr", &TestSmartPtr2ConstPtr, "ug4/testing");

	reg.add_class_<SmartPtrCls > ("SmartPtrCls", "ug4/testing")
			.add_constructor()
			.add_constructor<void(*)(std::string)>()
			.add_method("print_name", &SmartPtrCls::print_name)
			.add_method("create_data", &SmartPtrCls::create_data)
			.set_construct_as_smart_pointer(true);
}

} // end vrl::
}// end ug::


