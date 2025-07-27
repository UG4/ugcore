/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter, Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <iostream>
#include <sstream>
#include <boost/bind/bind.hpp>
#include "registry/registry.h"
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "common/util/message_hub.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_function_handle.h"
#endif

//	temporary include
#include "lib_grid/attachments/page_container.h"

using namespace std;

namespace ug
{
namespace bridge
{

/// \defgroup test_bridge Test Bridge
/// \ingroup misc_bridge
/// \{

// prints "Hello World"
void PrintHelloWorldToScreen()
{
	UG_LOG("Hello World !\n");
}

int Add(int a, int b)
{
	return a + b;
}

int Add(int a, int b, int c)
{
	return a + b + c;
}

string Add(const char* a, const char* b)
{
	string str = a;
	str.append(b);
	return str;
}

class Test
{
	public:
		Test() { UG_LOG("Test::Test() constructor used.\n")}
		Test(const char* msg)
		{
			UG_LOG("Test::Test(const char*) constructor used.\n")
			UG_LOG("Message is: '"<<msg<<"'.\n");
		}

		int add(int a, int b)
		{
			return a+b;
		}

		int add(int a, int b, int c){
			return a + b + c;
		}

		string add(const char* a, const char* b){
			string str = a;
			str.append(b);
			return str;
		}

		int print_name()
		{
			UG_LOG("Name is Test\n");
			return 1;
		}
		
		int print()			{UG_LOG("Test::print()\n"); return 0;}
		int print() const	{UG_LOG("Test::print() const\n"); return 1;}
};


///	calls non-const print on non-const object of class Test
int TestFunc(Test& t)
{
	return t.print();
}

///	calls const print on const object of class Test
int ConstTestFunc(const Test& t)
{
	return t.print();
}

///	returns a const reference to a non-const object of class Test
const Test& ToConst(Test& test)
{
	return test;
}

SmartPtr<Test> SmartTestImpl(){
	return SmartPtr<Test>(new Test());
}

ConstSmartPtr<Test> ConstSmartTestImpl(){
	return ConstSmartPtr<Test>(new Test());
}

int SmartTestFunc(SmartPtr<Test> test)
{
	UG_LOG("SmartTestFunc: ");
	return test->print();
}

int ConstSmartTestFunc(ConstSmartPtr<Test> test)
{
	UG_LOG("ConstSmartTestFunc: ");
	return test->print();
}

///	SmartTest is a test-class which shall only be used encapsulated in a smart-pointer
class SmartTest{
	public:
		SmartTest()	{}
		virtual ~SmartTest()		{UG_LOG("SmartTest destroyed...\n");}
		virtual void say_hello()	{UG_LOG("Hello, I'm SmartTest\n");}
};

class SmartTestDerived : public SmartTest
{
	public:
		SmartTestDerived()	{}
		virtual ~SmartTestDerived()	{UG_LOG("SmartTestDerived destroyed...\n");}
		virtual void say_hello()	{UG_LOG("Hello, I'm SmartTestDerived\n");}
};

typedef SmartPtr<SmartTest> SPSmartTest;
typedef SmartPtr<SmartTestDerived> SPSmartTestDerived;

void SmartTestArrived(SPSmartTest test)
{
	UG_LOG("HE, SMART TEST ARRIVED. HI SMART TEST.\n");
	UG_LOG("smart test answers: ");
	test->say_hello();
}



class Piece
{
	public:
		Piece()	: m_size(0)	{}
		Piece(int size) : m_size(size)	{}

		int size()
		{
			return m_size;
		}

	protected:
		int m_size;
};

class Cake
{
	public:
		Cake() : m_numPieces(16)	{}

		Piece& take_pieces(int size)
		{
			if(size > m_numPieces)
				size = m_numPieces;
			m_numPieces -= size;
			return *(new Piece(size));
		}

		int add_pieces(Piece& piece)
		{
			m_numPieces += piece.size();
			return m_numPieces;
		}

		int pieces_left()	{return m_numPieces;}

	protected:
		int m_numPieces;
};

// Some Base class
class Base
{
	public:
		virtual ~Base()	{}
		virtual void print() const
		{
			UG_LOG("Base::print() called\n");
		}
};

// Some Derived class
class Derived : public Base
{
	public:
		virtual ~Derived()	{}
		virtual void print() const
		{
			UG_LOG("Derived::print() called.\n");
		}
};

// Some Derived class
class FurtherDerived : public Derived
{
	public:
		virtual ~FurtherDerived()	{}
		virtual void print() const
		{
			UG_LOG("FurtherDerived::print() called.\n");
		}
};

const FurtherDerived* CreateConstFurtherDerived()
{
	return new FurtherDerived();
}

class Base0
{
	public:
		virtual ~Base0() {}
		virtual void virt_print_base0() const = 0;
		void print_base0() const {UG_LOG("Base0::print_base0() called.\n");}
};

class Base1
{
	public:
		virtual ~Base1() {}
		virtual void virt_print_base1() const = 0;
		void print_base1() const {UG_LOG("Base1::print_base1() called.\n");}
};

class Base2
{
	public:
		virtual ~Base2() {}
		virtual void virt_print_base2() const = 0;
		void print_base2() const {UG_LOG("Base2::print_base2() called.\n");}
};

class Base3
{
	public:
		virtual ~Base3() {}
		virtual void virt_print_base3() const = 0;
		void print_base3() const {UG_LOG("Base3::print_base3() called.\n");}
};

class Intermediate0 : public Base0, public Base1
{
	public:
		virtual ~Intermediate0() {}
		virtual void virt_print_intermediate0() const = 0;
		void print_intermediate0() const {UG_LOG("Intermediate0::print_intermediate0() called.\n");}
};

class Intermediate1 : public Base2, public Base3
{
	public:
		virtual ~Intermediate1() {}
		virtual void virt_print_intermediate1() const = 0;
		void print_intermediate1() const {UG_LOG("Intermediate1::print_intermediate1() called.\n");}
};

class MultipleDerived : public Intermediate0, public Intermediate1
{
	public:
		virtual ~MultipleDerived() {}
		void print_mulitple_derived(){UG_LOG("MultipleDerived::print() called\n");}

		virtual void virt_print_intermediate0() const	{UG_LOG("MultipleDerived::virt_print_intermediate0() called\n");}
		virtual void virt_print_intermediate1() const	{UG_LOG("MultipleDerived::virt_print_intermediate1() called\n");}

		virtual void virt_print_base0() const	{UG_LOG("MultipleDerived::virt_print_base0() called\n");}
		virtual void virt_print_base1() const	{UG_LOG("MultipleDerived::virt_print_base1() called\n");}
		virtual void virt_print_base2() const	{UG_LOG("MultipleDerived::virt_print_base2() called\n");}
		virtual void virt_print_base3() const	{UG_LOG("MultipleDerived::virt_print_base3() called\n");}
};

SmartPtr<MultipleDerived> SmartMultipleDerivedImpl(){
	return SmartPtr<MultipleDerived>(new MultipleDerived());
}

void PrintFunction(SmartPtr<Base3> b)
{
	b->virt_print_base3();
	b->print_base3();
}

void PrintFunction(SmartPtr<Base2> b)
{
	b->virt_print_base2();
	b->print_base2();
}

void PrintFunctionIntermediate(SmartPtr<Intermediate1> b)
{
	b->virt_print_intermediate1();
	b->print_intermediate1();
}

class ConstClass
{
	public:
		string const_method() const
		{
			return "const_method_called";
		}
};

void PrintFunctionIntermediate(Intermediate0& b)
{
	b.virt_print_intermediate0();
	b.print_intermediate0();
}

void PrintFunction(Base3& b)
{
	b.virt_print_base3();
	b.print_base3();
}

void PrintFunction(Base& b)
{
	b.print();
}

const char* StringTest()
{
	return "Jooik";
}

string StdStringTest()
{
	return string("stdJooik");
}

void PrintStringTest(const std::string& str)
{
	UG_LOG("PrintStringTest recieved: " << str << endl);
}


void TestPageContainer()
{
	UG_LOG("Testing page container whith default maxPageSize.\n");

	typedef PageContainer<double> PageCon;
	PageCon* ppc = new PageCon;
	PageCon& pc = *ppc;

	size_t testSize = 10000;

	UG_LOG("resizing...\n");
	pc.resize(testSize);

	UG_LOG("pushing doubles...\n");
	for(size_t i = 0; i < testSize; ++i){
		pc[i] = 3300. + (double)i;
	}

	UG_LOG("PageContainer content:");
	for(size_t i = testSize - 10; i < testSize; ++i){
		UG_LOG(" " << pc[i]);
	}
	UG_LOG(endl);
	UG_LOG("Releasing PageContainer...\n");
	delete ppc;
	UG_LOG("done.\n");
}

void PostRegisteredFunction()
{
	UG_LOG("PostRegisteredFunction successfully executed.\n");
}

void PostRegisterTest()
{
	bridge::Registry& reg = bridge::GetUGRegistry();
	reg.add_function("PostRegisteredFunction", &PostRegisteredFunction);
	reg.registry_changed();
}


class Unregistered{
	private:
		string	m_emptyString;
};

Unregistered* UnregisteredParameterTest(Unregistered& unregistered)
{
	return &unregistered;
}

class NonAllowedName1 {};
class NonAllowedName2 {};
class NonAllowedName3 {};
class NonAllowedName4 {};
class AllowName{public: void non_allowed_method(){};};

void NonAllowedFct1() {};
void NonAllowedFct2() {};
void NonAllowedFct3() {};


const char* REFINEMENT_MESSAGE_NAME = "MSG_GRID_ADAPTION";

class Message{
	//int msgType;
};



////////////////////////////////////////////////////////////////////////////////
//	The following test-case is used to check the MessageHub class
class TestMessage : public MessageHub::IMessage
{
	public:
		string	m_strMsg;
};

void TestMsgCallback(const TestMessage& msg){
	UG_LOG("Callback func received message: " << msg.m_strMsg << "\n");
}


class MessageHubTest
{
	public:
		MessageHubTest(){
			m_msgHub = SPMessageHub(new MessageHub());
			m_callbackId = m_msgHub->register_class_callback(this,
										  &ug::bridge::MessageHubTest::callback);

			m_msgHub->register_function_callback(&ug::bridge::TestMsgCallback);
		}

		void callback(const TestMessage& msg){
			UG_LOG("Received message: " << msg.m_strMsg << "\n");
		}

		void post_message(const char* message){
			if(m_msgHub.valid()){
				TestMessage msg;
				msg.m_strMsg = message;
				m_msgHub->post_message(msg);
			}
		}

	protected:
		SPMessageHub	m_msgHub;
		MessageHub::SPCallbackId m_callbackId;
};


template <typename T>
void PrintStdVector(std::vector<T> vec)
{
	UG_LOG("std::vector<T> passed: Size is: "<<vec.size()<<"\n");
	for(size_t i = 0; i < vec.size(); ++i)
		UG_LOG(i<<": "<<vec[i]<<"\n");
}

template <typename T>
void PrintConstStdVectorRef(const std::vector<T>& vec)
{
	UG_LOG("const std::vector<T>& passed: Size is: "<<vec.size()<<"\n");
	for(size_t i = 0; i < vec.size(); ++i)
		UG_LOG(i<<": "<<vec[i]<<"\n");
}

template <typename T>
void PrintStdVectorOfClass(std::vector<T> vec)
{
	UG_LOG("std::vector<T> passed: Size is: "<<vec.size()<<"\n");
	for(size_t i = 0; i < vec.size(); ++i)
		vec[i]->print();
}

template <typename T>
void ConstStdVectorRefOfClass(const std::vector<T>& vec)
{
	UG_LOG("const std::vector<T>& passed: Size is: "<<vec.size()<<"\n");
	for(size_t i = 0; i < vec.size(); ++i)
		vec[i]->print();
}

std::vector<std::string> ReturnStdVector_String(){
	std::vector<std::string> vec;
	vec.push_back("Entry1");
	vec.push_back("Entry2");
	vec.push_back("Entry3");
	return vec;
}

const std::vector<std::string>& ReturnConstStdVectorRef_String(){
	static std::vector<std::string> vec;
	vec.push_back("Entry1");
	vec.push_back("Entry2");
	vec.push_back("Entry3");
	return vec;
}

template <typename T>
std::vector<T> ReturnStdVector_Number(){
	std::vector<T> vec;
	vec.push_back((T)1.56);
	vec.push_back((T)144.87);
	vec.push_back((T)99.3);
	return vec;
}

template <typename T>
const std::vector<T>& ReturnConstStdVectorRef_Number(){
	static std::vector<T> vec;
	vec.push_back((T)1.56);
	vec.push_back((T)144.87);
	vec.push_back((T)99.3);
	return vec;
}

template <typename T>
std::vector<T> ReturnStdVectorOfClass(){
	std::vector<T> vec;
	vec.push_back(T(new Derived()));
	vec.push_back(T(new Base()));
	vec.push_back(T(new FurtherDerived()));
	return vec;
}

template <typename T>
const std::vector<T>& ReturnConstStdVectorRefOfClass(){
	static std::vector<T> vec;
	vec.push_back(new Derived());
	vec.push_back(new Base());
	vec.push_back(new FurtherDerived());
	return vec;
}

template <typename T>
const std::vector<T>& ReturnConstStdVectorRefOfClassSP(){
	static std::vector<T> vec;
	vec.push_back(make_sp(new Derived()));
	vec.push_back(make_sp(new Base()));
	vec.push_back(make_sp(new FurtherDerived()));
	return vec;
}

void NotAllowedParamPerValue(Piece P){}

#ifdef UG_FOR_LUA
void FunctionWithFunctionArgument(LuaFunctionHandle handle)
{
	UG_LOG("Test LuaFunctionHandle: handle.ref=" << handle.ref << std::endl);
	return;
}
#endif

void RegisterBridge_Test(Registry& reg, string parentGroup)
{
//	get group string
	stringstream groupString; groupString << parentGroup << "/Test";
	string grp = groupString.str();

	try
	{
	//	registering hello world
		reg.add_function("PrintHelloWorld", &PrintHelloWorldToScreen, grp);

#ifdef UG_FOR_LUA
	//	registering LuaFunctionHandle
		reg.add_function("FunctionWithFunctionArgument", &FunctionWithFunctionArgument, grp);
#endif

	//	registering add
		reg.add_function("add", static_cast<int (*)(int, int)>(
								&Add), grp, "c", "a#b");
		reg.add_function("add", static_cast<int (*)(int, int, int)>(
								&Add), grp, "d", "a#b#c");
		reg.add_function("add", static_cast<string (*)(const char*, const char*)>(
								&Add), grp, "c", "a#b");

	//	register class "Test"
		reg.add_class_<Test>("Test", grp)
			.add_constructor()
			.add_constructor<void (*)(const char*)>()
			.add_method("add", static_cast<int (Test::*)(int, int)>(&Test::add), "c", "a#b")
			.add_method("add", static_cast<int (Test::*)(int, int, int)>(&Test::add), "d", "a#b#c")
			.add_method("add", static_cast<string (Test::*)(const char*, const char*)>(&Test::add), "d", "a#b#c")
			.add_method("print_name", &Test::print_name)
			.add_method("print", static_cast<int(Test::*)()>(&Test::print))
			.add_method("print", static_cast<int(Test::*)() const>(&Test::print))
			.set_construct_as_smart_pointer(true);

	//	registering base class (without constructor)
		reg.add_class_<Base>("Base", grp)
			.add_constructor()
			.add_method("print", &Base::print)
			.set_construct_as_smart_pointer(true);

	//	registering derived class
		reg.add_class_<Derived, Base>("Derived", grp)
			.add_constructor()
			.add_method("print", &Derived::print)
			.set_construct_as_smart_pointer(true);

	//	registering derived class
		reg.add_class_<FurtherDerived, Derived>("FurtherDerived", grp)
			.add_constructor()
			.add_method("print", &FurtherDerived::print)
			.set_construct_as_smart_pointer(true);

		reg.add_function("CreateConstFurtherDerived", &CreateConstFurtherDerived, grp);

		reg.add_class_<Base0>("Base0", grp)
			.add_method("virt_print_base0", &Base0::virt_print_base0)
			.add_method("print_base0", &Base0::print_base0);

		reg.add_class_<Base1>("Base1", grp)
			.add_method("virt_print_base1", &Base1::virt_print_base1)
			.add_method("print_base1", &Base1::print_base1);

		reg.add_class_<Base2>("Base2", grp)
			.add_method("virt_print_base2", &Base2::virt_print_base2)
			.add_method("print_base2", &Base2::print_base2);

		reg.add_class_<Base3>("Base3", grp)
			.add_method("virt_print_base3", &Base3::virt_print_base3)
			.add_method("print_base3", &Base3::print_base3);

		reg.add_class_<Intermediate0, Base1, Base0>("Intermediate0", grp)
			.add_method("virt_print_intermediate0", &Intermediate0::virt_print_intermediate0)
			.add_method("print_intermediate0", &Intermediate0::print_intermediate0);

		reg.add_class_<Intermediate1, Base2, Base3>("Intermediate1", grp)
			.add_method("virt_print_intermediate1", &Intermediate1::virt_print_intermediate1)
			.add_method("print_intermediate1", &Intermediate1::print_intermediate1);

		reg.add_class_<MultipleDerived, Intermediate0, Intermediate1>("MultipleDerived", grp)
			.add_constructor()
			.add_method("print_mulitple_derived", &MultipleDerived::print_mulitple_derived)
			.set_construct_as_smart_pointer(true);

		reg.add_function("SmartMultipleDerivedImpl", SmartMultipleDerivedImpl, grp);

		reg.add_function("PrintFunctionIntermediate", static_cast<void (*)(SmartPtr<Intermediate1>)>(&PrintFunctionIntermediate), grp);
//		reg.add_function("PrintFunction", (void (*)(SmartPtr<Base3>))&PrintFunction, grp); // nasty example!
		reg.add_function("PrintFunction", static_cast<void (*)(SmartPtr<Base2>)>(&PrintFunction), grp);
		reg.add_function("PrintFunction", static_cast<void (*)(Base&)>(&PrintFunction), grp);
//		reg.add_function("PrintFunctionIntermediate", (void (*)(Intermediate1&))&PrintFunctionIntermediate, grp); // nasty example!
		reg.add_function("PrintFunctionIntermediate", static_cast<void (*)(Intermediate0&)>(&PrintFunctionIntermediate), grp);
		reg.add_function("PrintFunction", static_cast<void (*)(Base3&)>(&PrintFunction), grp);

		reg.add_class_<ConstClass>("ConstClass", grp)
			.add_constructor()
			.add_method("const_method", &ConstClass::const_method)
			.set_construct_as_smart_pointer(true);

		reg.add_class_<Piece>("Piece", grp)
			.add_constructor()
			.add_method("size", &Piece::size)
			.set_construct_as_smart_pointer(true);

		reg.add_class_<Cake>("Cake", grp)
			.add_constructor()
			.add_method("take_pieces", &Cake::take_pieces)
			.add_method("add_pieces", &Cake::add_pieces)
			.add_method("pieces_left", &Cake::pieces_left)
			.set_construct_as_smart_pointer(true);

		reg.add_function("TestFunc", TestFunc, grp)
			.add_function("SmartTestImpl", SmartTestImpl, grp)
			.add_function("ConstSmartTestImpl", ConstSmartTestImpl, grp)
			.add_function("SmartTestFunc", SmartTestFunc, grp)
			.add_function("ConstSmartTestFunc", ConstSmartTestFunc, grp)
			.add_function("ConstTestFunc", ConstTestFunc, grp)
			.add_function("ToConst", ToConst, grp)
			.add_function("StringTest", StringTest, grp)
			.add_function("StdStringTest", StdStringTest, grp)
			.add_function("PrintStringTest", PrintStringTest, grp)
			.add_function("TestPageContainer", TestPageContainer, grp);

		reg.add_class_<SmartTest>("SmartTest", grp)
			.add_constructor()
			.add_method("say_hello", &SmartTest::say_hello)
			.set_construct_as_smart_pointer(true);// always last!

		reg.add_class_<SmartTestDerived, SmartTest>("SmartTestDerived", grp)
			.add_constructor()
			.add_method("say_hello", &SmartTestDerived::say_hello)
			.set_construct_as_smart_pointer(true);// always last!

		reg.add_function("SmartTestArrived", &SmartTestArrived, grp);

		reg.add_class_<MessageHubTest>("MessageHubTest", grp)
			.add_constructor()
			.add_method("post_message", &MessageHubTest::post_message)
			.set_construct_as_smart_pointer(true);


	//	if the following registration is performed, the app should fail on startup,
	//	since the registered method takes an argument of an unregistered type.
		//reg.add_function("UnregisteredParameterTest", UnregisteredParameterTest);

	//	if the following registration is performed, the app should fail on startup,
	//	since the registered class/functon is ill named (only [a-zA-Z_][a-zA-Z_0-9]* allowed).
		//reg.add_class_<SmartTest>("SmartTest", grp);
		//reg.add_class_<NonAllowedName1>("BlimBlim$tshe");
		//reg.add_class_<NonAllowedName2>("5BlimBlimtshe");
		//reg.add_class_<NonAllowedName3>("Blim::Blimtshe");
		//reg.add_class_<NonAllowedName4>("Blim Blimtshe");
		//reg.add_class_<AllowName>("AllowName").add_method("$liadsf", &AllowName::non_allowed_method);
		//reg.add_function("Blub%brrr", &NonAllowedFct1);
		//reg.add_function("5Blubbrrr", &NonAllowedFct2);
		//reg.add_function("Blub_brr r", &NonAllowedFct3);
		//reg.add_function("NotAllowedParamPerValue", &NotAllowedParamPerValue);

		reg.add_function("PrintStdVectorBool", &PrintStdVector<bool>, grp);
		reg.add_function("PrintStdVectorInteger", &PrintStdVector<int>, grp);
		reg.add_function("PrintStdVectorSize_t", &PrintStdVector<size_t>, grp);
		reg.add_function("PrintStdVectorNumber", &PrintStdVector<number>, grp);
		reg.add_function("PrintStdVectorFloat", &PrintStdVector<float>, grp);
		reg.add_function("PrintStdVectorDouble", &PrintStdVector<double>, grp);
		reg.add_function("PrintStdVectorChar", &PrintStdVector<const char*>, grp);
		reg.add_function("PrintStdVectorString", &PrintStdVector<std::string>, grp);

		reg.add_function("PrintConstStdVectorRefBool", &PrintConstStdVectorRef<bool>, grp);
		reg.add_function("PrintConstStdVectorRefInteger", &PrintConstStdVectorRef<int>, grp);
		reg.add_function("PrintConstStdVectorRefSize_t", &PrintConstStdVectorRef<size_t>, grp);
		reg.add_function("PrintConstStdVectorRefNumber", &PrintConstStdVectorRef<number>, grp);
		reg.add_function("PrintConstStdVectorRefFloat", &PrintConstStdVectorRef<float>, grp);
		reg.add_function("PrintConstStdVectorRefDouble", &PrintConstStdVectorRef<double>, grp);
		reg.add_function("PrintConstStdVectorRefChar", &PrintConstStdVectorRef<const char*>, grp);
		reg.add_function("PrintConstStdVectorRefString", &PrintConstStdVectorRef<std::string>, grp);

		reg.add_function("PrintStdVector_BasePtr", &PrintStdVectorOfClass<Base*>, grp);
		reg.add_function("PrintStdVector_BaseConstPtr", &PrintStdVectorOfClass<const Base*>, grp);
		reg.add_function("PrintStdVector_BaseSmartPtr", &PrintStdVectorOfClass<SmartPtr<Base> >, grp);
		reg.add_function("PrintStdVector_BaseConstSmartPtr", &PrintStdVectorOfClass<ConstSmartPtr<Base> >, grp);

		reg.add_function("PrintConstStdVectorRef_BasePtr", &ConstStdVectorRefOfClass<Base*>, grp);
		reg.add_function("PrintConstStdVectorRef_BaseConstPtr", &ConstStdVectorRefOfClass<const Base*>, grp);
		reg.add_function("PrintConstStdVectorRef_BaseSmartPtr", &ConstStdVectorRefOfClass<SmartPtr<Base> >, grp);
		reg.add_function("PrintConstStdVectorRef_BaseConstSmartPtr", &ConstStdVectorRefOfClass<ConstSmartPtr<Base> >, grp);

		reg.add_function("ReturnConstStdVectorRef_String", &ReturnConstStdVectorRef_String, grp);
		reg.add_function("ReturnStdVector_String", &ReturnStdVector_String, grp);

		reg.add_function("ReturnConstStdVectorRef_Int", &ReturnConstStdVectorRef_Number<int>, grp);
		reg.add_function("ReturnStdVector_Int", &ReturnStdVector_Number<int>, grp);

		reg.add_function("ReturnConstStdVectorRef_Double", &ReturnConstStdVectorRef_Number<double>, grp);
		reg.add_function("ReturnStdVector_Double", &ReturnStdVector_Number<double>, grp);

		reg.add_function("ReturnConstStdVectorRef_BasePtr", &ReturnConstStdVectorRefOfClass<Base*>, grp);
		reg.add_function("ReturnStdVector_BasePtr", &ReturnStdVectorOfClass<Base*>, grp);

		reg.add_function("ReturnConstStdVectorRef_BaseConstPtr", &ReturnConstStdVectorRefOfClass<const Base*>, grp);
		reg.add_function("ReturnStdVector_BaseConstPtr", &ReturnStdVectorOfClass<const Base*>, grp);

		reg.add_function("ReturnConstStdVectorRef_BaseSmartPtr", &ReturnConstStdVectorRefOfClassSP<SmartPtr<Base> >, grp);
		reg.add_function("ReturnStdVector_BaseSmartPtr", &ReturnStdVectorOfClass<SmartPtr<Base> >, grp);

		reg.add_function("ReturnConstStdVectorRef_BaseConstSmartPtr", &ReturnConstStdVectorRefOfClassSP<ConstSmartPtr<Base> >, grp);
		reg.add_function("ReturnStdVector_BaseConstSmartPtr", &ReturnStdVectorOfClass<ConstSmartPtr<Base> >, grp);

	}
	UG_REGISTRY_CATCH_THROW(grp);
}

// end group test_bridge
/// \}

}//	end of namespace
}//	end of namespace
