//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "../registry.h"
#include "../ug_bridge.h"
#include <iostream>
#include <sstream>

using namespace std;

namespace ug
{
namespace bridge
{

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
		virtual void print() = 0;
};

// Some Derived class
class Derived : public Base
{
	public:
		virtual ~Derived()	{}
		virtual void print()
		{
			UG_LOG("Derived::print() called\n");
		}
};


class Base0
{
	public:
		virtual ~Base0() {}
		virtual void virt_print_base0() = 0;
		void print_base0() {UG_LOG("Base0::print_base0() called.\n");}
};

class Base1
{
	public:
		virtual ~Base1() {}
		virtual void virt_print_base1() = 0;
		void print_base1() {UG_LOG("Base1::print_base1() called.\n");}
};

class Base2
{
	public:
		virtual ~Base2() {}
		virtual void virt_print_base2() = 0;
		void print_base2() {UG_LOG("Base2::print_base2() called.\n");}
};

class Base3
{
	public:
		virtual ~Base3() {}
		virtual void virt_print_base3() = 0;
		void print_base3() {UG_LOG("Base3::print_base3() called.\n");}
};

class Intermediate0 : public Base0, public Base1
{
	public:
		virtual ~Intermediate0() {}
		virtual void virt_print_intermediate0() = 0;
		void print_intermediate0() {UG_LOG("Intermediate0::print_intermediate0() called.\n");}
};

class Intermediate1 : public Base2, public Base3
{
	public:
		virtual ~Intermediate1() {}
		virtual void virt_print_intermediate1() = 0;
		void print_intermediate1() {UG_LOG("Intermediate1::print_intermediate1() called.\n");}
};

class MultipleDerived : public Intermediate0, public Intermediate1
{
	public:
		virtual ~MultipleDerived() {}
		void print_mulitple_derived(){UG_LOG("MultipleDerived::print() called\n");}

		virtual void virt_print_intermediate0()	{UG_LOG("MultipleDerived::virt_print_intermediate0() called\n");}
		virtual void virt_print_intermediate1()	{UG_LOG("MultipleDerived::virt_print_intermediate1() called\n");}

		virtual void virt_print_base0()	{UG_LOG("MultipleDerived::virt_print_base0() called\n");}
		virtual void virt_print_base1()	{UG_LOG("MultipleDerived::virt_print_base1() called\n");}
		virtual void virt_print_base2()	{UG_LOG("MultipleDerived::virt_print_base2() called\n");}
		virtual void virt_print_base3()	{UG_LOG("MultipleDerived::virt_print_base3() called\n");}
};

SmartPtr<MultipleDerived> SmartMultipleDerivedImpl(){
	return SmartPtr<MultipleDerived>(new MultipleDerived());
}

void PrintFunction(SmartPtr<Base3> b)
{
	b->virt_print_base3();
	b->print_base3();
}

void PrintFunctionIntermediate(SmartPtr<Intermediate1> b)
{
	b->virt_print_intermediate1();
	b->print_intermediate1();
}

class ConstClass
{
	public:
		std::string const_method() const
		{
			return "const_method_called";
		}
};

void PrintFunctionIntermediate(Intermediate1& b)
{
	b.virt_print_intermediate1();
	b.print_intermediate1();
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

std::string StdStringTest()
{
	return std::string("stdJooik");
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

bool RegisterTestInterface(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		std::stringstream groupString; groupString << parentGroup << "/Test";
		std::string grp = groupString.str();

	//	registering hello world
		reg.add_function("PrintHelloWorld", &PrintHelloWorldToScreen);

	//	registering add
		reg.add_function("add", (int (*)(int, int))
								&Add, grp.c_str(), "c", "a#b");
		reg.add_function("add", (int (*)(int, int, int))
								&Add, grp.c_str(), "d", "a#b#c");
		reg.add_function("add", (string (*)(const char*, const char*))
								&Add, grp.c_str(), "c", "a#b");

	//	register class "Test"
		reg.add_class_<Test>("Test", grp.c_str())
			.add_constructor()
			.add_method("add", (int (Test::*)(int, int))&Test::add, "c", "a#b")
			.add_method("add", (int (Test::*)(int, int, int))&Test::add, "d", "a#b#c")
			.add_method("add", (string (Test::*)(const char*, const char*))&Test::add, "d", "a#b#c")
			.add_method("print_name", &Test::print_name)
			.add_method("print", (int(Test::*)()) &Test::print)
			.add_method("print", (int(Test::*)() const) &Test::print);

	//	registering base class (without constructor)
		reg.add_class_<Base>("Base", grp.c_str())
			.add_method("print", &Base::print);

	//	registering derived class
		reg.add_class_<Derived, Base>("Derived", grp.c_str())
			.add_constructor();


		reg.add_class_<Base0>("Base0", grp.c_str())
			.add_method("virt_print_base0", &Base0::virt_print_base0)
			.add_method("print_base0", &Base0::print_base0);

		reg.add_class_<Base1>("Base1", grp.c_str())
			.add_method("virt_print_base1", &Base1::virt_print_base1)
			.add_method("print_base1", &Base1::print_base1);

		reg.add_class_<Base2>("Base2", grp.c_str())
			.add_method("virt_print_base2", &Base2::virt_print_base2)
			.add_method("print_base2", &Base2::print_base2);

		reg.add_class_<Base3>("Base3", grp.c_str())
			.add_method("virt_print_base3", &Base3::virt_print_base3)
			.add_method("print_base3", &Base3::print_base3);

		reg.add_class_<Intermediate0, Base1, Base0>("Intermediate0", grp.c_str())
			.add_method("virt_print_intermediate0", &Intermediate0::virt_print_intermediate0)
			.add_method("print_intermediate0", &Intermediate0::print_intermediate0);

		reg.add_class_<Intermediate1, Base2, Base3>("Intermediate1", grp.c_str())
			.add_method("virt_print_intermediate1", &Intermediate1::virt_print_intermediate1)
			.add_method("print_intermediate1", &Intermediate1::print_intermediate1);

		reg.add_class_<MultipleDerived, Intermediate0, Intermediate1>("MultipleDerived", grp.c_str())
			.add_constructor()
			.add_method("print_mulitple_derived", &MultipleDerived::print_mulitple_derived);

		reg.add_function("SmartMultipleDerivedImpl", SmartMultipleDerivedImpl);

		reg.add_function("PrintFunctionIntermediate", (void (*)(SmartPtr<Intermediate1>))&PrintFunctionIntermediate);
		reg.add_function("PrintFunction", (void (*)(SmartPtr<Base3>))&PrintFunction);
		reg.add_function("PrintFunction", (void (*)(Base&))&PrintFunction);
		reg.add_function("PrintFunctionIntermediate", (void (*)(Intermediate1&))&PrintFunctionIntermediate);
		reg.add_function("PrintFunction", (void (*)(Base3&))&PrintFunction);

		reg.add_class_<ConstClass>("ConstClass", grp.c_str())
			.add_constructor()
			.add_method("const_method", &ConstClass::const_method);

		reg.add_class_<Piece>("Piece", grp.c_str())
			.add_constructor()
			.add_method("size", &Piece::size);

		reg.add_class_<Cake>("Cake", grp.c_str())
			.add_constructor()
			.add_method("take_pieces", &Cake::take_pieces)
			.add_method("add_pieces", &Cake::add_pieces)
			.add_method("pieces_left", &Cake::pieces_left);

		reg.add_function("TestFunc", TestFunc, grp.c_str())
			.add_function("SmartTestImpl", SmartTestImpl, grp.c_str())
			.add_function("ConstSmartTestImpl", ConstSmartTestImpl, grp.c_str())
			.add_function("SmartTestFunc", SmartTestFunc, grp.c_str())
			.add_function("ConstSmartTestFunc", ConstSmartTestFunc, grp.c_str())
			.add_function("ConstTestFunc", ConstTestFunc, grp.c_str())
			.add_function("ToConst", ToConst, grp.c_str())
			.add_function("StringTest", StringTest, grp.c_str())
			.add_function("StdStringTest", StdStringTest, grp.c_str());
			
		reg.add_function("PostRegisterTest", &PostRegisterTest);
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterTestInterface: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

}//	end of namespace
}//	end of namespace
