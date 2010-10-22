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

int Add(int a, int b)
{
	return a + b;
}

class Test
{
	public:
		int add(int a, int b)
		{
			return a+b;
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

class Base
{
	public:
		virtual ~Base()	{}
		virtual void print() = 0;
};

class Derived
{
	public:
		virtual ~Derived()	{}
		virtual void print()
		{
			UG_LOG("Derived::print() called\n");
		}
};

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

bool RegisterTestInterface(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		std::stringstream groupString; groupString << parentGroup << "/Test";
		std::string grp = groupString.str();

		reg.add_function("add", &Add, grp.c_str(), "c", "a,b");

		reg.add_class_<Test>("Test", grp.c_str())
			.add_constructor()
			.add_method("add", &Test::add, "c", "a,b")
			.add_method("print_name", &Test::print_name)
			.add_method("print", (int(Test::*)()) &Test::print)
			.add_method("print", (int(Test::*)() const) &Test::print);

		reg.add_class_<Piece>("Piece", grp.c_str())
			.add_constructor()
			.add_method("size", &Piece::size);

		reg.add_class_<Cake>("Cake", grp.c_str())
			.add_constructor()
			.add_method("take_pieces", &Cake::take_pieces)
			.add_method("add_pieces", &Cake::add_pieces)
			.add_method("pieces_left", &Cake::pieces_left);

		reg.add_class_<Base>("Base", grp.c_str())
			//.add_constructor()
			.add_method("print", &Base::print);
		reg.add_class_<Derived, Base>("Derived", grp.c_str())
			.add_constructor()
			.add_method("print", &Derived::print);
		reg.add_function("PrintFunction", &PrintFunction);

		reg.add_function("TestFunc", TestFunc, grp.c_str())
			.add_function("ConstTestFunc", ConstTestFunc, grp.c_str())
			.add_function("ToConst", ToConst, grp.c_str())
			.add_function("StringTest", StringTest, grp.c_str())
			.add_function("StdStringTest", StdStringTest, grp.c_str());
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
