//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "../ugbridge/registry.h"

using namespace std;

namespace ug{
namespace interface
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
};


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
		virtual void print() = 0;
};

class Derived
{
	public:
		virtual void print()
		{
			UG_LOG("Derived::print() called\n");
		}
};

void PrintFunction(Base& b)
{
	b.print();
}

void RegisterTestInterface(InterfaceRegistry& reg)
{
	reg.add_function("add", &Add, "c", "a,b");

	reg.add_class_<Test>("Test")
		.add_constructor()
		.add_method("add", &Test::add, "c", "a,b")
		.add_method("print_name", &Test::print_name);

	reg.add_class_<Piece>("Piece")
		.add_constructor()
		.add_method("size", &Piece::size);

	reg.add_class_<Cake>("Cake")
		.add_constructor()
		.add_method("take_pieces", &Cake::take_pieces)
		.add_method("add_pieces", &Cake::add_pieces)
		.add_method("pieces_left", &Cake::pieces_left);

	reg.add_class_<Base>("Base")
		//.add_constructor()
		.add_method("print", &Base::print);
	reg.add_class_<Derived, Base>("Derived")
		.add_constructor()
		.add_method("print", &Derived::print);
	reg.add_function("PrintFunction", &PrintFunction);

}

}//	end of namespace
}//	end of namespace
