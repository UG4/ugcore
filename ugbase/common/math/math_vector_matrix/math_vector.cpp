/*
 * math_vector.cpp
 *
 *  Created on: 08.07.2009
 *      Author: andreasvogel
 */

#include "math_vector.h"
#include  <iostream>

namespace ug{

///	formatted output of MathVector objects: (...,...)
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<1>& v)
{
	outStream << "(" << v.coord(0) << ")";
	return outStream;
}

///	formatted output of MathVector objects: (...,...)
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<2>& v)
{
	outStream << "(" << v.coord(0) << ", " << v.coord(1) << ")";
	return outStream;
}

///	formatted output of MathVector objects: (...,...)
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<3>& v)
{
	outStream << "(" << v.coord(0) << ", " << v.coord(1) << ", " << v.coord(2) << ")";
	return outStream;
}

///	formatted output of MathVector objects: (...,...)
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<4>& v)
{
	outStream << "(" << v.coord(0) << ", " << v.coord(1) << ", " << v.coord(2) << ", " << v.coord(3) << ")";
	return outStream;
}

///	plain text output of MathVector objects: space-separated coordinates
std::ostream& write_plain_txt (std::ostream& outStream, const ug::MathVector<1>& v)
{
	outStream << v.coord(0);
	return outStream;
}

///	plain text output of MathVector objects: space-separated coordinates
std::ostream& write_plain_txt (std::ostream& outStream, const ug::MathVector<2>& v)
{
	outStream << v.coord(0) << ' ' << v.coord(1);
	return outStream;
}

///	plain text output of MathVector objects: space-separated coordinates
std::ostream& write_plain_txt (std::ostream& outStream, const ug::MathVector<3>& v)
{
	outStream << v.coord(0) << ' ' << v.coord(1) << ' ' << v.coord(2);
	return outStream;
}

///	plain text output of MathVector objects: space-separated coordinates
std::ostream& write_plain_txt (std::ostream& outStream, const ug::MathVector<4>& v)
{
	outStream << v.coord(0) << ' ' << v.coord(1) << ' ' << v.coord(2) << ' ' << v.coord(3);
	return outStream;
}

///	plain text input of MathVectors: read space-separated coordinates
std::istream& read_plain_txt (std::istream& inStream, ug::MathVector<1>& v)
{
	inStream >> v.coord(0);
	return inStream;
}

///	plain text input of MathVectors: read space-separated coordinates
std::istream& read_plain_txt (std::istream& inStream, ug::MathVector<2>& v)
{
	inStream >> v.coord(0) >> v.coord(1);
	return inStream;
}

///	plain text input of MathVectors: read space-separated coordinates
std::istream& read_plain_txt (std::istream& inStream, ug::MathVector<3>& v)
{
	inStream >> v.coord(0) >> v.coord(1) >> v.coord(2);
	return inStream;
}

///	plain text input of MathVectors: read space-separated coordinates
std::istream& read_plain_txt (std::istream& inStream, ug::MathVector<4>& v)
{
	inStream >> std::ws >> v.coord(0) >> std::ws >> v.coord(1) >> std::ws >> v.coord(2) >> std::ws >> v.coord(3);
	return inStream;
}

}
