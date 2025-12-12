/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#include "math_vector.h"

#include  <iostream>

namespace ug {

///	formatted output of MathVector objects: (...,...)
std::ostream& operator << (std::ostream& outStream, const ug::MathVector<1>& v)
{
	outStream << "(" << v.coord(0) << ")";
	return outStream;
}

///	formatted output of MathVector objects: (...,...)
std::ostream& operator << (std::ostream& outStream, const ug::MathVector<2>& v)
{
	outStream << "(" << v.coord(0) << ", " << v.coord(1) << ")";
	return outStream;
}

///	formatted output of MathVector objects: (...,...)
std::ostream& operator << (std::ostream& outStream, const ug::MathVector<3>& v)
{
	outStream << "(" << v.coord(0) << ", " << v.coord(1) << ", " << v.coord(2) << ")";
	return outStream;
}

///	formatted output of MathVector objects: (...,...)
std::ostream& operator << (std::ostream& outStream, const ug::MathVector<4>& v)
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
