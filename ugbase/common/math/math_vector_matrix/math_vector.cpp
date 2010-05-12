/*
 * math_vector.cpp
 *
 *  Created on: 08.07.2009
 *      Author: andreasvogel
 */

#include "math_vector.h"
#include  <iostream>

std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<2>& v)
{
	outStream << "(" << v.coord(0) << ", " << v.coord(1) << ")";
	return outStream;
}

std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<3>& v)
{
	outStream << "(" << v.coord(0) << ", " << v.coord(1) << ", " << v.coord(2) << ")";
	return outStream;
}

std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<4>& v)
{
	outStream << "(" << v.coord(0) << ", " << v.coord(1) << ", " << v.coord(2) << ", " << v.coord(3) << ")";
	return outStream;
}
