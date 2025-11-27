/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef CONNECTION_H_
#define CONNECTION_H_

namespace ug{

template<typename T>
class AlgebraicConnection
{
public:
	size_t iIndex; // index to
	T dValue; // smallmatrix value;

	AlgebraicConnection() = default;
	AlgebraicConnection(size_t i, const T &v)
	: iIndex(i), dValue(v) {}

	void print() const {std::cout << *this;}

	friend std::ostream &operator << (std::ostream &output, const AlgebraicConnection&c)
	{
		output << "(" << c.iIndex << "-> ";
		output << c.dValue;
		output << ")";
		return output;
	}

	void operator = (const AlgebraicConnection &other)
	{
		iIndex = other.iIndex;
		dValue = other.dValue;
		// ø todo should normally return AlgebraicConnection& as a type
	}

	int operator < (const AlgebraicConnection &c) const
	{
		return iIndex < c.iIndex;
		// ø todo should normally return bool as a type to be compliant with stl
	}

	size_t &index()
	{
		return iIndex;
	}
	T &value()
	{
		return dValue;
	}
};

}
#endif