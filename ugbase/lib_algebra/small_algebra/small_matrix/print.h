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

#ifndef SMALL_ALGEBRA_PRINT_H_
#define SMALL_ALGEBRA_PRINT_H_

#include "densematrix.h"
#include <string>
#include <sstream>
#include "common/log.h"

namespace ug{
template<typename TStorage>
void MaplePrint(const DenseMatrix<TStorage> &A, const char *name)
{
	UG_LOG(name << " = " << MapleString(A));
}

template<typename TStorage>
std::string NonzeroString(const DenseMatrix<TStorage> &A, const char *brackets, const char *seperator, int base, const char *name)
{
	std::stringstream ss;
	ss << std::setprecision(9);
	for(size_t r=0; r<A.num_rows(); ++r)
	{
		for(size_t c=0; c<A.num_cols(); ++c)
			if(A(r, c) != 0.0)
				ss << name << brackets[0] << r+base << ", " << c+base << brackets[1] << " = " << A(r, c) << seperator;
		ss << "\n";
	}
	return ss.str();
}

template<typename TStorage>
std::string
MatlabString(const DenseMatrix<TStorage> &A, const char *name)
{
	std::stringstream ss;
	ss << name << " = zeros(" << A.num_rows() << ", " << A.num_cols() << ")\n";
	ss << NonzeroString(A, "()", "\n", 1, name);
	return ss.str();
}


template<typename TStorage>
std::string
JuliaString(const DenseMatrix<TStorage> &A, const char *name)
{
	std::stringstream ss;
	ss << name << " = zeros(" << A.num_rows() << ", " << A.num_cols() << ")\n";
	ss << NonzeroString(A, "[]", "\n", 1, name);
	return ss.str();
}

template<typename TStorage>
std::string
CPPString(const DenseMatrix<TStorage> &A, const char *name)
{
	std::stringstream ss;
	ss << name << " = Matrix(" << A.num_rows() << ", " << A.num_cols() << ")\n";
	ss << NonzeroString(A, "()", ";\n", 0, name);
	return ss.str();
}

}
#endif