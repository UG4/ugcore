/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include "algebra_type.h"

#include "common/common.h"

namespace ug {

AlgebraType::AlgebraType()
	: m_type(CPU), m_blockSize(VariableBlockSize)
{
}

AlgebraType::AlgebraType(Type type, int blockSize)
	: m_type(type), m_blockSize(blockSize)
{
	if(blockSize <= 0 && blockSize != VariableBlockSize)
		UG_THROW("BlockSize not allowed. Choose > 0 or VariableBlockSize");
}

AlgebraType::AlgebraType(const char* type, int blockSize)
	: m_blockSize(blockSize)
{
	if(blockSize <= 0 && blockSize != VariableBlockSize)
		UG_THROW("BlockSize not allowed. Choose > 0 or VariableBlockSize");

	std::string sType(type);

	if(sType == "CPU") m_type = CPU;
	else if(sType == "GPU") m_type = GPU;
	else if(sType == "CRS") { UG_THROW("Type CRS is deprecated, use CPU instead."); }
	else UG_THROW("Algebra Type '"<<sType<<"' not recognized. Available: CPU, GPU.");
}

AlgebraType::AlgebraType(const char* type)
	: m_blockSize(VariableBlockSize)
{
	std::string sType(type);

	if(sType == "CPU") m_type = CPU;
	else if(sType == "GPU") m_type = GPU;
	else if(sType == "CRS") { UG_THROW("Type CRS is deprecated, use CPU instead."); }
	else UG_THROW("Algebra Type '"<<sType<<"' not recognized. Available: CPU, CRS.");
}


inline std::ostream& operator << (std::ostream& out,	const AlgebraType& v)
{
	std::stringstream ss;
	if(v.blocksize() >= 0) ss << v.blocksize();
	else if(v.blocksize() == AlgebraType::VariableBlockSize) ss << "variable";
	else ss << "invalid";

	switch(v.type())
	{
		case AlgebraType::CPU: out << "(CPU, " << ss.str() << ")"; break;
		case AlgebraType::GPU: out << "(GPU, " << ss.str() << ")"; break;
		default: out << "(unknown, " << ss.str() << ")";
	}
	return out;
}

AlgebraType DefaultAlgebra::m_default = AlgebraType(AlgebraType::CPU, 1);

} // end namespace ug
