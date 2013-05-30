/*
 * algebra_type.cpp
 *
 *  Created on: 06.10.2011
 *      Author: andreasvogel
 */

#include "algebra_type.h"
#include "common/common.h"

namespace ug{

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
	else UG_THROW("Algebra Type '"<<sType<<"' not reconized. Available: CPU, GPU.");
}

AlgebraType::AlgebraType(const char* type)
	: m_blockSize(VariableBlockSize)
{
	std::string sType(type);

	if(sType == "CPU") m_type = CPU;
	else if(sType == "GPU") m_type = GPU;
	else if(sType == "CRS") { UG_THROW("Type CRS is deprecated, use CPU instead."); }
	else UG_THROW("Algebra Type '"<<sType<<"' not reconized. Available: CPU, CRS.");
}


inline std::ostream& operator<<(std::ostream& out,	const AlgebraType& v)
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
