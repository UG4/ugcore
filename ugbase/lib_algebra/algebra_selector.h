/*
 * algebra_selector.h
 *
 *  Created on: 03.11.2010
 *      Author: Martin Rupp
 */

#ifndef __H__LIB_ALGEBRA__ALGEBRA_SELECTOR_
#define __H__LIB_ALGEBRA__ALGEBRA_SELECTOR_

#include <sstream>
#include <iostream>
#include "common/log.h"
#include "common/assert.h"

namespace ug{

/**
 * \brief Algebra Library
 *
 *
 * \defgroup lib_algebra lib_algebra
 */

////////////////////////////////////////////////////////////////////////////////
//   Algebra Type
////////////////////////////////////////////////////////////////////////////////

class AlgebraType
{
	public:
	///	types of algebra
		enum Type
		{
			CPU = 0,
		};

	///	indicating variable block size
		enum {VariableBlockSize = -1};

	public:
	///	constructor
		AlgebraType(Type type, int blockSize) : m_type(type), m_blockSize(blockSize)
		{
			if(blockSize <= 0 && blockSize != VariableBlockSize)
				throw(UGFatalError("BlockSize not allowed. Choose > 0 or VariableBlockSize"));
		}

	///	returns the type
		int type() const {return m_type;}

	///	returns the blocksize
		int blocksize() const {return m_blockSize;}

	protected:
		int m_type;
		int m_blockSize;
};

/// writes the Identifier to the output stream
inline std::ostream& operator<<(std::ostream& out,	const AlgebraType& v)
{
	std::stringstream ss;
	if(v.blocksize() >= 0) ss << v.blocksize();
	else if(v.blocksize() == AlgebraType::VariableBlockSize) ss << "variable";
	else ss << "invalid";

	switch(v.type())
	{
		case AlgebraType::CPU: out << "(CPU, " << ss.str() << ")"; break;
		default: out << "(unknown, " << ss.str() << ")";
	}
	return out;
}

/// enum for the different types of algebra provided by the library
enum enum_AlgebraType
{
	eCPUAlgebra = 0,
	eCPUBlockAlgebra2x2 = 1,
	eCPUBlockAlgebra3x3 = 2,
	eCPUBlockAlgebra4x4 = 3,
	eCPUVariableBlockAlgebra = 4
};

////////////////////////////////////////////////////////////////////////////////
//   Algebra Selector
////////////////////////////////////////////////////////////////////////////////

/// interface for the handling of the selection of algebra type
class IAlgebraTypeSelector
{
public:
	virtual ~IAlgebraTypeSelector() {}
	virtual AlgebraType get_algebra_type() const = 0;
};

/// handles the selection of algebra type
class CPUAlgebraSelector : public IAlgebraTypeSelector
{
	public:
		CPUAlgebraSelector() { m_variable = false; m_blocksize = 1;}
		virtual ~CPUAlgebraSelector() {}
		void set_fixed_blocksize(size_t blocksize)
		{
			m_blocksize = blocksize; m_variable = false;
			if(m_blocksize > 4 || m_blocksize == 0)
			{
				UG_LOG("CPUAlgebraSelector: Fixed blocksize " << m_blocksize <<
					   " not supported. Using variable blocks.\n");
				m_variable=true;
			}

		}

		void set_variable_blocksize() {m_variable = true;}

		virtual AlgebraType get_algebra_type() const
		{
			if(m_variable) return AlgebraType(AlgebraType::CPU, AlgebraType::VariableBlockSize);
			else return AlgebraType(AlgebraType::CPU, m_blocksize);
		}

	private:
		bool m_variable;
		size_t m_blocksize;
};

} // end namespace ug

#endif //  __H__LIB_ALGEBRA__ALGEBRA_SELECTOR_
