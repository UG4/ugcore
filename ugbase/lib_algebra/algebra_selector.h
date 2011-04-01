/*
 * algebra_selector.h
 *
 *  Created on: 03.11.2010
 *      Author: Martin Rupp
 */

#ifndef __H__LIB_ALGEBRA__ALGEBRA_SELECTOR_
#define __H__LIB_ALGEBRA__ALGEBRA_SELECTOR_

#include "common/log.h"
#include "common/assert.h"

namespace ug{

/**
 * \brief Algebra Library
 *
 *
 * \defgroup lib_algebra lib_algebra
 */



/////////////////////////////////////////////
/////////////////////////////////////////////
//   Algebra Selector
/////////////////////////////////////////////
/////////////////////////////////////////////


/// enum for the different types of algebra provided by the library
enum enum_AlgebraType
{
	eCPUAlgebra = 0,
	eCPUBlockAlgebra2x2 = 1,
	eCPUBlockAlgebra3x3 = 2,
	eCPUBlockAlgebra4x4 = 3,
	eCPUVariableBlockAlgebra = 4
};

/// interface for the handling of the selection of algebra type
class IAlgebraTypeSelector
{
public:
	virtual ~IAlgebraTypeSelector() {}
	virtual int get_algebra_type() = 0;
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
			UG_LOG("Fixed blocksize " << m_blocksize << " not supported. Using variable blocks.\n");
			m_variable=true;
		}

	}
	void set_variable_blocksize()
	{
		m_variable = true;
	}

	virtual int get_algebra_type()
	{
		if(m_variable)
			return eCPUVariableBlockAlgebra;
		else switch(m_blocksize)
		{
		case 1: 	return eCPUAlgebra;
		case 2: 	return eCPUBlockAlgebra2x2;
		case 3: 	return eCPUBlockAlgebra3x3;
		case 4: 	return eCPUBlockAlgebra4x4;
		default:	return eCPUVariableBlockAlgebra;
		}
	}
private:
	bool m_variable;
	size_t m_blocksize;
};

} // end namespace ug

#endif //  __H__LIB_ALGEBRA__ALGEBRA_SELECTOR_
