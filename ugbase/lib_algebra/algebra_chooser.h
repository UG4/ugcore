/*
 * algebra_chooser.h
 *
 *  Created on: 03.11.2010
 *      Author: Martin Rupp
 */

#ifndef __H__LIB_ALGEBRA__ALGEBRA_CHOOSER_
#define __H__LIB_ALGEBRA__ALGEBRA_CHOOSER_

enum enum_AlgebraType
{
	eCPUAlgebra = 0,
	eCPUBlockAlgebra2x2 = 1,
	eCPUBlockAlgebra3x3 = 2,
	eCPUBlockAlgebra4x4 = 3,
	eCPUVariableBlockAlgebra = 4
};


class AlgebraTypeChooserInterface
{
public:
	virtual ~AlgebraTypeChooserInterface() {}
	virtual int get_algebra_type() = 0;
};

class CPUAlgebraChooser : public AlgebraTypeChooserInterface
{
public:
	CPUAlgebraChooser() { m_variable = false; m_blocksize = 1;}
	virtual ~CPUAlgebraChooser() {}
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

#endif //  __H__LIB_ALGEBRA__ALGEBRA_CHOOSER_
