/*
 * algebra_chooser.h
 *
 *  Created on: 03.11.2010
 *      Author: Martin Rupp
 */

#ifndef __H__LIB_ALGEBRA__ALGEBRA_CHOOSER_
#define __H__LIB_ALGEBRA__ALGEBRA_CHOOSER_

// vector and sparse_matrix
#include "cpu_algebra/vector.h"
#include "cpu_algebra/sparsematrix.h"

// parallel support
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallel_vector.h"
	#include "lib_algebra/parallelization/parallel_matrix.h"
#endif


namespace ug{

/**
 * \brief Algebra Library
 *
 *
 * \defgroup lib_algebra lib_algebra
 */

/////////////////////////////////////////////
/////////////////////////////////////////////
//   Algebra Types
/////////////////////////////////////////////
/////////////////////////////////////////////

/*  Define different algebra types.
 *  An Algebra should export the following typedef:
 *  - matrix_type
 *  - vector_type
 */

/////////////////////////////////////////////
//   CPU Algebra (Block 1x1 Algebra)
/////////////////////////////////////////////

class CPUAlgebra
{
public:
#ifdef UG_PARALLEL
		typedef ParallelMatrix<SparseMatrix<double> > matrix_type;
		typedef ParallelVector<Vector<double> > vector_type;
#else
		typedef SparseMatrix<double> matrix_type;
		typedef Vector<double> vector_type;
#endif
};

/////////////////////////////////////////////
//   CPU Fixed Block Algebra
/////////////////////////////////////////////
template<int TBlockSize>
class CPUBlockAlgebra
{
public:
#ifdef UG_PARALLEL
	typedef ParallelMatrix<SparseMatrix<DenseMatrix<FixedArray2<double, TBlockSize, TBlockSize> > > > matrix_type;
	typedef ParallelVector<Vector<DenseVector<FixedArray1<double, TBlockSize> > > > vector_type;
#else
	typedef  SparseMatrix<DenseMatrix<FixedArray2<double, TBlockSize, TBlockSize> > > matrix_type;
	typedef Vector<DenseVector<FixedArray1<double, TBlockSize> > > vector_type;
#endif
};


/////////////////////////////////////////////
//   CPU Variable Block Algebra
/////////////////////////////////////////////

class CPUVariableBlockAlgebra
{
public:
#ifdef UG_PARALLEL
	typedef ParallelMatrix<SparseMatrix<DenseMatrix<VariableArray2<double> > > > matrix_type;
	typedef ParallelVector<Vector<DenseVector<VariableArray1<double> > > > vector_type;
#else
	typedef  SparseMatrix<DenseMatrix<VariableArray2<double> > > matrix_type;
	typedef Vector<DenseVector<VariableArray1<double> > > vector_type;
#endif
};


/////////////////////////////////////////////
/////////////////////////////////////////////
//   Algebra Chooser
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
class AlgebraTypeChooserInterface
{
public:
	virtual ~AlgebraTypeChooserInterface() {}
	virtual int get_algebra_type() = 0;
};

/// handles the selection of algebra type
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

} // end namespace ug

#endif //  __H__LIB_ALGEBRA__ALGEBRA_CHOOSER_
