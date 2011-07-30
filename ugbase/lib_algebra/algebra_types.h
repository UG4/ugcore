/*
 * algebra_types.h
 *
 *  Created on: 01.04.2011
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__ALGEBRA_TYPES__
#define __H__LIB_ALGEBRA__ALGEBRA_TYPES__

// vector and sparse_matrix
#include "cpu_algebra/vector.h"
#include "cpu_algebra/sparsematrix.h"

// parallel support
#ifdef UG_PARALLEL
	#include "parallelization/parallel_vector.h"
	#include "parallelization/parallel_matrix.h"
#endif

#include "algebra_selector.h"

namespace ug{

/**
 * \brief Algebra Library
 *
 *
 * \defgroup lib_algebra lib_algebra
 */
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//   Algebra Types
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*  Define different algebra types.
 *  An Algebra should export the following typedef:
 *  - matrix_type
 *  - vector_type
 */

////////////////////////////////////////////////////////////////////////////////
//   CPU Algebra (Block 1x1 Algebra)
////////////////////////////////////////////////////////////////////////////////

struct CPUAlgebra
{
#ifdef UG_PARALLEL
		typedef ParallelMatrix<SparseMatrix<double> > matrix_type;
		typedef ParallelVector<Vector<double> > vector_type;
#else
		typedef SparseMatrix<double> matrix_type;
		typedef Vector<double> vector_type;
#endif

	static const int blockSize = 1;
};

////////////////////////////////////////////////////////////////////////////////
//   CPU Fixed Block Algebra
////////////////////////////////////////////////////////////////////////////////
template<int TBlockSize>
struct CPUBlockAlgebra
{
#ifdef UG_PARALLEL
	typedef ParallelMatrix<SparseMatrix<DenseMatrix<FixedArray2<double, TBlockSize, TBlockSize> > > > matrix_type;
	typedef ParallelVector<Vector<DenseVector<FixedArray1<double, TBlockSize> > > > vector_type;
#else
	typedef  SparseMatrix<DenseMatrix<FixedArray2<double, TBlockSize, TBlockSize> > > matrix_type;
	typedef Vector<DenseVector<FixedArray1<double, TBlockSize> > > vector_type;
#endif

	static const int blockSize = TBlockSize;
};

////////////////////////////////////////////////////////////////////////////////
//   CPU Variable Block Algebra
////////////////////////////////////////////////////////////////////////////////

struct CPUVariableBlockAlgebra
{
#ifdef UG_PARALLEL
	typedef ParallelMatrix<SparseMatrix<DenseMatrix<VariableArray2<double> > > > matrix_type;
	typedef ParallelVector<Vector<DenseVector<VariableArray1<double> > > > vector_type;
#else
	typedef  SparseMatrix<DenseMatrix<VariableArray2<double> > > matrix_type;
	typedef Vector<DenseVector<VariableArray1<double> > > vector_type;
#endif

	static const int blockSize = AlgebraType::VariableBlockSize;
};

////////////////////////////////////////////////////////////////////////////////
// AlgebraEnum to Algebra Class Type
////////////////////////////////////////////////////////////////////////////////

template<enum_AlgebraType type>
class AlgebraEnumToType;

template <>
class AlgebraEnumToType<eCPUAlgebra> : public CPUAlgebra { };

template <>
class AlgebraEnumToType<eCPUBlockAlgebra2x2> : public CPUBlockAlgebra<2> { };

template <>
class AlgebraEnumToType<eCPUBlockAlgebra3x3> : public CPUBlockAlgebra<3> { };

template <>
class AlgebraEnumToType<eCPUBlockAlgebra4x4> : public CPUBlockAlgebra<4> { };

template <>
class AlgebraEnumToType<eCPUVariableBlockAlgebra> : public CPUVariableBlockAlgebra { };

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__ALGEBRA_TYPES__ */
