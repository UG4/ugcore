/*
 * cpu_algebra_types.h
 *
 *  Created on: 01.04.2011
 *      Author: mrupp
 */

#ifndef __H__UG__LIB_ALGEBRA__CPU_ALGEBRA_TYPES__
#define __H__UG__LIB_ALGEBRA__CPU_ALGEBRA_TYPES__

#include "algebra_type.h"

// vector and sparse_matrix
#include "cpu_algebra/vector.h"
#include "cpu_algebra/sparsematrix.h"

#ifdef UG_GPU
//#include "gpu_algebra/gpuvector.h"
//#include "gpu_algebra/gpusparsematrix.h"
#endif

// parallel support
#ifdef UG_PARALLEL
	#include "parallelization/parallel_vector.h"
	#include "parallelization/parallel_matrix.h"
#endif

#include "small_algebra/small_algebra.h"
namespace ug{

////////////////////////////////////////////////////////////////////////////////
//   CPU Algebra
////////////////////////////////////////////////////////////////////////////////

/*  Define different algebra types.
 *  An Algebra should export the following typedef:
 *  - matrix_type
 *  - vector_type
 */

////////////////////////////////////////////////////////////////////////////////
//   CPU Algebra (Block 1x1 Algebra)
////////////////////////////////////////////////////////////////////////////////

/**
 * \defgroup cpu_algebra CPU Algebra
 * \ingroup lib_algebra
 * \{
 */

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
	static AlgebraType get_type()
	{
		return AlgebraType(AlgebraType::CPU, 1);
	}
};

// end group cpu_algebra
/// \}

/**
 * \defgroup crs_algebra CRS Algebra
 * \ingroup lib_algebra
 * \{
 */
/*#ifdef UG_GPU
struct GPUAlgebra
{
#ifdef UG_PARALLEL
		typedef ParallelMatrix<GPUSparseMatrix<double> > matrix_type;
		typedef ParallelVector<GPUVector<double> > vector_type;
#else
		typedef GPUSparseMatrix<double> matrix_type;
		typedef GPUVector<double> vector_type;
#endif

	static const int blockSize = 1;
	static AlgebraType get_type()
	{
		return AlgebraType(AlgebraType::GPU, 1);
	}
};
#endif
// end group crs_algebra
/// \}
*/
////////////////////////////////////////////////////////////////////////////////
//   CPU Fixed Block Algebra
////////////////////////////////////////////////////////////////////////////////
/// \addtogroup cpu_algebra
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
	static AlgebraType get_type()
	{
		return AlgebraType(AlgebraType::CPU, TBlockSize);
	}
};

////////////////////////////////////////////////////////////////////////////////
//   CPU Variable Block Algebra
////////////////////////////////////////////////////////////////////////////////

/// \addtogroup cpu_algebra
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
	static AlgebraType get_type()
	{
		return AlgebraType(AlgebraType::CPU, AlgebraType::VariableBlockSize);
	}
};

}
#endif /* __H__UG__LIB_ALGEBRA__CPU_ALGEBRA_TYPES__ */
