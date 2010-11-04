/*
 * lib_algebra.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__LIB_ALGEBRA__
#define __H__LIB_ALGEBRA__LIB_ALGEBRA__

#include <iomanip>

// other ug4 modules
#include "common/common.h"


// operator interface
#include "operator/operator.h"

#include "common/operations.h"

#include "algebra_chooser.h"

/////////////////////////////////////////////
/////////////////////////////////////////////
//   Algebra
/////////////////////////////////////////////
/////////////////////////////////////////////


/////////////////////////////////////////////
//   small_algebra
/////////////////////////////////////////////

#include "small_algebra/small_algebra.h"

#ifdef LAPACK_AVAILABLE
#ifdef BLAS_AVAILABLE
#include "small_algebra/lapack/lapack.h"
#endif
#endif

/** Define different algebra types.
 *  An Algebra should export the following typedef:
 *  - matrix_type
 *  - vector_type
 *  - index_type
 */

/////////////////////////////////////////////
//   CPU Algebra
/////////////////////////////////////////////

#include "cpu_algebra/algebra_misc.h"
#include "cpu_algebra/vector.h"
#include "cpu_algebra/sparsematrix.h"
#include "cpu_algebra/core_smoothers.h"

#ifdef LAPACK_AVAILABLE
#ifdef BLAS_AVAILABLE
#include "cpu_algebra/lapack_lu.h"
#endif
#endif

// parallel support
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug
{

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

} // end namespace ug


#endif /* __H__LIB_ALGEBRA__LIB_ALGEBRA__ */
