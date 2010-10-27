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
class CPUAlgebra
	{
	public:
		// matrix type
#ifdef UG_PARALLEL
		typedef ParallelMatrix<SparseMatrix<double> > matrix_type;
#else
		typedef SparseMatrix<double> matrix_type;
#endif
		// vector type
#ifdef UG_PARALLEL
		typedef ParallelVector<Vector<double> > vector_type;
#else
		typedef Vector<double> vector_type;
#endif
	};
}


/////////////////////////////////////////////
//   ublas algebra
/////////////////////////////////////////////

#include "ublas_algebra/ublas_matrix.h"
#include "ublas_algebra/ublas_vector.h"
#include "ublas_algebra/ublas_linearsolver.h"

namespace ug {

class UblasAlgebra{
	public:
		// matrix type
		typedef UblasMatrix matrix_type;

		// vector type
#ifdef UG_PARALLEL
		typedef ParallelVector<UblasVector> vector_type;
#else
		typedef UblasVector vector_type;
#endif
};

} // namespace ug

/////////////////////////////////////////////
//   Hypre Algebra
/////////////////////////////////////////////


#if HYPRELIB_DIR

#include "hypre_algebra/hyprematrix.h"
#include "hypre_algebra/hyprevector.h"
#include "hypre_algebra/hyprelinearsolver.h"

namespace ug{
class HypreAlgebra{
	public:
		// matrix type
		typedef HypreMatrix matrix_type;

		// vector type
		typedef HypreVector vector_type;

		typedef HYPREboomerAMG linear_solver_type;

};
}

#endif


#endif /* __H__LIB_ALGEBRA__LIB_ALGEBRA__ */
