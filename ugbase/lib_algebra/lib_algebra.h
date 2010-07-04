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

#include "local_matrix_vector/flex_local_matrix_vector.h"

// library intern includes
#include "lib_algebra/multi_index/multi_indices.h"

// parallel support
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

/////////////////////////////////////////////
/////////////////////////////////////////////
//   Algebra
/////////////////////////////////////////////
/////////////////////////////////////////////


/** Define different algebra types.
 *  An Algebra should export the following typedef:
 *  - matrix_type
 *  - vector_type
 *  - index_type
 */

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

		// index_type
		typedef MultiIndex<1> index_type;
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

		// index_type
		typedef MultiIndex<1> index_type;

		typedef HYPREboomerAMG linear_solver_type;

};
}

#endif

/////////////////////////////////////////////
//   Martin Algebra
/////////////////////////////////////////////

//#define USE_MARTIN_ALGEBRA
#ifdef USE_MARTIN_ALGEBRA
#include "martin_algebra/vector.h"

#ifdef LAPACK_AVAILABLE
#ifdef BLAS_AVAILABLE
#include "martin_algebra/lapack_lu.h"
#endif
#endif

namespace ug
{
class MartinAlgebra
	{
	public:
		// matrix type
		typedef SparseMatrix<double> matrix_type;

		// vector type
#ifdef UG_PARALLEL
		typedef ParallelVector<Vector<double> > vector_type;
#else
		typedef Vector<double> vector_type;
#endif
		// index_type
		typedef MultiIndex<1> index_type;

		//	typedef HYPREboomerAMG linear_solver_type;
	};
}

#endif

#endif /* __H__LIB_ALGEBRA__LIB_ALGEBRA__ */
