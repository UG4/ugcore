/*
 * lib_algebra.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__
#define __H__LIB_ALGEBRA__

// other ug4 modules
#include "common/common.h"

#include "local_matrix_vector/flex_local_matrix_vector.h"

// library intern includes
#include "lib_algebra/multi_index/multi_indices.h"

#include "hypre_algebra/hyprematrix.h"
#include "hypre_algebra/hyprevector.h"
#include "hypre_algebra/hyprelinearsolver.h"

#include "arne_algebra/arnematrix.h"
#include "arne_algebra/arnevector.h"
#include "arne_algebra/arnelinearsolver.h"

namespace ug {

/** Define different algebra types.
 *  An Algebra should export the following typedef:
 *  - matrix_type
 *  - vector_type
 *  - index_type
 */
class ArneAlgebra{
	public:
		// matrix type
		typedef ArneMatrix matrix_type;

		// vector type
		typedef ArneVector vector_type;

		// index_type
		typedef MultiIndex<1> index_type;

		typedef ArneJacobi linear_solver_type;
};

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

#endif /* __H__LIB_ALGEBRA__ */
