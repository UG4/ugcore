#ifndef __H__UG__SMALL_ALGEBRA__
#define __H__UG__SMALL_ALGEBRA__

/**
 * \defgroup small_algebra Small Algebra
 * \ingroup lib_algebra
 */

#include "blocks.h"


#ifdef LAPACK_AVAILABLE
#ifdef BLAS_AVAILABLE
#include "lapack/lapack.h"
#endif
#else
#include "no_lapack/no_lapack.h"
#endif

#include "small_matrix/densematrix_inverse.h"


#endif /* __H__UG__SMALL_ALGEBRA__ */
