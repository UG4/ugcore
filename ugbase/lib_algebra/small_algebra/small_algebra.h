#include "blocks.h"


#ifdef LAPACK_AVAILABLE
#ifdef BLAS_AVAILABLE
#include "lapack/lapack.h"
#endif
#else
#include "no_lapack/no_lapack.h"
#endif

#include "small_matrix/densematrix_inverse.h"

