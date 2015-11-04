
#ifndef __H__LIB_ALGEBRA__LIB_ALGEBRA_IMPL__
#define __H__LIB_ALGEBRA__LIB_ALGEBRA_IMPL__

//////////////////////////////////////////////////////////////////////////////
//
// This file is intended to include implementations of lib_algebra
//
//////////////////////////////////////////////////////////////////////////////


#include "lib_algebra.h"


// preconditioner
#include "lib_algebra/operator/preconditioner/jacobi.h"
#include "lib_algebra/operator/preconditioner/gauss_seidel.h"
#include "lib_algebra/operator/preconditioner/ilu.h"
#include "lib_algebra/operator/preconditioner/ilut.h"
#include "lib_algebra/operator/preconditioner/vanka.h"
#include "lib_algebra/operator/preconditioner/line_smoothers.h"

// solver
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/operator/linear_solver/cg.h"
#include "lib_algebra/operator/linear_solver/bicgstab.h"
#include "lib_algebra/operator/linear_solver/lu.h"
#ifdef UG_PARALLEL
#include "lib_algebra/operator/linear_solver/feti.h"
	#ifdef UG_HLIBPRO
	#include "lib_algebra/operator/linear_solver/hlibpro.h"
	#endif
#endif

// operator util
#include "lib_algebra/operator/preconditioner/iterator_product.h"
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/vector_writer.h"

#endif /* __H__LIB_ALGEBRA__LIB_ALGEBRA_IMPL__ */
