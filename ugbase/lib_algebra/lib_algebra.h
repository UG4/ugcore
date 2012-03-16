/*
 * lib_algebra.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__LIB_ALGEBRA__
#define __H__LIB_ALGEBRA__LIB_ALGEBRA__

//////////////////////////////////////////////////////////////////////////////
//
// This file is intended to include all parts of lib_algebra
//
//////////////////////////////////////////////////////////////////////////////

#include "common/common.h"

// algebra chooser and types
#include "algebra_type.h"

// common
#include "common/operations.h"
#include "common/stl_debug.h"

// cpu_algebra
#include "cpu_algebra/algebra_misc.h"
#include "cpu_algebra/vector.h"
#include "cpu_algebra/sparsematrix.h"
#include "cpu_algebra/core_smoothers.h"

// parallel support
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

// small algebra
#include "small_algebra/small_algebra.h"

// operator interfaces
#include "lib_algebra/operator/operator_base_interface.h"
#include "lib_algebra/operator/operator_interface.h"
#include "lib_algebra/operator/operator_inverse_interface.h"
#include "lib_algebra/operator/operator_iterator_interface.h"

// preconditioner
#include "lib_algebra/operator/preconditioner/jacobi.h"
#include "lib_algebra/operator/preconditioner/gauss_seidel.h"
#include "lib_algebra/operator/preconditioner/ilu.h"
#include "lib_algebra/operator/preconditioner/ilut.h"

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
#include "lib_algebra/operator/operator_iterator_product.h"
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/vector_writer.h"

#endif /* __H__LIB_ALGEBRA__LIB_ALGEBRA__ */
