/*
 * operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR__
#define __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR__

// operator interface
#include "operator_interface.h"

// generic solver routines
#include "operator_util.h"

// special operators
#include "linear_operator/interpolation_operator.h"
#include "linear_operator/projection_operator.h"
#include "linear_operator/transfer_operator.h"
#include "linear_operator/assembled_linear_operator.h"
#include "linear_operator/multi_grid_solver/mg_solver.h"
#include "linear_operator/linear_solver.h"
#include "linear_operator/cg_solver.h"

#ifdef USE_MARTIN_ALGEBRA
#ifdef LAPACK_AVAILABLE
#ifdef BLAS_AVAILABLE
#include "linear_operator/lapack_lu_operator.h"
#endif
#endif
#endif

// non linear operators
#include "non_linear_operator/assembled_non_linear_operator.h"
#include "non_linear_operator/newton_solver/newton.h"

#endif /* __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR__ */
