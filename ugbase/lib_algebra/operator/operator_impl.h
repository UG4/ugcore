/*
 * operator_impl.h
 *
 *  Created on: 23.03.2011
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_IMPL__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_IMPL__


///////////////////////////
// Implementations
///////////////////////////


// products and sums of iterator operators
#include "operator_iterator_product.h"


#include "linear_solver/lu.h"

// solvers
#include "linear_solver/cg.h"
#include "linear_solver/bicgstab.h"
#include "linear_solver/linear_solver.h"

// preconditioner
#include "preconditioner/ilu.h"
#include "preconditioner/ilut.h"
#include "preconditioner/jacobi.h"
#include "preconditioner/gauss_seidel.h"


#ifdef LAPACK_AVAILABLE
//#include "eigensolver/pinvit.h"
#endif

// feti solver
#include "linear_solver/dirichletdirichlet.h"
#include "linear_solver/feti.h"

// HLIBpro based solver
#ifdef USE_HLIBPRO
#include "linear_solver/hlib_operator.h"
#endif


#endif /* __H__LIB_ALGEBRA__OPERATOR__OPERATOR_IMPL__ */
