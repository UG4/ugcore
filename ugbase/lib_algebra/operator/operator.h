/*
 * operator.h
 *
 *  Created on: 01.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR__

///////////////////////////
// Interfaces
///////////////////////////

// operator interface without types
#include "operator_base_interface.h"

// operator interface using template of concrete type
#include "operator_interface.h"

// inverse operators
#include "operator_inverse_interface.h"

// iterator operators (Preconditioners)
#include "operator_iterator_interface.h"

// convergence check
#include "convergence_check.h"

///////////////////////////
// Implementations
///////////////////////////



#include "linear_solver/lu_operator.h"

// solvers
#include "linear_solver/cg_solver.h"
#include "linear_solver/bicgstab_solver.h"
#include "linear_solver/linear_solver.h"

// preconditioner
#include "preconditioner/ilu.h"
#include "preconditioner/ilut.h"
#include "preconditioner/jacobi.h"
#include "preconditioner/gauss_seidel.h"

//#include "preconditioner/amg/amg.h"

#ifdef LAPACK_AVAILABLE
//#include "eigensolver/pinvit.h"
#endif

// feti solver
#include "linear_solver/feti.h"

#endif /* __H__LIB_ALGEBRA__OPERATOR__OPERATOR__ */
