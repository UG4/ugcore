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

// debug writer
#include "debug_writer.h"
// vector writer
#include "vector_writer.h"

#include "operator_util.h"

// matrix operator methods
#include "matrix_operator_functions.h"


// for compatibilty. remove this later.
#include "operator_impl.h"

#endif /* __H__LIB_ALGEBRA__OPERATOR__OPERATOR__ */
