/*
 * operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR__
#define __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR__

// generic solver routines
#include "operator_util.h"

// special operators
#include "linear_operator/interpolation_operator.h"
#include "linear_operator/projection_operator.h"
#include "linear_operator/prolongation_operator.h"
#include "linear_operator/assembled_linear_operator.h"
#include "linear_operator/multi_grid_solver/mg_solver.h"

// non linear operators
#include "non_linear_operator/assembled_non_linear_operator.h"
#include "non_linear_operator/newton_solver/newton.h"
#include "non_linear_operator/line_search.h"

#endif /* __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR__ */
