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
// This file is intended to include all parts of lib_discretization
//
//////////////////////////////////////////////////////////////////////////////

// algebra chooser and types
#include "algebra_chooser.h"

// common
#include "common/operations.h"
#include "common/stl_debug.h"

// cpu_algebra
#include "cpu_algebra/algebra_misc.h"
#include "cpu_algebra/vector.h"
#include "cpu_algebra/sparsematrix.h"
#include "cpu_algebra/core_smoothers.h"

// operator base interface
#include "operator/operator.h"

// parallel support
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

// small algebra
#include "small_algebra/small_algebra.h"

#endif /* __H__LIB_ALGEBRA__LIB_ALGEBRA__ */
