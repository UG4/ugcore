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

/**
 * \defgroup lib_algebra lib_algebra
 * \brief Algebra Library
 */

#include "common/common.h"

// algebra chooser and types
#include "algebra_type.h"

// common
/**
 * \defgroup lib_algebra_common lib_algebra Common
 * \ingroup lib_algebra
 * \brief Utilities for libAlgebra
 */
#include "common/operations.h"
#include "common/stl_debug.h"

// cpu_algebra
#include "cpu_algebra/algebra_misc.h"
#include "cpu_algebra/vector.h"
#include "cpu_algebra/sparsematrix.h"

// parallel support
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

// small algebra
#include "small_algebra/small_algebra.h"


// operator interfaces
#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/operator_iterator.h"

#endif /* __H__LIB_ALGEBRA__LIB_ALGEBRA__ */
