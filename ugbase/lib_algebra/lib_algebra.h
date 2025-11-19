/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
//#include "lib_algebra/operator/interface/operator.h"
//#include "lib_algebra/operator/interface/operator_inverse.h"
//#include "lib_algebra/operator/interface/operator_iterator.h"

#endif