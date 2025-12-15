/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

/**
 * \file ugmath_types.h
 * \brief type definitions for ugmath
 *
 * \author Sebastian Reiter
 * \date  y08 m11 d12
 *
 * This file contains type definitions for ug math.
 */

#ifndef __H__UGMATH_TYPES__
#define __H__UGMATH_TYPES__

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	This header defines common matrix and vector types

////////////////////////////////////////////////////////////////////////
//	include vector and matrix classes
#include "math_vector_matrix/math_vector.h"
#include "math_vector_matrix/math_matrix.h"
// #include "math_vector_matrix/math_symmetric_matrix.h"
//#include "math_vector_matrix/math_tensor.h"

////////////////////////////////////////////////////////////////////////

namespace ug {

/// \addtogroup vectors
/// \{

/// a 1d vector
using vector1 = MathVector<1, number>;

/// a 2d vector
using vector2 = MathVector<2, number>;

/// a 3d vector
using vector3 = MathVector<3, number>;

/// a 4d vector
using vector4 = MathVector<4, number>;

/// a 2x2 matrix
using matrix22 = MathMatrix<2,2, number>;

/// a 3x3 matrix
using matrix33 = MathMatrix<3,3, number>;

/// a 4x4 matrix
using matrix44 = MathMatrix<4,4, number>;

// end group vectors
/// \}

}

#endif
