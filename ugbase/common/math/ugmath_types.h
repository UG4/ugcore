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
 * \brief typedefs for ugmath
 *
 * \author Sebastian Reiter
 * \date  y08 m11 d12
 *
 * This file contains typedefs for ug math.
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
#include "math_vector_matrix/math_symmetric_matrix.h"
#include "math_vector_matrix/math_tensor.h"

////////////////////////////////////////////////////////////////////////
//	useful typedefs
namespace ug
{

/// \addtogroup vectors
/// \{

/// a 1d vector
typedef MathVector<1, number> vector1;

/// a 2d vector
typedef MathVector<2, number> vector2;

/// a 3d vector
typedef MathVector<3, number> vector3;

/// a 4d vector
typedef MathVector<4, number> vector4;

/// a 2x2 matrix
typedef MathMatrix<2,2, number> matrix22;

/// a 3x3 matrix
typedef MathMatrix<3,3, number> matrix33;

/// a 4x4 matrix
typedef MathMatrix<4,4, number> matrix44;

// end group vectors
/// \}

}

#endif
