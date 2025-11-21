/*
 * Copyright (c) 2011:  G-CSC, Goethe University Frankfurt
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

//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference triangle.


#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__TRIANGLE__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__TRIANGLE__

#include "gauss_quad.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Implemented GaussQuadrature for ReferenceTriangle
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceTriangle, 1>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 1>, 2, 1, 1>{};

template <>
class GaussQuadrature<ReferenceTriangle, 2>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 2>, 2, 2, 3>{};

template <>
class GaussQuadrature<ReferenceTriangle, 3>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 3>, 2, 3, 4>{};

template <>
class GaussQuadrature<ReferenceTriangle, 4>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 4>, 2, 4, 6>{};

template <>
class GaussQuadrature<ReferenceTriangle, 5>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 5>, 2, 5, 7>{};

template <>
class GaussQuadrature<ReferenceTriangle, 6>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 6>, 2, 6, 12>{};

template <>
class GaussQuadrature<ReferenceTriangle, 7>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 7>, 2, 7, 12>{};

template <>
class GaussQuadrature<ReferenceTriangle, 8>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 8>, 2, 8, 16>{};

template <>
class GaussQuadrature<ReferenceTriangle, 9>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 9>, 2, 9, 19>{};

template <>
class GaussQuadrature<ReferenceTriangle, 10>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 10>, 2, 10, 25>{};

template <>
class GaussQuadrature<ReferenceTriangle, 11>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 11>, 2, 11, 28>{};

template <>
class GaussQuadrature<ReferenceTriangle, 12>
: public GaussQuadBase<GaussQuadrature<ReferenceTriangle, 12>, 2, 12, 33>{};

////////////////////////////////////////////////////////////////////////////////
// template wrapper for orders implemented by higher order
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceTriangle, 0> 
: public GaussQuadrature<ReferenceTriangle, 1>{};

} // namespace ug

#endif
