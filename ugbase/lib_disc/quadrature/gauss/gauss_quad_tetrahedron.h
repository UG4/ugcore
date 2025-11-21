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
//  It provides the Gauss Quadratures for a reference tetrahedron.


#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__TETRAHEDRON__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__TETRAHEDRON__

#include "gauss_quad.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Implemented GaussQuadrature for ReferenceTetrahedron
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceTetrahedron, 0>
: public GaussQuadBase<GaussQuadrature<ReferenceTetrahedron, 0>, 3, 0, 7>{};

template <>
class GaussQuadrature<ReferenceTetrahedron, 1>
: public GaussQuadBase<GaussQuadrature<ReferenceTetrahedron, 1>, 3, 1, 1>{};

template <>
class GaussQuadrature<ReferenceTetrahedron, 2>
: public GaussQuadBase<GaussQuadrature<ReferenceTetrahedron, 2>, 3, 2, 4>{};

template <>
class GaussQuadrature<ReferenceTetrahedron, 3>
: public GaussQuadBase<GaussQuadrature<ReferenceTetrahedron, 3>, 3, 3, 8>{};

template <>
class GaussQuadrature<ReferenceTetrahedron, 5>
: public GaussQuadBase<GaussQuadrature<ReferenceTetrahedron, 5>, 3, 5, 15>{};

template <>
class GaussQuadrature<ReferenceTetrahedron, 6>
: public GaussQuadBase<GaussQuadrature<ReferenceTetrahedron, 6>, 3, 6, 24>{};

template <>
class GaussQuadrature<ReferenceTetrahedron, 7>
: public GaussQuadBase<GaussQuadrature<ReferenceTetrahedron, 7>, 3, 7, 31>{};

template <>
class GaussQuadrature<ReferenceTetrahedron, 8>
: public GaussQuadBase<GaussQuadrature<ReferenceTetrahedron, 8>, 3, 8, 43>{};

////////////////////////////////////////////////////////////////////////////////
// template wrapper for orders implemented by higher order
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceTetrahedron, 4> 
: public GaussQuadrature<ReferenceTetrahedron, 5>{};

} // namespace ug

#endif