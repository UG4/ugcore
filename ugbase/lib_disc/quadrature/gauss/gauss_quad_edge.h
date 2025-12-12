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
//  It provides the Gauss Quadratures for a reference edge.


#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__EDGE__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__EDGE__

#include "gauss_quad.h"
#include "lib_disc/reference_element/reference_element.h"

namespace ug {

////////////////////////////////////////////////////////////////////////////////
// Implemented GaussQuadrature for ReferenceEdge
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceEdge, 1>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 1>, 1, 1, 1>{};

template <>
class GaussQuadrature<ReferenceEdge, 3>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 3>, 1, 3, 2>{};

template <>
class GaussQuadrature<ReferenceEdge, 5>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 5>, 1, 5, 3>{};

template <>
class GaussQuadrature<ReferenceEdge, 7>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 7>, 1, 7, 4>{};

template <>
class GaussQuadrature<ReferenceEdge, 9>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 9>, 1, 9, 5>{};

template <>
class GaussQuadrature<ReferenceEdge, 11>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 11>, 1, 11, 6>{};

template <>
class GaussQuadrature<ReferenceEdge, 13>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 13>, 1, 13, 7>{};

template <>
class GaussQuadrature<ReferenceEdge, 15>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 15>, 1, 15, 8>{};

template <>
class GaussQuadrature<ReferenceEdge, 17>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 17>, 1, 17, 9>{};

template <>
class GaussQuadrature<ReferenceEdge, 19>
: public GaussQuadBase<GaussQuadrature<ReferenceEdge, 19>, 1, 19, 10>{};

////////////////////////////////////////////////////////////////////////////////
// template wrapper for orders implemented by higher order
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceEdge, 18> 
: public GaussQuadrature<ReferenceEdge, 19>{};

template <>
class GaussQuadrature<ReferenceEdge, 16> 
: public GaussQuadrature<ReferenceEdge, 17>{};

template <>
class GaussQuadrature<ReferenceEdge, 14> 
: public GaussQuadrature<ReferenceEdge, 15>{};

template <>
class GaussQuadrature<ReferenceEdge, 12> 
: public GaussQuadrature<ReferenceEdge, 13>{};

template <>
class GaussQuadrature<ReferenceEdge, 10> 
: public GaussQuadrature<ReferenceEdge, 11>{};

template <>
class GaussQuadrature<ReferenceEdge, 8> 
: public GaussQuadrature<ReferenceEdge, 9>{};

template <>
class GaussQuadrature<ReferenceEdge, 6> 
: public GaussQuadrature<ReferenceEdge, 7>{};

template <>
class GaussQuadrature<ReferenceEdge, 4> 
: public GaussQuadrature<ReferenceEdge, 5>{};

template <>
class GaussQuadrature<ReferenceEdge, 2> 
: public GaussQuadrature<ReferenceEdge, 3>{};

template <>
class GaussQuadrature<ReferenceEdge, 0> 
: public GaussQuadrature<ReferenceEdge, 1>{};

}; // namespace ug

#endif