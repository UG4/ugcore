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
//  It provides the Gauss Quadratures for a reference quadrilateral.


#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__QUADRILATERAL__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__QUADRILATERAL__

#include "gauss_quad.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Implemented GaussQuadrature for ReferenceQuadrilateral
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceQuadrilateral, 1>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 1>, 2, 1, 1>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 2>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 2>, 2, 2, 4>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 3>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 3>, 2, 3, 4>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 4>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 4>, 2, 4, 6>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 5>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 5>, 2, 5, 7>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 6>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 6>, 2, 6, 10>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 7>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 7>, 2, 7, 12>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 8>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 8>, 2, 8, 16>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 9>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 9>, 2, 9, 17>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 11>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 11>, 2, 11, 24>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 13>
: public GaussQuadBase<GaussQuadrature<ReferenceQuadrilateral, 13>, 2, 13, 33>{};

////////////////////////////////////////////////////////////////////////////////
// template wrapper for orders implemented by higher order
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceQuadrilateral, 12> 
: public GaussQuadrature<ReferenceQuadrilateral, 13>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 10> 
: public GaussQuadrature<ReferenceQuadrilateral, 11>{};

template <>
class GaussQuadrature<ReferenceQuadrilateral, 0> 
: public GaussQuadrature<ReferenceQuadrilateral, 1>{};

}; // namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__QUADRILATERAL__ */

