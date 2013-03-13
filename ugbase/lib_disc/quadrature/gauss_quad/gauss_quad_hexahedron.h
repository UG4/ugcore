//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference hexahedron.


#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__HEXAHEDRON__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__HEXAHEDRON__

#include "gauss_quad.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Implemented GaussQuadrature for ReferenceHexahedron
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceHexahedron, 2>
: public GaussQuadBase<GaussQuadrature<ReferenceHexahedron, 2>, 3, 2, 8>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 3>
: public GaussQuadBase<GaussQuadrature<ReferenceHexahedron, 3>, 3, 3, 6>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 5>
: public GaussQuadBase<GaussQuadrature<ReferenceHexahedron, 5>, 3, 5, 14>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 7>
: public GaussQuadBase<GaussQuadrature<ReferenceHexahedron, 7>, 3, 7, 31>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 8>
: public GaussQuadBase<GaussQuadrature<ReferenceHexahedron, 8>, 3, 8, 47>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 9>
: public GaussQuadBase<GaussQuadrature<ReferenceHexahedron, 9>, 3, 9, 58>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 11>
: public GaussQuadBase<GaussQuadrature<ReferenceHexahedron, 11>, 3, 11, 90>{};

////////////////////////////////////////////////////////////////////////////////
// template wrapper for orders implemented by higher order
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceHexahedron, 10> 
: public GaussQuadrature<ReferenceHexahedron, 11>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 6> 
: public GaussQuadrature<ReferenceHexahedron, 7>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 4> 
: public GaussQuadrature<ReferenceHexahedron, 5>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 1> 
: public GaussQuadrature<ReferenceHexahedron, 2>{};

template <>
class GaussQuadrature<ReferenceHexahedron, 0> 
: public GaussQuadrature<ReferenceHexahedron, 1>{};

}; // namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__HEXAHEDRON__ */

