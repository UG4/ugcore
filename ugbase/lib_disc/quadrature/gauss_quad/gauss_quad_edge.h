//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference edge.


#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__EDGE__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__EDGE__

#include "gauss_quad.h"

namespace ug{

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

#endif /* __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__EDGE__ */

