//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference pyramid.


#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__PYRAMID__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__PYRAMID__

#include "gauss_quad.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Implemented GaussQuadrature for ReferencePyramid
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferencePyramid, 2>
: public GaussQuadBase<GaussQuadrature<ReferencePyramid, 2>, 3, 2, 8>{};

////////////////////////////////////////////////////////////////////////////////
// template wrapper for orders implemented by higher order
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferencePyramid, 1> 
: public GaussQuadrature<ReferencePyramid, 2>{};

template <>
class GaussQuadrature<ReferencePyramid, 0> 
: public GaussQuadrature<ReferencePyramid, 2>{};

}; // namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__PYRAMID__ */

