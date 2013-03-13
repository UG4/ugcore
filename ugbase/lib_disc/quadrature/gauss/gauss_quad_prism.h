//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference prism.


#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__PRISM__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__PRISM__

#include "gauss_quad.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Implemented GaussQuadrature for ReferencePrism
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferencePrism, 0>
: public GaussQuadBase<GaussQuadrature<ReferencePrism, 0>, 3, 0, 6>{};

template <>
class GaussQuadrature<ReferencePrism, 2>
: public GaussQuadBase<GaussQuadrature<ReferencePrism, 2>, 3, 2, 6>{};

////////////////////////////////////////////////////////////////////////////////
// template wrapper for orders implemented by higher order
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferencePrism, 1> 
: public GaussQuadrature<ReferencePrism, 2>{};

}; // namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__PRISM__ */

