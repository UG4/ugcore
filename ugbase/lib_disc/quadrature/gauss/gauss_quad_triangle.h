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

}; // namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__TRIANGLE__ */

