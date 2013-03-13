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

