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

}; // namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__TETRAHEDRON__ */

