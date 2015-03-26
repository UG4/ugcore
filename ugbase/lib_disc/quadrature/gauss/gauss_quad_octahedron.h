#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__OCTAHEDRON__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__OCTAHEDRON__

#include "gauss_quad.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Implemented GaussQuadrature for ReferenceOctahedron
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceOctahedron, 2>
: public GaussQuadBase<GaussQuadrature<ReferenceOctahedron, 2>, 3, 2, 16>{};

////////////////////////////////////////////////////////////////////////////////
// template wrapper for orders implemented by higher order
////////////////////////////////////////////////////////////////////////////////
template <>
class GaussQuadrature<ReferenceOctahedron, 1>
: public GaussQuadrature<ReferenceOctahedron, 2>{};

template <>
class GaussQuadrature<ReferenceOctahedron, 0>
: public GaussQuadrature<ReferenceOctahedron, 2>{};

}; // namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__OCTAHEDRON__ */

