/*
 * gauss_quad.h
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__
#define __H__LIBDISCRETIZATION__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__

#include "common/common.h"

namespace ug{

/// fixed order gauss quadrature
template <typename TRefElem, int order>
class GaussQuadrature;

} // namespace ug

// include implementation
#include "gauss_quad_edge.h"
#include "gauss_quad_triangle.h"
#include "gauss_quad_quadrilateral.h"
#include "gauss_quad_tetrahedron.h"
#include "gauss_quad_pyramid.h"
#include "gauss_quad_prism.h"
#include "gauss_quad_hexahedron.h"


#endif /* __H__LIBDISCRETIZATION__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__ */
