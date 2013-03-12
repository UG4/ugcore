/*
 * gauss_quad.h
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__

#include "common/common.h"
#include "../quadrature.h"
#include "../quadrature_provider.h"

namespace ug{

/// fixed order gauss quadrature
template <typename TRefElem, int order>
class GaussQuadrature;

/// flexible order gauss quadrature
/**
 * Providing gauss quadrature for an reference element
 * \tparam 		TRefElem		Reference Element Type
 */
template <typename TRefElem>
class FlexGaussQuadrature
	: public QuadratureRule<TRefElem::dim>
{
	public:
	///	Constructor
		FlexGaussQuadrature(int order);

	///	Destructor
		~FlexGaussQuadrature() {}
};

// registering function
template <typename TRefElem>
bool RegisterGaussQuadRule(QuadratureRuleProvider<TRefElem::dim>& factory);
template <> bool RegisterGaussQuadRule<ReferenceVertex>(QuadratureRuleProvider<ReferenceVertex::dim>& factory);
template <> bool RegisterGaussQuadRule<ReferenceEdge>(QuadratureRuleProvider<ReferenceEdge::dim>& factory);
template <> bool RegisterGaussQuadRule<ReferenceTriangle>(QuadratureRuleProvider<ReferenceTriangle::dim>& factory);
template <> bool RegisterGaussQuadRule<ReferenceQuadrilateral>(QuadratureRuleProvider<ReferenceQuadrilateral::dim>& factory);
template <> bool RegisterGaussQuadRule<ReferenceTetrahedron>(QuadratureRuleProvider<ReferenceTetrahedron::dim>& factory);
template <> bool RegisterGaussQuadRule<ReferencePrism>(QuadratureRuleProvider<ReferencePrism::dim>& factory);
template <> bool RegisterGaussQuadRule<ReferencePyramid>(QuadratureRuleProvider<ReferencePyramid::dim>& factory);
template <> bool RegisterGaussQuadRule<ReferenceHexahedron>(QuadratureRuleProvider<ReferenceHexahedron::dim>& factory);

// registering function
template <int dim>
inline bool RegisterGaussQuadRuleDim(QuadratureRuleProvider<dim>& factory);

// implementation 0d
template <>
inline bool RegisterGaussQuadRuleDim(QuadratureRuleProvider<0>& factory)
{
	bool bRet = true;
	bRet &= RegisterGaussQuadRule<ReferenceVertex>(factory);
	return bRet;
}
// implementation 1d
template <>
inline bool RegisterGaussQuadRuleDim(QuadratureRuleProvider<1>& factory)
{
	bool bRet = true;
	bRet &= RegisterGaussQuadRule<ReferenceEdge>(factory);
	return bRet;
}
// implementation 2d
template <>
inline bool RegisterGaussQuadRuleDim(QuadratureRuleProvider<2>& factory)
{
	bool bRet = true;
	bRet &= RegisterGaussQuadRule<ReferenceTriangle>(factory);
	bRet &= RegisterGaussQuadRule<ReferenceQuadrilateral>(factory);
	return bRet;
}
// implementation 3d
template <>
inline bool RegisterGaussQuadRuleDim(QuadratureRuleProvider<3>& factory)
{
	bool bRet = true;
	bRet &= RegisterGaussQuadRule<ReferenceTetrahedron>(factory);
	bRet &= RegisterGaussQuadRule<ReferencePrism>(factory);
	bRet &= RegisterGaussQuadRule<ReferencePyramid>(factory);
	bRet &= RegisterGaussQuadRule<ReferenceHexahedron>(factory);
	return bRet;
}

} // namespace ug

// include implementation
#include "gauss_quad_edge.h"
#include "gauss_quad_triangle.h"
#include "gauss_quad_quadrilateral.h"
#include "gauss_quad_tetrahedron.h"
#include "gauss_quad_pyramid.h"
#include "gauss_quad_prism.h"
#include "gauss_quad_hexahedron.h"


#endif /* __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__ */
