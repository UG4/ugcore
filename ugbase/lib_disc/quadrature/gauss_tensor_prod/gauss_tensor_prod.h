/*
 * gauss_tensor_prod.h
 *
 *  Created on: 15.03.2013
 *      Author: lisagrau, andreasvogel
 */

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_TENSOR_PROD__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_TENSOR_PROD__

/*
 * In this file, quadrature rules for arbitray order for elements with dim > 1
 * are implemented. The basic idea relies on a transformation of the considered
 * domain to the unit cube [0,1]^d and then using Gauss quadratures to carry out
 * the integrals. Note, that on some elements the transformation introduces
 * jacobian determinants of type (1-x)^{\alpha}(1+x)^{\beta}. In order to
 * integrate those, the gauss-jacobi quadratures are used.
 */

#include "../quadrature.h"

namespace ug
{

/**
 * The following class provides QuadratureRules for triangles
 * by using Gauss Jacobi Quadrature and Gauss Legendre Qadrature
 */
class GaussQuadratureTriangle : public QuadratureRule<2> {

	public:
	/// constructor
		GaussQuadratureTriangle(size_t order);

	///	destructor
		~GaussQuadratureTriangle();
};

/**
 * The following class provides QuadratureRules for quadrilaterals
 * by using Gauss Legendre Qadrature
 */
class GaussQuadratureQuadrilateral : public QuadratureRule<2> {

	public:
	///	constructor
		GaussQuadratureQuadrilateral(size_t order);

	///	destructor
		~GaussQuadratureQuadrilateral();
};

/**
 * The following class provides QuadratureRules for hexahedrons
 * by using Gauss Jacobi Quadrature and Gauss Legendre Qadrature
 */
class GaussQuadratureHexahedron : public QuadratureRule<3> {

	public:
	///	constructor
		GaussQuadratureHexahedron(size_t order);

	///	destructor
		~GaussQuadratureHexahedron();
};

/**
 * The following class provides QuadratureRules for tetrahedrons
 * by using Gauss Jacobi Quadrature and Gauss Legendre Qadrature
 */
class GaussQuadratureTetrahedron : public QuadratureRule<3> {

	public:
	///	constructor
		GaussQuadratureTetrahedron(size_t order);

	///	destructor
		~GaussQuadratureTetrahedron();
};

/**
 * The following class provides QuadratureRules for prisms
 * by using Gauss Jacobi Quadrature and Gauss Legendre Qadrature
 */
class GaussQuadraturePrism : public QuadratureRule<3> {

	public:
	///	constructor
		GaussQuadraturePrism(size_t order);

	///	destructor
		~GaussQuadraturePrism();
};

/**
 * The following class provides QuadratureRules for pyramids
 * by using Gauss Jacobi Quadrature and Gauss Legendre Qadrature
 *
 * The idea of this quadrature is to divide the pyramid into two tetrahdrons given
 * by {x,0,1},{y,0,x},{z,0,1-x} and {x,0,1},{y,x,1},{z,0,1-y} and carry out
 * integration on those using GaussQuadratureTetrehedron
 */
class GaussQuadraturePyramid : public QuadratureRule<3> {

	public:
	///	constructor
		GaussQuadraturePyramid(size_t order);

	///	destructor
		~GaussQuadraturePyramid();
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__GAUSS_TENSOR_PROD__ */
