/*
 * gauss_tensor_prod.h
 *
 *  Created on: 15.03.2013
 *      Author: lisagrau
 */

#ifndef GAUSS_QUADRATURE_HEXAHEDRON_H_
#define GAUSS_QUADRATURE_HEXAHEDRON_H_

/**
 * if higher orders for the following classes are needed, the mathematica script
 * will generate them, further information concerning the quadrature rules can be
 * found on wikipedia
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
		//constructor
		GaussQuadratureTriangle(int order);

		//destructor
		~GaussQuadratureTriangle();
	};

	/**
	 * The following class provides QuadratureRules for quadrilaterals
	 * by using Gauss Legendre Qadrature
	 */

	class GaussQuadratureQuadrilateral : public QuadratureRule<2> {

	public:
		//constructor
		GaussQuadratureQuadrilateral(int order);

		//destructor
		~GaussQuadratureQuadrilateral();
	};

	/**
	 * The following class provides QuadratureRules for hexahedrons
	 * by using Gauss Jacobi Quadrature and Gauss Legendre Qadrature
	 */

	class GaussQuadratureHexahedron : public QuadratureRule<3> {

	public:
		//constructor
		GaussQuadratureHexahedron(int order);

		//destructor
		~GaussQuadratureHexahedron();
	};

	/**
	 * The following class provides QuadratureRules for tetrahedrons
	 * by using Gauss Jacobi Quadrature and Gauss Legendre Qadrature
	 */

	class GaussQuadratureTetrahedron : public QuadratureRule<3> {

	public:
		//constructor
		GaussQuadratureTetrahedron(int order);

		//destructor
		~GaussQuadratureTetrahedron();
	};

	/**
	 * The following class provides QuadratureRules for prisms
	 * by using Gauss Jacobi Quadrature and Gauss Legendre Qadrature
	 */

	class GaussQuadraturePrism : public QuadratureRule<3> {

	public:
		//constructor
		GaussQuadraturePrism(int order);

		//destructor
		~GaussQuadraturePrism();
	};

	/**
	 * The following class provides QuadratureRules for pyramids
	 * by using Gauss Jacobi Quadrature and Gauss Legendre Qadrature
	 */

	class GaussQuadraturePyramid : public QuadratureRule<3> {

	public:
		//constructor
		GaussQuadraturePyramid(int order);

		//destructor
		~GaussQuadraturePyramid();
	};

} // namespace ug

#endif /* GAUSS_QUADRATURE_HEXAHEDRON_H_ */
