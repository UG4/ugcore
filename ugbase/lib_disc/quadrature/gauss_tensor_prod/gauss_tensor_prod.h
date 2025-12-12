/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Lisa Grau, Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_TENSOR_PROD__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_TENSOR_PROD__

/*
 * In this file, quadrature rules for arbitrary order for elements with dim > 1
 * are implemented. The basic idea relies on a transformation of the considered
 * domain to the unit cube [0,1]^d and then using Gauss quadratures to carry out
 * the integrals. Note, that on some elements the transformation introduces
 * jacobian determinants of type (1-x)^{\alpha}(1+x)^{\beta}. In order to
 * integrate those, the gauss-jacobi quadratures are used.
 */

#include "../quadrature.h"

namespace ug {

/**
 * The following class provides QuadratureRules for triangles
 * by using Gauss Jacobi Quadrature and Gauss Legendre Quadrature
 */
class GaussQuadratureTriangle : public QuadratureRule<2> {

	public:
	/// constructor
		explicit GaussQuadratureTriangle(size_t order);

	///	destructor
		~GaussQuadratureTriangle() override;
};

/**
 * The following class provides QuadratureRules for quadrilaterals
 * by using Gauss Legendre Quadrature
 */
class GaussQuadratureQuadrilateral : public QuadratureRule<2> {

	public:
	///	constructor
		explicit GaussQuadratureQuadrilateral(size_t order);

	///	destructor
		~GaussQuadratureQuadrilateral() override;
};

/**
 * The following class provides QuadratureRules for hexahedrons
 * by using Gauss Jacobi Quadrature and Gauss Legendre Quadrature
 */
class GaussQuadratureHexahedron : public QuadratureRule<3> {

	public:
	///	constructor
		explicit GaussQuadratureHexahedron(size_t order);

	///	destructor
		~GaussQuadratureHexahedron() override;
};

/**
 * The following class provides QuadratureRules for tetrahedrons
 * by using Gauss Jacobi Quadrature and Gauss Legendre Quadrature
 */
class GaussQuadratureTetrahedron : public QuadratureRule<3> {

	public:
	///	constructor
		explicit GaussQuadratureTetrahedron(size_t order);

	///	destructor
		~GaussQuadratureTetrahedron() override;
};

/**
 * The following class provides QuadratureRules for prisms
 * by using Gauss Jacobi Quadrature and Gauss Legendre Quadrature
 */
class GaussQuadraturePrism : public QuadratureRule<3> {

	public:
	///	constructor
		explicit GaussQuadraturePrism(size_t order);

	///	destructor
		~GaussQuadraturePrism() override;
};

/**
 * The following class provides QuadratureRules for pyramids
 * by using Gauss Jacobi Quadrature and Gauss Legendre Quadrature
 *
 * The idea of this quadrature is to divide the pyramid into two tetrahdrons given
 * by {x,0,1},{y,0,x},{z,0,1-x} and {x,0,1},{y,x,1},{z,0,1-y} and carry out
 * integration on those using GaussQuadratureTetrehedron
 */
class GaussQuadraturePyramid : public QuadratureRule<3> {

	public:
	///	constructor
		explicit GaussQuadraturePyramid(size_t order);

	///	destructor
		~GaussQuadraturePyramid() override;
};

/**
 * The following class provides QuadratureRules for octahedra
 * by using Gauss Jacobi Quadrature and Gauss Legendre Quadrature
 *
 * The idea of this quadrature is to divide the octahedron into four tetrahdrons
 * and carry out integration on those using GaussQuadratureTetrehedron.
 */
class GaussQuadratureOctahedron : public QuadratureRule<3> {

	public:
	///	constructor
		explicit GaussQuadratureOctahedron(size_t order);

	///	destructor
		~GaussQuadratureOctahedron() override;
};

} // namespace ug

#endif