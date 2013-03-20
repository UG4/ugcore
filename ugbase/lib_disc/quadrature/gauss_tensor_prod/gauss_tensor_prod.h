/*
 * gauss_tensor_prod.h
 *
 *  Created on: 15.03.2013
 *      Author: lisagrau
 */

#ifndef GAUSS_QUADRATURE_HEXAHEDRON_H_
#define GAUSS_QUADRATURE_HEXAHEDRON_H_


#include "../quadrature.h"

namespace ug
{
class GaussQuadratureTriangle : public QuadratureRule<2> {

public:
	//constructor
	GaussQuadratureTriangle(int order);

	//destructor
	~GaussQuadratureTriangle();
};

class GaussQuadratureQuadrilateral : public QuadratureRule<2> {

public:
	//constructor
	GaussQuadratureQuadrilateral(int order);

	//destructor
	~GaussQuadratureQuadrilateral();
};

class GaussQuadratureHexahedron : public QuadratureRule<3> {

public:
	//constructor
	GaussQuadratureHexahedron(int order);

	//destructor
	~GaussQuadratureHexahedron();
};

class GaussQuadratureTetrahedron : public QuadratureRule<3> {

public:
	//constructor
	GaussQuadratureTetrahedron(int order);

	//destructor
	~GaussQuadratureTetrahedron();
};

class GaussQuadraturePrism : public QuadratureRule<3> {

public:
	//constructor
	GaussQuadraturePrism(int order);

	//destructor
	~GaussQuadraturePrism();
};
} // namespace ug

#endif /* GAUSS_QUADRATURE_HEXAHEDRON_H_ */
