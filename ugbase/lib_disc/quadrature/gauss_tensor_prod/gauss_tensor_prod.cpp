/*
 * gauss_tensor_prod.cpp
 *
 *  Created on: 15.03.2013
 *      Author: lisagrau
 */


#include "../quadrature.h"
#include "common/util/provider.h"
#include "gauss_tensor_prod.h"
#include "../gauss_legendre/gauss_legendre.h"
#include "../gauss_jacobi/gauss_jacobi10.h"
#include "../gauss_jacobi/gauss_jacobi20.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"

namespace ug {

GaussQuadratureTriangle::GaussQuadratureTriangle(size_t order)
{
	GaussLegendre quadRule = GaussLegendre(order);
	GaussJacobi10 quadRule10 = GaussJacobi10(order);

	m_order = std::min(quadRule.order(), quadRule10.order());
	m_numPoints = quadRule10.size() * quadRule.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	size_t cnt = 0;
	for(size_t i = 0; i < quadRule.size(); i++){
		for(size_t j = 0; j < quadRule10.size(); j++, cnt++){
			pvPoint[cnt][0] = quadRule10.point(j)[0];
			pvPoint[cnt][1] = (1.0 - quadRule10.point(j)[0] ) * quadRule.point(i)[0];
			pvWeight[cnt] = quadRule.weight(i) * quadRule10.weight(j);
		}
	}
};

GaussQuadratureTriangle::~GaussQuadratureTriangle()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

GaussQuadratureQuadrilateral::GaussQuadratureQuadrilateral(size_t order)
{
	GaussLegendre quadRule = GaussLegendre(order);

	m_order = std::min(quadRule.order(), quadRule.order());
	m_numPoints = quadRule.size() * quadRule.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	size_t cnt  = 0;
	for(size_t i = 0; i < quadRule.size(); i ++) {
		for(size_t j = 0; j < quadRule.size(); j++, cnt++) {
			pvPoint[cnt][0] = quadRule.point(i)[0];
			pvPoint[cnt][1] = quadRule.point(j)[0];
			pvWeight[cnt] = quadRule.weight(i) * quadRule.weight(j);
		}
	}

};

GaussQuadratureQuadrilateral::~GaussQuadratureQuadrilateral()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

GaussQuadratureHexahedron::GaussQuadratureHexahedron(size_t order)
{
	GaussLegendre quadRule = GaussLegendre(order);

	m_order = std::min(quadRule.order(), std::min(quadRule.order(), quadRule.order()));
	m_numPoints = quadRule.size() * quadRule.size() * quadRule.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	size_t cnt  = 0;
	for(size_t i = 0; i < quadRule.size(); i ++) {
		for(size_t j = 0; j < quadRule.size(); j++) {
			for(size_t k = 0; k < quadRule.size(); k++, cnt++) {
				pvPoint[cnt][0] = quadRule.point(i)[0];
				pvPoint[cnt][1] = quadRule.point(j)[0];
				pvPoint[cnt][2] = quadRule.point(k)[0];
				pvWeight[cnt] = quadRule.weight(i) * quadRule.weight(j) * quadRule.weight(k);
			}
		}
	}

};

GaussQuadratureHexahedron::~GaussQuadratureHexahedron()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

GaussQuadratureTetrahedron::GaussQuadratureTetrahedron(size_t order)
{
	GaussLegendre quadRule = GaussLegendre(order);
	GaussJacobi10 quadRule10 = GaussJacobi10(order);
	GaussJacobi20 quadRule20 = GaussJacobi20(order);

	m_order = std::min(quadRule.order(), std::min(quadRule10.order(), quadRule20.order()));
	m_numPoints = quadRule.size() * quadRule10.size() * quadRule20.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	size_t cnt = 0;
	for(size_t i = 0; i < quadRule20.size(); i++) {
		for(size_t j = 0; j < quadRule10.size(); j++) {
			for(size_t k = 0; k < quadRule.size(); k++, cnt++) {
				pvPoint[cnt][0] = quadRule20.point(i)[0];
				pvPoint[cnt][1] = (1.0 - quadRule20.point(i)[0] ) * quadRule10.point(j)[0];
				pvPoint[cnt][2] = (1.0 - quadRule20.point(i)[0]) * (1.0 - quadRule10.point(j)[0]) * quadRule.point(k)[0];
				pvWeight[cnt] = quadRule20.weight(i) * quadRule10.weight(j) * quadRule.weight(k);
			}
		}
	}

};

GaussQuadratureTetrahedron::~GaussQuadratureTetrahedron()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

GaussQuadraturePrism::GaussQuadraturePrism(size_t order)
{
	GaussLegendre quadRule = GaussLegendre(order);
	GaussJacobi10 quadRule10 = GaussJacobi10(order);

	m_order = std::min(quadRule.order(), quadRule10.order());
	m_numPoints = quadRule10.size() * quadRule.size() * quadRule.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	size_t cnt = 0;
	for(size_t i = 0; i < quadRule10.size(); i++) {
		for(size_t j = 0; j < quadRule.size(); j++) {
			for(size_t k = 0; k < quadRule.size(); k++, cnt++) {
				pvPoint[cnt][0] = quadRule10.point(i)[0];
				pvPoint[cnt][1] = (1.0 - quadRule10.point(i)[0]) * quadRule.point(j)[0];
				pvPoint[cnt][2] = quadRule.point(k)[0];
				pvWeight[cnt] = quadRule10.weight(i) * quadRule.weight(j) * quadRule.weight(k);
			}
		}
	}

};

GaussQuadraturePrism::~GaussQuadraturePrism()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

GaussQuadraturePyramid::GaussQuadraturePyramid(size_t order)
{
	GaussQuadratureTetrahedron quadRule = GaussQuadratureTetrahedron(order);

	m_order = quadRule.order();
	m_numPoints = quadRule.size() * 2;
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	MathVector<3> Tet1Co[4];
	Tet1Co[0] = MathVector<3>(0,0,0);
	Tet1Co[1] = MathVector<3>(1,1,0);
	Tet1Co[2] = MathVector<3>(0,0,1);
	Tet1Co[3] = MathVector<3>(0,1,0);

	MathVector<3> Tet2Co[4];
	Tet2Co[0] = MathVector<3>(0,0,0);
	Tet2Co[1] = MathVector<3>(1,0,0);
	Tet2Co[2] = MathVector<3>(1,1,0);
	Tet2Co[3] = MathVector<3>(0,0,1);

	DimReferenceMapping<3, 3>& map1 =
			ReferenceMappingProvider::get<3,3>(ROID_TETRAHEDRON, Tet1Co);

	size_t cnt = 0;
	for(size_t i = 0; i < quadRule.size(); i++, cnt++) {
		map1.local_to_global(pvPoint[cnt], quadRule.point(i));
		pvWeight[cnt] = quadRule.weight(i) * map1.sqrt_gram_det(quadRule.point(i));
	}


	DimReferenceMapping<3, 3>& map2 =
			ReferenceMappingProvider::get<3,3>(ROID_TETRAHEDRON, Tet2Co);

	for(size_t j = 0; j < quadRule.size(); j++, cnt++) {
		map2.local_to_global(pvPoint[cnt], quadRule.point(j));
		pvWeight[cnt] = quadRule.weight(j) * map2.sqrt_gram_det(quadRule.point(j));
	}
};

GaussQuadraturePyramid::~GaussQuadraturePyramid()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}
}

